#!/usr/bin/env python3
"""
C2H2 ULTRA-BOOST OPTIMIZATION (Version 2)

USER INSIGHT: In CH-suppression solution, C2H2 is too low (3.4e11, 0.11×)
and 97% is lost to walls. Need to push C2H2 to ~1e13 to feed C2 production!

Strategy:
1. TARGET C2H2 at 1e13 (3× higher than baseline 3e12, 30× higher than current 3.4e11)
2. REDUCE C2H2 wall losses (stick_C2H2, loss_C2H2) as tunable parameters
3. BOOST H (needed for C2H2+H→C2 pathway)
4. MAINTAIN CH suppression (keep high CH wall sticking from previous success)
5. Use H_drift_gain to decouple H from chemistry

New Weights:
- H: 30 (boost for C2H2→C2 pathway)
- CH: 80 (MAX - maintain suppression)
- C2: 60 (primary target)
- C2H2: 40 (VERY HIGH - push to 1e13)

Goal: H, CH, C2 all within 0.6-1.4×, C2H2 at ~1e13, charge balance <100%
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import time
import signal
from contextlib import contextmanager
import os
import json

from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database, get_tunable_rates_for_target
from build_reactions import build_reactions

# Experimental targets - C2H2 BOOSTED to 1e13!
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
    'C2H2': 1.0e13,  # ← BOOSTED! (was 3e12, now 1e13 to feed C2)
}

# Create results directory
os.makedirs('optimization_results_C2H2_ultra_boost_v2', exist_ok=True)

# Global counter for iterations
best_result = {'objective': 1e10, 'params': None, 'densities': None}


class TimeoutException(Exception):
    pass


@contextmanager
def time_limit(seconds):
    """Context manager to timeout function calls."""
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


# Modified PlasmaODE that reads H_drift_gain from params
class PlasmaODE_H_Tunable:
    """ODE function with tunable H sheath source."""

    def __init__(self, params):
        self.species = params['species']
        self.ns = len(self.species)
        self.R = params['R']
        self.nr = len(self.R)
        self.k = params['k']
        self.tags = params['tags']

        # Read H_drift_gain from params (TUNABLE!)
        self.H_drift_gain = params.get('H_drift_gain', 3.2e17)

        # Cache species indices
        self.e_idx = self.species.index('e')
        self.Ar_idx = self.species.index('Ar')
        self.CH4_idx = self.species.index('CH4')
        self.H_idx = self.species.index('H')

        # Build stoichiometry arrays
        self._build_stoich()

    def _build_stoich(self):
        """Build stoichiometry arrays."""
        from scipy import sparse

        rows_prod, cols_prod, vals_prod = [], [], []
        rows_react, cols_react, vals_react = [], [], []

        for rxn_idx, reaction in enumerate(self.R):
            prod_species = np.where(reaction.products > 0)[0]
            for sp_idx in prod_species:
                rows_prod.append(sp_idx)
                cols_prod.append(rxn_idx)
                vals_prod.append(reaction.products[sp_idx])

            react_species = np.where(reaction.reactants > 0)[0]
            for sp_idx in react_species:
                rows_react.append(sp_idx)
                cols_react.append(rxn_idx)
                vals_react.append(reaction.reactants[sp_idx])

        self.prod_matrix = sparse.csr_matrix(
            (vals_prod, (rows_prod, cols_prod)),
            shape=(self.ns, self.nr)
        )
        self.react_matrix = sparse.csr_matrix(
            (vals_react, (rows_react, cols_react)),
            shape=(self.ns, self.nr)
        )

    def __call__(self, t, y):
        """ODE function dy/dt = f(t, y)"""
        y = np.maximum(y, 1e-10)

        rates = np.zeros(self.nr)
        for rxn_idx, reaction in enumerate(self.R):
            rate_constant = self.k.get(reaction.rate, 0.0)

            if 'drift' in self.tags[rxn_idx]:
                ne = y[self.e_idx]
                rates[rxn_idx] = rate_constant * ne

            elif 'loss' in self.tags[rxn_idx] or 'stick' in self.tags[rxn_idx]:
                react_species = np.where(reaction.reactants > 0)[0]
                if len(react_species) > 0:
                    rates[rxn_idx] = rate_constant * y[react_species[0]]

            else:
                rate = rate_constant
                react_species = np.where(reaction.reactants > 0)[0]
                for sp_idx in react_species:
                    rate *= y[sp_idx] ** reaction.reactants[sp_idx]
                rates[rxn_idx] = rate

        production = self.prod_matrix.dot(rates)
        consumption = self.react_matrix.dot(rates)
        dydt = production - consumption

        # Add H sheath source (TUNABLE!)
        dydt[self.H_idx] += self.H_drift_gain

        return dydt


def analyze_chemistry(y, species, params):
    """Analyze chemical fluxes for key species."""
    analysis = {}

    for target_sp in ['H', 'CH', 'C2', 'C2H2']:
        if target_sp not in species:
            continue

        target_idx = species.index(target_sp)
        prod_fluxes = []
        loss_fluxes = []

        for rxn_idx, reaction in enumerate(params['R']):
            rate_constant = params['k'].get(reaction.rate, 0.0)

            if 'drift' in params['tags'][rxn_idx]:
                rate = rate_constant * y[species.index('e')]
            elif 'loss' in params['tags'][rxn_idx] or 'stick' in params['tags'][rxn_idx]:
                react_species = np.where(reaction.reactants > 0)[0]
                if len(react_species) > 0:
                    rate = rate_constant * y[react_species[0]]
                else:
                    rate = 0.0
            else:
                rate = rate_constant
                react_species = np.where(reaction.reactants > 0)[0]
                for sp_idx in react_species:
                    rate *= y[sp_idx] ** reaction.reactants[sp_idx]

            if reaction.products[target_idx] > 0:
                prod_fluxes.append((reaction.rate, rate * reaction.products[target_idx]))
            if reaction.reactants[target_idx] > 0:
                loss_fluxes.append((reaction.rate, rate * reaction.reactants[target_idx]))

        prod_fluxes.sort(key=lambda x: x[1], reverse=True)
        loss_fluxes.sort(key=lambda x: x[1], reverse=True)

        analysis[target_sp] = {
            'density': float(y[target_idx]),
            'production': sum(f[1] for f in prod_fluxes),
            'loss': sum(f[1] for f in loss_fluxes),
            'top_production': prod_fluxes[:10],
            'top_loss': loss_fluxes[:10]
        }

    return analysis


def select_tunable_rates():
    """Select rates to tune, focusing on C2H2 wall losses."""
    # Start with existing tunable rates
    base_rates = get_tunable_rates_for_target('C2')

    # base_rates is a dict, so get the keys
    if isinstance(base_rates, dict):
        base_rate_names = list(base_rates.keys())
    else:
        base_rate_names = list(base_rates)

    # Add critical C2H2 loss rates (reduce these to boost C2H2!)
    c2h2_rates = [
        'loss_C2H2_11_19',      # ← KEY! C2H2 gas-phase loss
        'stick_C2H2_9_11',      # ← KEY! C2H2 wall sticking
        'loss_C2H2Star_11_25',  # C2H2* loss
    ]

    # Add CH suppression rates (keep these high!)
    ch_rates = [
        'stick_CH_9_3',         # ← Keep HIGH for CH suppression
        'loss_CH_11_9',         # ← Keep HIGH for CH suppression
    ]

    # Add C2H2→C2 conversion pathway
    c2_production = [
        'C2H2_H_C2_H2_H_cm3_7_50',  # ← KEY! C2H2 + H → C2 + H2 + H
        'e_C2H2_C2_H2_cm3_1_16',     # e + C2H2 → C2 + H2
    ]

    # Add C2 loss rates
    c2_rates = [
        'loss_C2_11_3',
        'stick_C2_9_9',
        'C2_H_CH_C_cm3_7_6',    # C2 + H → CH + C (suppress this!)
    ]

    # Combine all
    all_rates = list(set(base_rate_names + c2h2_rates + ch_rates + c2_production + c2_rates))

    return all_rates


def run_simulation_with_logging(rate_values, E_field, ne, Te, H_drift_gain, params_base, log_file=None):
    """Run plasma simulation with given parameters."""
    try:
        params = params_base.copy()
        params['E_field'] = E_field
        params['ne'] = ne
        params['Te'] = Te
        params['H_drift_gain'] = H_drift_gain

        k = define_rates_tunable(params)
        db = get_complete_rate_database()

        for name, val in rate_values.items():
            if name in k:
                if name in db:
                    val = np.clip(val, db[name].min, db[name].max)
                k[name] = val

        for name, rate_db in db.items():
            if name in k:
                if k[name] < rate_db.min:
                    k[name] = rate_db.min
                elif k[name] > rate_db.max:
                    k[name] = rate_db.max

        params['k'] = k
        params['R'], params['tags'] = build_reactions(params)

        species = params['species']
        ns = len(species)
        y0 = np.ones(ns) * 1e3

        def set_density(name, value):
            try:
                idx = species.index(name)
                y0[idx] = value
            except ValueError:
                pass

        set_density('e', ne)
        set_density('Ar', 0.85 * 9.66e15)
        set_density('CH4', 0.15 * 9.66e15)
        set_density('ArPlus', 1e7)
        set_density('CH4Plus', 1e5)
        set_density('CH3Plus', 1e5)
        set_density('CH5Plus', 1e3)
        set_density('ArHPlus', 5e5)
        set_density('CH3Minus', 5e4)
        set_density('H2', 1e12)
        set_density('ArStar', 5e6)
        set_density('H', 1e11)
        set_density('C2', 5e7)
        set_density('CH', 5e4)
        set_density('C2H4', 5e7)
        set_density('C2H6', 1e6)
        set_density('CH2', 1e11)
        set_density('C2H2', 1e12)
        set_density('C2H5', 1e6)
        set_density('CH3', 5e7)
        set_density('C', 5e7)

        ode_func = PlasmaODE_H_Tunable(params)

        try:
            with time_limit(30):
                sol = solve_ivp(
                    ode_func,
                    (0, 100),
                    y0,
                    method='BDF',
                    rtol=1e-5,
                    atol=1e-6,
                    max_step=10.0
                )
        except TimeoutException:
            return None

        if not sol.success:
            return None

        y_final = sol.y[:, -1]

        def get_density(name):
            try:
                idx = species.index(name)
                return y_final[idx]
            except ValueError:
                return 0.0

        results = {
            'H': get_density('H'),
            'CH': get_density('CH'),
            'C2': get_density('C2'),
            'C2H2': get_density('C2H2')
        }

        # Calculate ion densities for charge balance
        positive_ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                        'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus',
                        'H3Plus', 'CHPlus', 'H2Plus']

        n_i_total = sum(get_density(ion) for ion in positive_ions)
        results['n_i_total'] = n_i_total
        results['charge_imbalance'] = abs(n_i_total - ne) / ne

        # Calculate C2H2+H→C2 flux
        try:
            c2h2_idx = species.index('C2H2')
            h_idx = species.index('H')
            rate_key = 'C2H2_H_C2_H2_H_cm3_7_50'
            if rate_key in k:
                flux = k[rate_key] * y_final[c2h2_idx] * y_final[h_idx]
                results['C2H2_H_C2_flux'] = flux
            else:
                results['C2H2_H_C2_flux'] = 0.0
        except:
            results['C2H2_H_C2_flux'] = 0.0

        # If logging requested, save detailed analysis
        if log_file:
            all_densities = {species[i]: float(y_final[i]) for i in range(ns)}
            chemistry = analyze_chemistry(y_final, species, params)

            log_data = {
                'Te': Te,
                'Ne': ne,
                'E_field': E_field,
                'H_drift_gain': H_drift_gain,
                'rate_values': {k: float(v) for k, v in rate_values.items()},
                'target_densities': {k: results[k] for k in ['H', 'CH', 'C2', 'C2H2']},
                'charge_balance': {
                    'n_i_total': float(n_i_total),
                    'n_e': float(ne),
                    'imbalance_percent': float(results['charge_imbalance'] * 100)
                },
                'C2H2_H_C2_flux': float(results['C2H2_H_C2_flux']),
                'all_densities': all_densities,
                'chemistry_analysis': chemistry
            }

            with open(log_file, 'w') as f:
                json.dump(log_data, f, indent=2)

        return results

    except Exception as e:
        return None


def objective_function(x, param_names, params_base):
    """Objective function: weighted error from targets."""
    global best_result

    if not hasattr(objective_function, 'counter'):
        objective_function.counter = 0
        objective_function.start_time = time.time()

    objective_function.counter += 1

    if objective_function.counter % 10 == 0:
        elapsed = time.time() - objective_function.start_time
        print(f"  [{objective_function.counter} evaluations, {elapsed/60:.1f} min elapsed]")

    # Extract parameters
    rate_values = {name: val for name, val in zip(param_names[:-4], x[:-4])}
    E_field = x[-4]
    Te = x[-3]
    ne = x[-2]
    H_drift_gain = x[-1]

    params_updated = params_base.copy()
    params_updated['ne'] = ne

    results = run_simulation_with_logging(rate_values, E_field, ne, Te, H_drift_gain, params_updated)

    if results is None:
        return 1e10

    # C2H2-BOOST WEIGHTS
    weights = {
        'H': 30.0,    # Boost H for C2H2→C2 pathway
        'CH': 80.0,   # MAXIMUM - maintain CH suppression
        'C2': 60.0,   # Primary target
        'C2H2': 40.0  # VERY HIGH - push C2H2 to 1e13!
    }

    species_error = 0.0
    for species in ['H', 'CH', 'C2', 'C2H2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        species_error += weights[species] * rel_error ** 2

    # Charge balance penalty
    charge_penalty = 20.0 * results['charge_imbalance']**2

    total_error = species_error + charge_penalty

    # Save best result
    if total_error < best_result['objective']:
        best_result['objective'] = total_error
        best_result['params'] = {
            'Te': Te,
            'Ne': ne,
            'E_field': E_field,
            'H_drift_gain': H_drift_gain,
            'rates': dict(rate_values)
        }
        best_result['densities'] = results

        log_file = f'optimization_results_C2H2_ultra_boost_v2/best_f{total_error:.1f}_Hdrift{H_drift_gain:.2e}.json'
        run_simulation_with_logging(rate_values, E_field, ne, Te, H_drift_gain, params_updated, log_file)

        print(f"\n  *** NEW BEST: f(x) = {total_error:.2f} at evaluation {objective_function.counter}")
        print(f"      Te: {Te:.2f} eV")
        print(f"      E-field: {E_field:.1f} V/cm")
        print(f"      Ne: {ne:.2e} cm⁻³")
        print(f"      H_drift_gain: {H_drift_gain:.2e} cm⁻³/s")
        print(f"      Charge imbalance: {results['charge_imbalance']*100:.2f}%")
        print(f"      H: {results['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {results['H']/TARGETS['H']:.2f}x)")
        print(f"      CH: {results['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {results['CH']/TARGETS['CH']:.2f}x)")
        print(f"      C2: {results['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {results['C2']/TARGETS['C2']:.2f}x)")
        print(f"      C2H2: {results['C2H2']:.2e} (target: {TARGETS['C2H2']:.2e}, ratio: {results['C2H2']/TARGETS['C2H2']:.2f}x)")
        print(f"      C2H2+H→C2 flux: {results['C2H2_H_C2_flux']:.2e} cm⁻³/s")
        print(f"      Species error: {species_error:.2f}, Charge penalty: {charge_penalty:.2f}")
        print(f"      Saved to: {log_file}\n")

    return total_error


def main():
    print("=" * 80)
    print(" C2H2 ULTRA-BOOST OPTIMIZATION (Version 2)")
    print("=" * 80)
    print("\nUSER INSIGHT:")
    print("  • In CH-suppression solution: C2H2 = 3.4e11 (0.11×), 97% lost to walls")
    print("  • Need to push C2H2 to ~1e13 to feed C2 production via C2H2+H→C2")
    print("\nSTRATEGY:")
    print("  • TARGET C2H2 at 1e13 (3× baseline, 30× current)")
    print("  • REDUCE C2H2 wall losses (tune stick_C2H2, loss_C2H2)")
    print("  • BOOST H (needed for C2H2→C2 pathway)")
    print("  • MAINTAIN CH suppression (keep high CH wall sticking)")
    print("\nNEW WEIGHTS:")
    print("  • H: 30 (boost for pathway)")
    print("  • CH: 80 (MAX - maintain suppression)")
    print("  • C2: 60 (primary target)")
    print("  • C2H2: 40 (VERY HIGH - push to 1e13!)")

    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()
    print(f"  Selected {len(tunable_rates)} key rates")

    db = get_complete_rate_database()

    print("\nSelected rates (top 15):")
    print("-" * 80)
    for i, name in enumerate(tunable_rates[:15], 1):
        if name in db:
            rate = db[name]
            range_factor = rate.max / rate.min if rate.min > 0 else 1.0
            flag = ""
            if 'C2H2' in name and ('loss' in name or 'stick' in name):
                flag = " ← REDUCE for C2H2 boost!"
            elif 'CH_' in name and ('loss' in name or 'stick' in name):
                flag = " ← Keep HIGH for CH suppression!"
            elif 'C2H2_H_C2' in name:
                flag = " ← KEY PATHWAY!"
            print(f"{i:2d}. {name:45s} [{rate.min:.2e}, {rate.max:.2e}] ({range_factor:.1f}x){flag}")
    print(f"    ... and {len(tunable_rates) - 15} more")

    params_base = {
        'E_field': 50,
        'L_discharge': 0.45,
        'ne': 3.3e9,
        'Te': 1.0,
        'H_drift_gain': 3.2e17,
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                    'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                    'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                    'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 4949.6, 'CHPlus': 4949.6, 'H3Plus': 4949.6,
            'H2Plus': 4949.6, 'CH3Minus': 1000, 'HMinus': 1000
        }
    }

    # Build optimization bounds
    bounds = []
    param_names = []

    for name in tunable_rates:
        if name in db:
            rate = db[name]
            bounds.append((rate.min, rate.max))
            param_names.append(name)

    # Add global parameters
    bounds.append((100, 400))      # E_field (V/cm)
    bounds.append((0.7, 8.0))      # Te (eV)
    bounds.append((1e8, 5e9))      # Ne (cm^-3)
    bounds.append((1e16, 5e18))    # H_drift_gain (cm^-3/s)

    print(f"\n Optimization parameters:")
    print(f"  Tunable rates: {len(param_names)}")
    print(f"  E field: [100, 400] V/cm")
    print(f"  Te: [0.7, 8] eV")
    print(f"  Ne: [1e8, 5e9] cm⁻³")
    print(f"  H_drift_gain: [1e16, 5e18] cm⁻³/s")
    print(f"  Total parameters: {len(bounds)}")
    print()

    print(f" Targets:")
    print(f"  H:    {TARGETS['H']:.2e} cm⁻³  (weight: 30.0)")
    print(f"  CH:   {TARGETS['CH']:.2e} cm⁻³  (weight: 80.0) ← MAX!")
    print(f"  C2:   {TARGETS['C2']:.2e} cm⁻³  (weight: 60.0)")
    print(f"  C2H2: {TARGETS['C2H2']:.2e} cm⁻³  (weight: 40.0) ← BOOSTED to 1e13!")
    print()
    print("  + Charge balance penalty: 20.0×")
    print()

    print("=" * 80)
    print(" Running Optimization (20 iterations, pop=8)")
    print("=" * 80)
    print()
    print("Goal: Push C2H2 to ~1e13 to feed C2 production")
    print("      H, CH, C2 all within 0.6-1.4× range")
    print("      Charge balance <100%")
    print()

    # Run optimization
    result = differential_evolution(
        lambda x: objective_function(x, param_names, params_base),
        bounds,
        maxiter=20,
        popsize=8,
        seed=42,
        disp=True,
        workers=1,
        updating='deferred',
        atol=0.01,
        tol=0.01
    )

    print()
    print("=" * 80)
    print(" OPTIMIZATION COMPLETE")
    print("=" * 80)
    print()
    print(f"Best objective value: {best_result['objective']:.2f}")
    print(f"Total evaluations: {objective_function.counter}")
    print()
    print("Final solution:")
    if best_result['densities']:
        for species in ['H', 'CH', 'C2', 'C2H2']:
            ratio = best_result['densities'][species] / TARGETS[species]
            status = '✓' if 0.6 <= ratio <= 1.4 else '✗'
            print(f"  {species:5s}: {best_result['densities'][species]:.2e} / {TARGETS[species]:.2e} = {ratio:.2f}× {status}")
        print()
        print(f"Charge imbalance: {best_result['densities']['charge_imbalance']*100:.2f}%")
        print(f"C2H2+H→C2 flux: {best_result['densities']['C2H2_H_C2_flux']:.2e} cm⁻³/s")


if __name__ == '__main__':
    main()
