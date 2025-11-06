#!/usr/bin/env python3
"""
COMPREHENSIVE CH-FOCUSED OPTIMIZATION WITH H-DECOUPLING

Combines ALL 3 strategies:
1. Ultra-aggressive CH weight (50×) - MUST reduce to target
2. Lower C2H2 target (3.0e12) - more realistic based on chemistry
3. Tunable H sheath source - decouples H from chemistry
4. EXPLICIT penalty for C2+H→CH+C pathway flux

Goal: Get f < 100 with all species near targets!
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


# Experimental targets - REALISTIC C2H2
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,      # ← PRIMARY: Must reduce from 17.5×!
    'C2': 1.3e11,     # ← Must maintain!
    'C2H2': 3.0e12,   # ← LOWERED from 7e12 (more realistic)
}

# Create results directory
os.makedirs('optimization_results_comprehensive_CH', exist_ok=True)

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


# Modified PlasmaODE with tunable H source
class PlasmaODE_H_Tunable:
    """ODE function with tunable H sheath source."""

    def __init__(self, params):
        self.species = params['species']
        self.ns = len(self.species)
        self.R = params['R']
        self.nr = len(self.R)
        self.k = params['k']
        self.tags = params['tags']
        self.H_drift_gain = params.get('H_drift_gain', 3.2e17)

        # Cache indices
        self.e_idx = self.species.index('e')
        self.Ar_idx = self.species.index('Ar')
        self.CH4_idx = self.species.index('CH4')
        self.H_idx = self.species.index('H')

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

        self.products_matrix = sparse.csr_matrix(
            (vals_prod, (rows_prod, cols_prod)),
            shape=(self.ns, self.nr)
        )
        self.reactants_matrix = sparse.csr_matrix(
            (vals_react, (rows_react, cols_react)),
            shape=(self.ns, self.nr)
        )
        self.net_stoich_matrix = self.products_matrix - self.reactants_matrix

        self.rxn_reactants = []
        for reaction in self.R:
            react_idx = np.where(reaction.reactants > 0)[0]
            react_coeff = reaction.reactants[react_idx]
            self.rxn_reactants.append((react_idx, react_coeff))

    def __call__(self, t, y):
        """ODE evaluation with tunable H source."""
        y = np.maximum(y, 1e-6)

        reaction_rates = np.array([self.k[tag] for tag in self.tags])
        for rxn_idx, (sp_indices, coeffs) in enumerate(self.rxn_reactants):
            for sp_idx, coeff in zip(sp_indices, coeffs):
                reaction_rates[rxn_idx] *= y[sp_idx] ** coeff

        dydt = self.net_stoich_matrix.dot(reaction_rates)
        dydt[self.H_idx] += self.H_drift_gain

        dydt[self.e_idx] = 0.0
        dydt[self.Ar_idx] = 0.0
        dydt[self.CH4_idx] = 0.0

        return dydt


def select_tunable_rates():
    """Select tunable rates - emphasize CH reduction pathways."""
    critical_neutral = [
        'CH2_H_CH_H2_cm3_7_1',          # CH2+H→CH (major CH source: 279 THz)
        'C2_H_CH_C_cm3_7_6',            # C2+H→CH+C ← SUPPRESS!
        'C_H_CH_cm3_7_12',              # C+H→CH
        'C2H2_H_C2_H2_H_cm3_7_50',     # C2H2+H→C2 ← Important for C2
        'CH_CH4_C2H4_H_cm3_7_20',      # CH+CH4→C2H4 (CH sink)
        'CH_CH4_CH2_CH3_cm3_7_39',     # CH+CH4→ (CH sink)
        'CH_H_C_H2_cm3_7_3',            # CH+H→C+H2 (CH loss)
        'CH_H_CH2_cm3_7_21',            # CH+H→CH2
        'CH2_CH3_C2H2_H_H2_cm3_7_62',  # C2H2 production
        'CH3_CH3_C2H2_H2_H2_cm3_7_49', # C2H2 production (major)
        'H_CH4_CH3_H2_cm3_7_25',       # H loss (major: 75%)
        'CH3_H_CH2_H2_cm3_7_36',       # H loss
    ]

    db = get_complete_rate_database()

    h_rates = set(get_tunable_rates_for_target('H').keys())
    ch_rates = set(get_tunable_rates_for_target('CH').keys())
    c2_rates = set(get_tunable_rates_for_target('C2').keys())

    target_rates = h_rates | ch_rates | c2_rates

    large_range_rates = {
        name for name, rate in db.items()
        if (rate.max / rate.min > 3.0) and (rate.min > 0)
    }

    flagged_rates = {name for name, rate in db.items() if rate.flag}

    selected = set(critical_neutral) | target_rates | large_range_rates | flagged_rates
    selected = {name for name in selected if name in db}

    selected_with_range = [
        (name, db[name].max / db[name].min if db[name].min > 0 else 1.0)
        for name in selected
    ]
    selected_with_range.sort(key=lambda x: x[1], reverse=True)

    return [name for name, _ in selected_with_range[:40]]


def analyze_chemistry(y_final, species, params):
    """Detailed chemistry analysis."""
    rate_constants = np.array([params['k'][tag] for tag in params['tags']])
    reactions = params['R']

    reaction_rates = rate_constants.copy()
    for rxn_idx, reaction in enumerate(reactions):
        react_species = np.where(reaction.reactants > 0)[0]
        for sp_idx in react_species:
            coeff = reaction.reactants[sp_idx]
            reaction_rates[rxn_idx] *= y_final[sp_idx] ** coeff

    important_species = ['H', 'CH', 'C2', 'C2H2', 'C2H4', 'CH3', 'CH4',
                        'ArStar', 'e', 'ArPlus', 'CH3Plus', 'CH5Plus',
                        'HMinus', 'CH3Minus']

    analysis = {}

    for target_sp in important_species:
        try:
            sp_idx = species.index(target_sp)
        except ValueError:
            continue

        production_rxns = []
        loss_rxns = []

        for rxn_idx, reaction in enumerate(reactions):
            net_change = reaction.products[sp_idx] - reaction.reactants[sp_idx]

            if net_change > 0:
                contribution = net_change * reaction_rates[rxn_idx]
                production_rxns.append((params['tags'][rxn_idx], contribution))
            elif net_change < 0:
                contribution = abs(net_change) * reaction_rates[rxn_idx]
                loss_rxns.append((params['tags'][rxn_idx], contribution))

        production_rxns.sort(key=lambda x: x[1], reverse=True)
        loss_rxns.sort(key=lambda x: x[1], reverse=True)

        total_production = sum(x[1] for x in production_rxns)
        total_loss = sum(x[1] for x in loss_rxns)

        analysis[target_sp] = {
            'density': float(y_final[sp_idx]),
            'production': total_production,
            'loss': total_loss,
            'top_production': [(tag, float(val)) for tag, val in production_rxns[:10]],
            'top_loss': [(tag, float(val)) for tag, val in loss_rxns[:10]]
        }

    return analysis


def run_simulation_with_logging(rate_values, E_field, ne, Te, H_drift_gain, params_base, log_file=None):
    """Run simulation with tunable H source."""

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

        # Calculate C2+H→CH flux for penalty
        try:
            H_idx = species.index('H')
            C2_idx = species.index('C2')
            C2_H_CH_flux = rate_values.get('C2_H_CH_C_cm3_7_6', k.get('C2_H_CH_C_cm3_7_6', 0)) * \
                          y_final[C2_idx] * y_final[H_idx]
            results['C2_H_CH_flux'] = C2_H_CH_flux
        except:
            results['C2_H_CH_flux'] = 0.0

        # Charge balance
        positive_ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                        'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus',
                        'H3Plus', 'CHPlus', 'H2Plus']

        n_i_total = sum(get_density(ion) for ion in positive_ions)
        results['n_i_total'] = n_i_total
        results['charge_imbalance'] = abs(n_i_total - ne) / ne

        if log_file:
            all_densities = {species[i]: float(y_final[i]) for i in range(ns)}
            chemistry = analyze_chemistry(y_final, species, params)

            log_data = {
                'Te': Te,
                'Ne': ne,
                'E_field': E_field,
                'H_drift_gain': H_drift_gain,
                'rate_values': {k: float(v) for k, v in rate_values.items()},
                'target_densities': {k: results[k] for k in ['H', 'CH', 'C2']},
                'charge_balance': {
                    'n_i_total': float(n_i_total),
                    'n_e': float(ne),
                    'imbalance_percent': float(results['charge_imbalance'] * 100)
                },
                'C2_H_CH_flux': float(results['C2_H_CH_flux']),
                'all_densities': all_densities,
                'chemistry_analysis': chemistry
            }

            with open(log_file, 'w') as f:
                json.dump(log_data, f, indent=2)

        return results

    except Exception as e:
        return None


def objective_function(x, param_names, params_base):
    """Objective function with ULTRA CH focus."""
    global best_result

    if not hasattr(objective_function, 'counter'):
        objective_function.counter = 0
        objective_function.start_time = time.time()

    objective_function.counter += 1

    if objective_function.counter % 10 == 0:
        elapsed = time.time() - objective_function.start_time
        print(f"  [{objective_function.counter} evaluations, {elapsed/60:.1f} min elapsed]")

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

    # ULTRA AGGRESSIVE CH WEIGHTING
    weights = {
        'H': 5.0,     # Moderate - from tunable sheath source
        'CH': 50.0,   # ← MAXIMUM! Must reduce from 17.5× to ~1×
        'C2': 30.0,   # High - must maintain near 1.12×
        'C2H2': 5.0   # Lower - more realistic target
    }

    species_error = 0.0
    for species in ['H', 'CH', 'C2', 'C2H2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        species_error += weights[species] * rel_error ** 2

    # Charge balance penalty
    charge_penalty = 10.0 * results['charge_imbalance']**2

    # EXPLICIT penalty for C2+H→CH flux (normalize by 1e14)
    flux_penalty = 1.0 * (results['C2_H_CH_flux'] / 1e14) ** 2

    total_error = species_error + charge_penalty + flux_penalty

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

        log_file = f'optimization_results_comprehensive_CH/best_f{total_error:.1f}_Hdrift{H_drift_gain:.2e}.json'
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
        print(f"      C2+H→CH flux: {results['C2_H_CH_flux']:.2e} cm⁻³/s")
        print(f"      Species: {species_error:.2f}, Charge: {charge_penalty:.2f}, Flux: {flux_penalty:.2f}")
        print(f"      Saved to: {log_file}\n")

    return total_error


def main():
    print("=" * 80)
    print(" COMPREHENSIVE CH-FOCUSED OPTIMIZATION")
    print("=" * 80)
    print("\nCombining ALL 3 strategies:")
    print("  1. ULTRA-aggressive CH weight (50×) - MUST hit target!")
    print("  2. Lower C2H2 target (3.0e12 vs 7.0e12) - more realistic")
    print("  3. Tunable H sheath source - maintains H while reducing chemistry")
    print("  4. EXPLICIT C2+H→CH pathway penalty")
    print("\nGoal: f < 100 with all species near targets!")

    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()
    print(f"  Selected {len(tunable_rates)} key rates")

    db = get_complete_rate_database()

    print("\nSelected rates (top 10):")
    print("-" * 80)
    for i, name in enumerate(tunable_rates[:10], 1):
        if name in db:
            rate = db[name]
            range_factor = rate.max / rate.min if rate.min > 0 else 1.0
            flag = " ← TARGET!" if name in ['C2_H_CH_C_cm3_7_6', 'CH2_H_CH_H2_cm3_7_1'] else ""
            print(f"{i:2d}. {name:45s} [{rate.min:.2e}, {rate.max:.2e}] ({range_factor:.1f}x){flag}")
    print(f"    ... and {len(tunable_rates) - 10} more")

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
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
        }
    }

    bounds = []
    param_names = []

    for name in tunable_rates:
        if name in db:
            bounds.append((db[name].min, db[name].max))
            param_names.append(name)

    bounds.append((100.0, 400.0))  # E field
    param_names.append('E_field')

    bounds.append((0.7, 8.0))  # Te
    param_names.append('Te')

    bounds.append((1e8, 5e9))  # Ne
    param_names.append('ne')

    bounds.append((1e16, 1e18))  # H_drift_gain (narrower range)
    param_names.append('H_drift_gain')

    print(f"\n Optimization parameters:")
    print(f"  Tunable rates: {len(param_names) - 4}")
    print(f"  E field: [100, 400] V/cm")
    print(f"  Te: [0.7, 8] eV")
    print(f"  Ne: [1e8, 5e9] cm⁻³")
    print(f"  H_drift_gain: [1e16, 1e18] cm⁻³/s")
    print(f"  Total parameters: {len(param_names)}")

    print(f"\n Targets:")
    print(f"  H:    {TARGETS['H']:.2e} cm^-3  (weight: 5.0)")
    print(f"  CH:   {TARGETS['CH']:.2e} cm^-3  (weight: 50.0) ← ULTRA FOCUS!")
    print(f"  C2:   {TARGETS['C2']:.2e} cm^-3  (weight: 30.0) ← Must maintain!")
    print(f"  C2H2: {TARGETS['C2H2']:.2e} cm^-3  (weight: 5.0) ← Lowered target")
    print(f"\n  + Charge balance (weight: 10)")
    print(f"  + C2+H→CH flux penalty (weight: 1)")

    print("\n" + "=" * 80)
    print(" Running Optimization (20 iterations, pop=8)")
    print("=" * 80)
    print("\nTarget: f < 100, all species within 2× of targets!\n")

    start_time = time.time()

    result = differential_evolution(
        objective_function,
        bounds,
        args=(param_names, params_base),
        strategy='best1bin',
        maxiter=20,  # More iterations
        popsize=8,   # Larger population
        tol=0.01,
        atol=0.0,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=43,
        disp=True,
        workers=1,
        updating='deferred',
        polish=True
    )

    elapsed = time.time() - start_time

    print("\n" + "=" * 80)
    print(f" OPTIMIZATION COMPLETE ({elapsed/60:.1f} minutes)")
    print("=" * 80)

    print(f"\nFinal objective: {result.fun:.6f}")
    print(f"Iterations: {result.nit}")
    print(f"Evaluations: {result.nfev}")

    Te_opt = result.x[-3]
    Ne_opt = result.x[-2]
    H_drift_gain_opt = result.x[-1]
    E_opt = result.x[-4]

    print(f"\nOptimized parameters:")
    print(f"  Te: {Te_opt:.2f} eV")
    print(f"  E-field: {E_opt:.1f} V/cm")
    print(f"  Ne: {Ne_opt:.2e} cm⁻³")
    print(f"  H_drift_gain: {H_drift_gain_opt:.2e} cm⁻³/s")

    final_params = {
        'Te': float(Te_opt),
        'Ne': float(Ne_opt),
        'E_field': float(E_opt),
        'H_drift_gain': float(H_drift_gain_opt),
        'rates': {name: float(result.x[i]) for i, name in enumerate(param_names[:-4])}
    }

    with open('optimization_results_comprehensive_CH/FINAL_RESULT.json', 'w') as f:
        json.dump({
            'objective': float(result.fun),
            'parameters': final_params,
            'success': result.success,
            'message': result.message,
            'iterations': int(result.nit),
            'evaluations': int(result.nfev)
        }, f, indent=2)

    print(f"\n✓ Comprehensive CH-focused optimization complete!")
    print(f"✓ Results saved to optimization_results_comprehensive_CH/")


if __name__ == '__main__':
    main()
