#!/usr/bin/env python3
"""
C2H2-BOOST WITH HOT CONDITIONED WALLS

EXPERIMENTAL CONTEXT (from user):
- Gas T ~ 570 K (300°C)
- Cathode likely 200-300°C hot
- Likely has hydrocarbon layer deposited
- Wall sticking should be MUCH LOWER than room-temperature clean walls

Strategy:
1. REDUCE C2H2 wall sticking by 5-10× to allow C2H2 accumulation
2. REDUCE other hydrocarbon wall sticking (hot walls!)
3. MAINTAIN moderate CH wall sticking (still need some suppression)
4. TARGET C2H2 at ~1e13 to drive C2 production via C2H2+H→C2
5. Use H_drift_gain to decouple H from chemistry

Hot Wall Sticking (vs baseline):
- stick_C2H2: [50, 500] (was [500, 2000]) ← 10× lower!
- stick_C2H4: [50, 500] (was [500, 2000]) ← 10× lower!
- stick_C2H6: [50, 400] (was [320, 1600]) ← 6× lower!
- stick_CH: [500, 3000] (was [1250, 6250]) ← Keep moderate for suppression
- stick_C2: [300, 2000] (was [1250, 6250]) ← Lower but controlled

Weights:
- H: 30 (boost for C2H2→C2 pathway)
- CH: 70 (high but not maximum - hot walls help)
- C2: 60 (primary target)
- C2H2: 40 (very high - target ~1e13)

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
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions

# Experimental targets - C2H2 BOOSTED to 1e13!
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
    'C2H2': 1.0e13,  # ← Target based on user's hunch
}

# Create results directory
os.makedirs('optimization_results_C2H2_hot_walls', exist_ok=True)

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
            rate_constant = self.k[self.tags[rxn_idx]]

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
            rate_constant = params['k'][params['tags'][rxn_idx]]

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


def get_hot_wall_bounds():
    """
    Get modified rate bounds for hot (200-300°C) conditioned walls.

    Hot walls have LOWER sticking coefficients:
    - Thermal desorption is significant
    - Hydrocarbon layer reduces sticking probability
    - Radicals have more chance to react in gas phase
    """

    # Get baseline database
    db = get_complete_rate_database()

    # Override wall sticking bounds for hot conditioned walls
    hot_wall_overrides = {
        # C2H2 - KEY! Reduce wall loss to allow accumulation
        'stick_C2H2_9_11': (50.0, 500.0),        # 10× lower than baseline [500, 2000]

        # Other hydrocarbons - hot walls reduce sticking
        'stick_C2H4_9_12': (50.0, 500.0),        # 10× lower than baseline [500, 2000]
        'stick_C2H6_9_14': (50.0, 400.0),        # 6× lower than baseline [320, 1600]
        'stick_C2H5_9_17': (100.0, 1000.0),      # 5× lower than baseline [500, 2000]

        # CH - Keep moderate sticking for suppression
        'stick_CH_9_3': (500.0, 3000.0),         # 2-3× lower than baseline [1250, 6250]
        'stick_CH3_9_2': (600.0, 3000.0),        # 2-5× lower than baseline [1200, 5820]
        'stick_CH2_9_13': (500.0, 2000.0),       # 2-3× lower than baseline [1000, 3000]

        # C2 - Reduce but control
        'stick_C2_9_9': (300.0, 2000.0),         # 3-4× lower than baseline [1250, 6250]

        # H - Hot walls reduce H sticking
        'stick_H_9_1': (100.0, 1000.0),          # 4-10× lower than baseline [389, 3890]

        # C - Atomic carbon sticking
        'stick_C_9_10': (500.0, 3000.0),         # 2-6× lower than baseline [1000, 6000]
    }

    # Build modified database
    modified_bounds = {}
    for rate_name in hot_wall_overrides:
        if rate_name in db:
            min_val, max_val = hot_wall_overrides[rate_name]
            modified_bounds[rate_name] = (min_val, max_val)

    return modified_bounds


def select_tunable_rates():
    """Select rates to tune, focusing on wall sticking for hot walls."""

    # Get chemistry reactions - always tune these
    base_rates = [
        'C2H2_H_C2_H2_H_cm3_7_50',   # C2H2 + H → C2 (KEY PATHWAY!)
        'e_C2H2_C2_H2_cm3_1_16',      # e + C2H2 → C2
        'C2_H_CH_C_cm3_7_6',          # C2 + H → CH (suppress)
        'C_CH_C2_H_cm3_7_4',          # C + CH → C2
        'CH_C_C2_H_cm3_7_9',          # CH + C → C2
        'C2H2_C_C2_CH2_cm3_7_19',     # C2H2 + C → C2 + CH2
    ]

    # Add gas-phase losses
    loss_rates = [
        'loss_C2H2_11_19',            # C2H2 gas loss
        'loss_C2_11_3',               # C2 gas loss
        'loss_CH_11_9',               # CH gas loss
    ]

    # Add all hot-wall sticking rates
    hot_wall_sticking = list(get_hot_wall_bounds().keys())

    # Combine
    all_rates = list(set(base_rates + loss_rates + hot_wall_sticking))

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
        hot_wall_bounds = get_hot_wall_bounds()

        # Apply rate values with hot-wall bounds
        for name, val in rate_values.items():
            if name in k:
                # Use hot-wall bounds if available, otherwise use database bounds
                if name in hot_wall_bounds:
                    min_val, max_val = hot_wall_bounds[name]
                    val = np.clip(val, min_val, max_val)
                elif name in db:
                    val = np.clip(val, db[name].min, db[name].max)
                k[name] = val

        # Enforce bounds for all rates
        for name in k:
            if name in hot_wall_bounds:
                min_val, max_val = hot_wall_bounds[name]
                if k[name] < min_val:
                    k[name] = min_val
                elif k[name] > max_val:
                    k[name] = max_val
            elif name in db:
                rate_db = db[name]
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

    # HOT-WALL C2H2-BOOST WEIGHTS
    weights = {
        'H': 30.0,    # Boost H for C2H2→C2 pathway
        'CH': 70.0,   # High but not max (hot walls help)
        'C2': 60.0,   # Primary target
        'C2H2': 40.0  # Very high - target ~1e13!
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

        log_file = f'optimization_results_C2H2_hot_walls/best_f{total_error:.1f}_Hdrift{H_drift_gain:.2e}.json'
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
    print(" C2H2-BOOST WITH HOT CONDITIONED WALLS")
    print("=" * 80)
    print("\nEXPERIMENTAL CONTEXT:")
    print("  • Gas T ~ 570 K (300°C)")
    print("  • Cathode likely 200-300°C hot")
    print("  • Likely has hydrocarbon layer deposited")
    print("  • Wall sticking should be MUCH LOWER than baseline")
    print("\nHOT-WALL STICKING COEFFICIENTS (vs baseline):")
    print("  • stick_C2H2: [50, 500] (was [500, 2000]) - 10× lower!")
    print("  • stick_C2H4: [50, 500] (was [500, 2000]) - 10× lower!")
    print("  • stick_C2H6: [50, 400] (was [320, 1600]) - 6× lower!")
    print("  • stick_CH: [500, 3000] (was [1250, 6250]) - keep moderate")
    print("  • stick_C2: [300, 2000] (was [1250, 6250]) - lower but controlled")
    print("\nSTRATEGY:")
    print("  • Reduced C2H2 wall losses → C2H2 can accumulate to ~1e13")
    print("  • More C2H2+H→C2 conversion drives C2 production")
    print("  • Moderate CH wall sticking maintains suppression")
    print("\nWEIGHTS:")
    print("  • H: 30 (boost for pathway)")
    print("  • CH: 70 (high but not max - hot walls help)")
    print("  • C2: 60 (primary target)")
    print("  • C2H2: 40 (very high - target ~1e13!)")

    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()
    print(f"  Selected {len(tunable_rates)} key rates")

    db = get_complete_rate_database()
    hot_wall_bounds = get_hot_wall_bounds()

    print("\nSelected rates with HOT-WALL bounds:")
    print("-" * 80)
    count = 0
    for name in tunable_rates:
        count += 1
        if name in hot_wall_bounds:
            min_val, max_val = hot_wall_bounds[name]
            range_factor = max_val / min_val if min_val > 0 else 1.0
            flag = " ← HOT WALL! (reduced)"
            print(f"{count:2d}. {name:45s} [{min_val:.2e}, {max_val:.2e}] ({range_factor:.1f}x){flag}")
        elif name in db:
            rate = db[name]
            range_factor = rate.max / rate.min if rate.min > 0 else 1.0
            flag = ""
            if 'C2H2_H_C2' in name:
                flag = " ← KEY PATHWAY!"
            print(f"{count:2d}. {name:45s} [{rate.min:.2e}, {rate.max:.2e}] ({range_factor:.1f}x){flag}")
        if count >= 20:
            print(f"    ... and {len(tunable_rates) - 20} more")
            break

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
        if name in hot_wall_bounds:
            min_val, max_val = hot_wall_bounds[name]
            bounds.append((min_val, max_val))
            param_names.append(name)
        elif name in db:
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
    print(f"  CH:   {TARGETS['CH']:.2e} cm⁻³  (weight: 70.0)")
    print(f"  C2:   {TARGETS['C2']:.2e} cm⁻³  (weight: 60.0)")
    print(f"  C2H2: {TARGETS['C2H2']:.2e} cm⁻³  (weight: 40.0) ← Target ~1e13!")
    print()
    print("  + Charge balance penalty: 20.0×")
    print()

    print("=" * 80)
    print(" Running Optimization (20 iterations, pop=8)")
    print("=" * 80)
    print()
    print("Goal: Hot walls allow C2H2 to accumulate to ~1e13")
    print("      More C2H2+H→C2 conversion drives C2 production")
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
