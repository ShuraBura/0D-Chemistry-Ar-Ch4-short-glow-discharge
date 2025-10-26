#!/usr/bin/env python3
"""
V7b: BALANCED OPTIMIZATION - Species + Charge Balance

STRATEGY: Balance between species accuracy AND charge neutrality.
This seeks the optimal trade-off between targets and physics constraints.

KEY IMPROVEMENTS FROM V6:
1. NEW REACTIONS:
   - ArH+ + e → Ar + H (fixes 7.8% of ions missing recombination)
   - e + CH4 → CH3 + H⁻ (dissociative attachment)

2. EXPANDED TUNING:
   - All primary, secondary, AND tertiary chemistry pathways
   - Based on pathway analysis: H → CH3 → C2H2 → C2 cascade
   - 150+ tunable reactions (vs 120 in V6)

3. MODERATE CHARGE BALANCE PENALTY:
   - CHARGE_BALANCE_WEIGHT = 150 (between V5's 75 and V6's 300)
   - Balanced focus on species + charge

4. MORE ITERATIONS:
   - 25 iterations (vs 20 in V6) for better convergence

Expected result: Good species accuracy with improved charge balance.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import time
import signal
from contextlib import contextmanager
import os
import json

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database, get_tunable_rates_for_target
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

TARGET_ION_E_RATIO = 4.5
CHARGE_BALANCE_WEIGHT = 150.0  # *** BALANCED (between V5's 75 and V6's 300) ***

os.makedirs('optimization_results_v7b_balanced', exist_ok=True)

iteration_counter = 0
best_result = {'objective': 1e10, 'params': None, 'densities': None}

class TimeoutException(Exception):
    pass

@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

def select_tunable_rates():
    """V7b: COMPREHENSIVE + secondary/tertiary pathways"""
    db = get_complete_rate_database()
    selected = set()

    # Base: All rates with range > 1.5x
    selected |= {name for name, rate in db.items() if (rate.max / rate.min > 1.5) and (rate.min > 0)}

    # All wall/volume losses
    selected |= {name for name in db.keys() if name.startswith('stick_')}
    selected |= {name for name in db.keys() if name.startswith('loss_')}

    # All electron-impact (except ionization)
    ionization_keywords = ['CH4Plus', 'CH3Plus', 'ArPlus', 'CH5Plus', 'CHPlus', 'CH2Plus',
                          'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus']
    selected |= {name for name in db.keys() if name.startswith('e_') and
                 not any(ion in name for ion in ionization_keywords)}

    # ALL dissociative recombination (16 reactions now with ArH+)
    dissoc_recomb = [
        'e_CH4Plus_CH3_H_cm3_6_4', 'e_CH4Plus_CH2_H2_cm3_6_9',
        'e_CH4Plus_CH_H2_H_cm3_6_11', 'e_CH4Plus_C_2H2_cm3_6_13',
        'CH5Plus_e_CH4_H_cm3_6_3', 'CH5Plus_e_CH3_H2_cm3_6_8',
        'CH5Plus_e_CH2_H2_H_cm3_6_10', 'CH5Plus_e_CH3_2H_cm3_6_12',
        'ArPlus_e_Ar_cm3_6_1', 'CH3Plus_e_CH3_cm3_6_2',
        'C2H5Plus_e_C2H4_H_cm3_6_14', 'C2H4Plus_e_C2H2_H2_cm3_6_15',
        'C2H3Plus_e_C2H2_H_cm3_6_16', 'C2HPlus_e_C2_H_cm3_6_18',
        'C2H5Plus_e_C2H4_H_cm3_6_28',
        'ArHPlus_e_Ar_H_cm3_6_29',  # NEW IN V7!
        'e_CH4_CH3_HMinus_cm3_8_1',  # NEW IN V7!
    ]
    selected |= set(dissoc_recomb)

    # SECONDARY/TERTIARY PATHWAYS from analysis
    secondary_tertiary = [
        # H chemistry cascade
        'H_CH4_CH3_H2_cm3_7_25',      # 77% of H loss, 99.5% of CH3 production
        'CH3_CH3_C2H2_H2_H2_cm3_7_30', # 88% of C2H2 production
        'CH2_CH3_C2H2_H_H2_cm3_7_56',  # 5% of C2H2 production
        'C2H2_H_C2_H2_H_cm3_7_50',     # 96% of C2 production

        # CH pathways
        'CH2_H_CH_H2_cm3_7_1',         # 52% CH production
        'C2_H_CH_C_cm3_7_6',          # 22% CH production, 20% C2 loss
        'C_H_CH_cm3_7_12',            # 10% CH production
        'CH_CH4_C2H4_H_cm3_7_20',     # 34% CH loss
        'CH_CH4_CH2_CH3_cm3_7_39',    # 23% CH loss
        'CH_CH3_C2H3_H_cm3_7_10',     # CH+CH3 reactions
        'CH_CH3_C2H4_cm3_7_5',

        # C2 pathways
        'CH2_CH2_C2_H2_H2_cm3_7_58',  # 2% C2 production
        'C2H2_C_C2_CH2_cm3_7_51',     # 1% C2 production
        'CH3_H_CH2_H2_cm3_7_36',      # CH3 + H reactions

        # Ion chemistry (tertiary)
        'ArPlus_CH4_CH3Plus_H_cm3_5_9', # Ion-molecule
        'CH4Plus_H2_CH5Plus_H_cm3_5_6',
        'ArPlus_CH4_Ar_CH4Plus_cm3_5_5',

        # Ar* pathways
        'ArStar_CH4_CH3_H_cm3_3_1',
        'e_Ar_ArStar_cm3_1_7',
    ]
    selected |= set(secondary_tertiary)

    # Target-specific rates
    selected |= set(get_tunable_rates_for_target('H').keys())
    selected |= set(get_tunable_rates_for_target('CH').keys())
    selected |= set(get_tunable_rates_for_target('C2').keys())

    # Flagged rates
    selected |= {name for name, rate in db.items() if rate.flag}

    selected = {name for name in selected if name in db}
    selected_names = sorted(selected)

    print(f"  V7a selection: {len(selected_names)} rates")
    return selected_names

def analyze_chemistry(y_final, species, params):
    rate_constants = np.array([params['k'][tag] for tag in params['tags']])
    reactions = params['R']
    reaction_rates = rate_constants.copy()
    for rxn_idx, reaction in enumerate(reactions):
        react_species = np.where(reaction.reactants > 0)[0]
        for sp_idx in react_species:
            reaction_rates[rxn_idx] *= y_final[sp_idx] ** reaction.reactants[sp_idx]

    important_species = ['H', 'CH', 'C2', 'CH3', 'C2H2', 'CH2', 'e', 'ArPlus', 'CH4Plus', 'CH5Plus', 'ArHPlus']
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
                production_rxns.append((params['tags'][rxn_idx], net_change * reaction_rates[rxn_idx]))
            elif net_change < 0:
                loss_rxns.append((params['tags'][rxn_idx], abs(net_change) * reaction_rates[rxn_idx]))

        production_rxns.sort(key=lambda x: x[1], reverse=True)
        loss_rxns.sort(key=lambda x: x[1], reverse=True)

        analysis[target_sp] = {
            'density': float(y_final[sp_idx]),
            'production': sum(x[1] for x in production_rxns),
            'loss': sum(x[1] for x in loss_rxns),
            'top_production': [(tag, float(val)) for tag, val in production_rxns[:10]],
            'top_loss': [(tag, float(val)) for tag, val in loss_rxns[:10]]
        }

    return analysis

def run_simulation_with_logging(rate_values, E_field, ne, params_base, log_file=None):
    try:
        params = params_base.copy()
        params['E_field'] = E_field
        params['ne'] = ne
        k = define_rates(params)
        db = get_complete_rate_database()

        for name, val in rate_values.items():
            if name in k and name in db:
                k[name] = np.clip(val, db[name].min, db[name].max)

        for name, rate_db in db.items():
            if name in k:
                k[name] = np.clip(k[name], rate_db.min, rate_db.max)

        params['k'] = k
        params['R'], params['tags'] = build_reactions(params)

        species = params['species']
        ns = len(species)
        y0 = np.ones(ns) * 1e3

        def set_density(name, value):
            try:
                y0[species.index(name)] = value
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
        set_density('HMinus', 1e3)

        ode_func = PlasmaODE_Optimized(params)

        try:
            with time_limit(30):
                sol = solve_ivp(ode_func, (0, 100), y0, method='BDF',
                               rtol=1e-5, atol=1e-6, max_step=10.0)
        except TimeoutException:
            return None

        if not sol.success:
            return None

        y_final = sol.y[:, -1]
        all_densities = {species[i]: float(y_final[i]) for i in range(ns)}

        results = {
            'H': all_densities.get('H', 0),
            'CH': all_densities.get('CH', 0),
            'C2': all_densities.get('C2', 0),
            'all_densities': all_densities
        }

        if log_file:
            chemistry = analyze_chemistry(y_final, species, params)
            log_data = {
                'Ne': ne,
                'E_field': E_field,
                'rate_values': {k: float(v) for k, v in rate_values.items()},
                'target_densities': {'H': results['H'], 'CH': results['CH'], 'C2': results['C2']},
                'all_densities': all_densities,
                'chemistry_analysis': chemistry
            }
            with open(log_file, 'w') as f:
                json.dump(log_data, f, indent=2)

        return results
    except Exception:
        return None

def objective_function(x, param_names, params_base):
    global iteration_counter, best_result

    if not hasattr(objective_function, 'counter'):
        objective_function.counter = 0
        objective_function.start_time = time.time()

    objective_function.counter += 1

    if objective_function.counter % 10 == 0:
        elapsed = time.time() - objective_function.start_time
        print(f"  [{objective_function.counter} evaluations, {elapsed/60:.1f} min elapsed]")

    rate_values = {name: val for name, val in zip(param_names[:-2], x[:-2])}
    E_field = x[-2]
    ne = x[-1]

    params_updated = params_base.copy()
    params_updated['ne'] = ne

    results = run_simulation_with_logging(rate_values, E_field, ne, params_updated)

    if results is None:
        return 1e10

    # BALANCED: SPECIES ERROR + CHARGE PENALTY
    weights = {'H': 1.0, 'CH': 30.0, 'C2': 3.0}

    error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        error += weights[species] * rel_error ** 2

    # ADD CHARGE BALANCE PENALTY (moderate weight = 150)
    all_dens = results['all_densities']
    positive_ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                     'H3Plus', 'CHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                     'C2H3Plus', 'C2HPlus']
    total_positive = sum(all_dens.get(ion, 0) for ion in positive_ions)
    ion_e_ratio = total_positive / ne if ne > 0 else 0

    # Charge balance penalty
    charge_error = abs(ion_e_ratio - TARGET_ION_E_RATIO) / TARGET_ION_E_RATIO
    error += CHARGE_BALANCE_WEIGHT * charge_error

    if error < best_result['objective']:
        best_result['objective'] = error
        best_result['params'] = {'Ne': ne, 'E_field': E_field, 'rates': dict(rate_values)}
        best_result['densities'] = results

        log_file = f'optimization_results_v7b_balanced/best_iteration_{iteration_counter:04d}_f{error:.1f}.json'
        run_simulation_with_logging(rate_values, E_field, ne, params_updated, log_file)

        print(f"\n  *** NEW BEST: f(x) = {error:.2f} at evaluation {objective_function.counter}")
        print(f"      H: {results['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {results['H']/TARGETS['H']:.2f}x)")
        print(f"      CH: {results['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {results['CH']/TARGETS['CH']:.2f}x)")
        print(f"      C2: {results['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {results['C2']/TARGETS['C2']:.2f}x)")
        print(f"      Ion/e: {ion_e_ratio:.2f}x (target: {TARGET_ION_E_RATIO:.1f}x, weight: {CHARGE_BALANCE_WEIGHT})")
        print(f"      Saved to: {log_file}\n")

    return error

def main():
    print("=" * 80)
    print(" V7b OPTIMIZATION - BALANCED SPECIES + CHARGE")
    print("=" * 80)
    print("\nSTRATEGY:")
    print("  ✓ BALANCED charge balance penalty (CHARGE_WEIGHT = 150)")
    print("  ✓ Trade-off between H, CH, C2 targets AND charge neutrality")
    print("  ✓ New reactions: ArH+ + e, e + CH4 → CH3 + H⁻")
    print("  ✓ Expanded tuning: 150+ reactions (primary+secondary+tertiary)")
    print("  ✓ 25 iterations for convergence")
    print("\nExpected: Good species accuracy + improved charge balance")

    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()

    db = get_complete_rate_database()

    params_base = {
        'E_field': 50, 'L_discharge': 0.45, 'ne': 1.4e9,
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                    'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                    'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                    'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'C2H2Star'],
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'CH3Minus': 3000, 'HMinus': 3000
        }
    }

    bounds = []
    param_names = []
    for name in tunable_rates:
        if name in db:
            bounds.append((db[name].min, db[name].max))
            param_names.append(name)

    bounds.append((5.0, 300.0))
    param_names.append('E_field')
    bounds.append((1.0e9, 5.0e9))
    param_names.append('ne')

    print(f"\nOptimization parameters: {len(param_names)} total")
    print(f"  Tunable rates: {len(param_names) - 2}")

    print("\n" + "=" * 80)
    print(" Running V7b Optimization (25 iterations, pop=6)")
    print("=" * 80)

    start_time = time.time()

    result = differential_evolution(
        objective_function,
        bounds,
        args=(param_names, params_base),
        strategy='best1bin',
        maxiter=25,
        popsize=6,
        tol=0.01,
        atol=0.0,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=42,
        disp=True,
        workers=1,
        updating='deferred',
        polish=True
    )

    elapsed = time.time() - start_time
    print(f"\nV7b COMPLETE ({elapsed/60:.1f} minutes)")
    print(f"Final objective: {result.fun:.6f}")

    with open('optimization_results_v7b_balanced/FINAL_RESULT.json', 'w') as f:
        json.dump({
            'objective': result.fun,
            'parameters': {'Ne': result.x[-1], 'E_field': result.x[-2],
                          'rates': {name: result.x[i] for i, name in enumerate(param_names[:-2])}},
            'success': result.success,
            'iterations': result.nit,
            'evaluations': result.nfev
        }, f, indent=2)

if __name__ == '__main__':
    main()
