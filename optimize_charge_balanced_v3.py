#!/usr/bin/env python3
"""
CHARGE-BALANCED optimization with expanded tunable rates.

KEY IMPROVEMENTS:
1. CHARGE BALANCE constraint: target ion/e ratio = 4.5x (range: 3-6x)
2. Lower E field minimum (5 V/cm) to reduce ion loss to walls
3. Includes critical neutral reactions:
   - CH2+H→CH (48.7% of CH production)
   - C2+H→CH (24% of CH production) ← CRITICAL!
   - C+H→CH (12.9% of CH production)
   - C2H2+H→C2 (96.5% of C2 production)

Physics: Lower E → less ion drift to walls → higher ion density → charge balance!
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


# Experimental targets
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

# Charge balance target (ion/electron ratio)
TARGET_ION_E_RATIO = 4.5  # Middle of 3-6x range
CHARGE_BALANCE_WEIGHT = 200.0  # Penalty weight for charge imbalance (INCREASED from 5.0!)

# Create results directory
os.makedirs('optimization_results_charge_balanced_v3', exist_ok=True)

# Global counter for iterations
iteration_counter = 0
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


def select_tunable_rates():
    """
    V3: MASSIVELY EXPANDED tunable rates for more degrees of freedom!

    Strategy: Tune 100+ rates to give optimizer flexibility to satisfy
    BOTH charge balance AND species targets simultaneously.

    Includes:
    - ALL wall sticking rates (28 rates)
    - ALL volume loss rates (25 rates)
    - Most electron-impact (except ionization: ~26 rates)
    - All rates with range > 2x (~77 rates)
    - Critical neutral reactions for H, CH, C2
    - Ion loss reactions (not production)
    """
    db = get_complete_rate_database()

    # Start with empty set
    selected = set()

    # 1. ALL wall sticking rates (retain flexibility for all species)
    wall_sticking = {name for name in db.keys() if name.startswith('stick_')}
    selected |= wall_sticking

    # 2. ALL volume loss rates (uncertain parameters)
    volume_loss = {name for name in db.keys() if name.startswith('loss_')}
    selected |= volume_loss

    # 3. ALL electron-impact rates EXCEPT ionization
    #    Ionization excluded per user guidance: e+X→X+
    ionization_keywords = ['CH4Plus', 'CH3Plus', 'ArPlus', 'CH5Plus',
                          'CHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                          'C2H3Plus', 'C2HPlus', 'H3Plus']

    electron_impact = {
        name for name in db.keys()
        if name.startswith('e_') and
        not any(ion in name for ion in ionization_keywords)
    }
    selected |= electron_impact

    # 4. All rates with range > 2x (moderately tunable)
    large_range = {
        name for name, rate in db.items()
        if (rate.max / rate.min > 2.0) and (rate.min > 0)
    }
    selected |= large_range

    # 5. Target-specific rates
    h_rates = set(get_tunable_rates_for_target('H').keys())
    ch_rates = set(get_tunable_rates_for_target('CH').keys())
    c2_rates = set(get_tunable_rates_for_target('C2').keys())
    selected |= (h_rates | ch_rates | c2_rates)

    # 6. Flagged rates (uncertain in literature)
    flagged_rates = {name for name, rate in db.items() if rate.flag}
    selected |= flagged_rates

    # 7. Critical neutral reactions (known to be important)
    critical_neutral = [
        'CH2_H_CH_H2_cm3_7_1',          # CH2+H→CH (48.7% of CH production)
        'C2_H_CH_C_cm3_7_6',            # C2+H→CH (24% of CH production) ← KEY!
        'C_H_CH_cm3_7_12',              # C+H→CH (12.9% of CH production)
        'C2H2_H_C2_H2_H_cm3_7_50',     # C2H2+H→C2 (96.5% of C2 production)
        'CH_CH4_C2H4_H_cm3_7_20',      # CH+CH4→C2H4 (CH loss)
        'CH_CH4_CH2_CH3_cm3_7_39',     # CH+CH4→CH2+CH3 (CH loss)
        'CH_H_C_H2_cm3_7_3',            # CH+H→C+H2 (CH loss)
    ]
    selected |= set(critical_neutral)

    # 8. Critical ion LOSS (not production!) for charge balance
    critical_ion_loss = [
        'stick_ArPlus_9_4',              # Ar+ wall loss
        'stick_HMinus_9_18',             # H- wall loss
        'stick_CH3Minus_9_15',           # CH3- wall loss
    ]
    selected |= set(critical_ion_loss)

    # Filter to only rates that exist in database
    selected = {name for name in selected if name in db}

    # Sort by range (largest first) for reporting
    selected_with_range = [
        (name, db[name].max / db[name].min if db[name].min > 0 else 1.0)
        for name in selected
    ]
    selected_with_range.sort(key=lambda x: x[1], reverse=True)

    selected_names = [name for name, _ in selected_with_range]

    print(f"  V3 expanded selection: {len(selected_names)} rates")
    print(f"    - Wall sticking: {len(wall_sticking)}")
    print(f"    - Volume loss: {len(volume_loss)}")
    print(f"    - Electron-impact (no ionization): {len(electron_impact)}")
    print(f"    - Range > 2x: {len(large_range)}")

    return selected_names


def analyze_chemistry(y_final, species, params):
    """Detailed chemistry analysis for important species."""

    # Get reaction rates
    rate_constants = np.array([params['k'][tag] for tag in params['tags']])
    reactions = params['R']

    reaction_rates = rate_constants.copy()
    for rxn_idx, reaction in enumerate(reactions):
        react_species = np.where(reaction.reactants > 0)[0]
        for sp_idx in react_species:
            coeff = reaction.reactants[sp_idx]
            reaction_rates[rxn_idx] *= y_final[sp_idx] ** coeff

    # Analyze important species
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


def run_simulation_with_logging(rate_values, E_field, ne, params_base, log_file=None):
    """Run simulation and optionally log detailed results."""

    try:
        params = params_base.copy()
        params['E_field'] = E_field
        params['ne'] = ne

        k = define_rates(params)
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

        ode_func = PlasmaODE_Optimized(params)

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

        # Get all species densities for charge balance calculation
        all_densities = {species[i]: float(y_final[i]) for i in range(ns)}

        results = {
            'H': get_density('H'),
            'CH': get_density('CH'),
            'C2': get_density('C2'),
            'all_densities': all_densities  # Include all for charge balance
        }

        # If logging requested, save detailed analysis
        if log_file:
            chemistry = analyze_chemistry(y_final, species, params)

            log_data = {
                'Ne': ne,
                'E_field': E_field,
                'rate_values': {k: float(v) for k, v in rate_values.items()},
                'target_densities': {
                    'H': results['H'],
                    'CH': results['CH'],
                    'C2': results['C2']
                },
                'all_densities': all_densities,
                'chemistry_analysis': chemistry
            }

            with open(log_file, 'w') as f:
                json.dump(log_data, f, indent=2)

        return results

    except Exception as e:
        return None


def objective_function(x, param_names, params_base):
    """Objective function with comprehensive logging."""
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

    # Extract target species errors
    weights = {
        'H': 1.0,
        'CH': 20.0,
        'C2': 3.0
    }

    error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        error += weights[species] * rel_error ** 2

    # ADD CHARGE BALANCE PENALTY
    all_dens = results['all_densities']

    # Calculate total positive ions
    positive_ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                     'H3Plus', 'CHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                     'C2H3Plus', 'C2HPlus']
    total_positive = sum(all_dens.get(ion, 0) for ion in positive_ions)

    # Ion/electron ratio (target: 4.5x, range 3-6x)
    ion_e_ratio = total_positive / ne if ne > 0 else 0

    # Charge balance error (penalize deviation from target ratio)
    charge_error = abs(ion_e_ratio - TARGET_ION_E_RATIO) / TARGET_ION_E_RATIO

    # Add weighted charge penalty to objective
    error += CHARGE_BALANCE_WEIGHT * charge_error

    # Save best result with detailed logging
    if error < best_result['objective']:
        best_result['objective'] = error
        best_result['params'] = {
            'Ne': ne,
            'E_field': E_field,
            'rates': dict(rate_values)
        }
        best_result['densities'] = results

        # Save detailed log for best result
        log_file = f'optimization_results_charge_balanced_v3/best_iteration_{iteration_counter:04d}_f{error:.1f}.json'
        run_simulation_with_logging(rate_values, E_field, ne, params_updated, log_file)

        print(f"\n  *** NEW BEST: f(x) = {error:.2f} at evaluation {objective_function.counter}")
        print(f"      H: {results['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {results['H']/TARGETS['H']:.2f}x)")
        print(f"      CH: {results['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {results['CH']/TARGETS['CH']:.2f}x)")
        print(f"      C2: {results['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {results['C2']/TARGETS['C2']:.2f}x)")
        print(f"      Ion/e ratio: {ion_e_ratio:.2f}x (target: {TARGET_ION_E_RATIO:.1f}x)")
        print(f"      Saved to: {log_file}\n")

    return error


def create_warm_start_population(param_names, bounds, db, popsize=6):
    """
    Create initial population biased for charge balance.

    Strategy for charge balance:
    - E ~ 10-50 V/cm (LOWER to retain ions, not lose them to walls!)
    - Ne ~ 1.5-2e9 (middle range)
    - Ion wall sticking ~ minimum (retain ions)
    - Key neutral rates toward minimums: C2H2+H→C2, C2+H→CH
    """
    print("\n  Creating warm-start population biased for CHARGE BALANCE...")

    population = np.zeros((popsize, len(bounds)))

    # Key reactions to bias toward minimum (to reduce C2H2→C2→CH pathway)
    bias_to_min = {
        'C2H2_H_C2_H2_H_cm3_5_4',  # C2H2 + H → C2 (91% of C2 production)
        'C2_H_CH_C_cm3_5_4',        # C2 + H → CH (59% of CH production)
    }

    # Ion wall sticking - bias toward MINIMUM to retain ions
    ion_wall_sticking = {
        'stick_ArPlus_9_4',
        'stick_HMinus_9_18',
        'stick_CH3Minus_9_15',
    }

    for i in range(popsize):
        for j, (name, (low, high)) in enumerate(zip(param_names, bounds)):
            if name == 'E_field':
                # Bias toward LOWER E (10-50 V/cm) to retain ions!
                population[i, j] = np.random.uniform(max(10, low), min(50, high))
            elif name == 'ne':
                # Bias toward lower Ne (1.3e9 - 2.0e9)
                population[i, j] = np.random.uniform(max(1.3e9, low), 2.0e9)
            elif name in bias_to_min:
                # Bias toward minimum for key reactions
                range_width = high - low
                population[i, j] = np.random.uniform(low, low + 0.3 * range_width)
            elif name in ion_wall_sticking:
                # Bias ion wall sticking toward MINIMUM to retain ions!
                range_width = high - low
                population[i, j] = np.random.uniform(low, low + 0.2 * range_width)
            else:
                # Random within bounds for other rates
                population[i, j] = np.random.uniform(low, high)

    print("  ✓ Population biased: E~10-50 V/cm (LOW!), Ne~1.5e9, ion sticking minimal")
    return population


def main():
    global iteration_counter

    print("=" * 80)
    print(" CHARGE-BALANCED OPTIMIZATION V3 - MASSIVELY EXPANDED RATES")
    print("=" * 80)
    print("\nV3 STRATEGY (per user suggestion - tune MORE rates!):")
    print("  ✓ 100+ tunable rates (was 40 in V2)")
    print("  ✓ ALL wall sticking rates (28 rates)")
    print("  ✓ ALL volume loss rates (25 rates)")
    print("  ✓ Most electron-impact (no ionization, ~26 rates)")
    print("  ✓ All rates with range > 2x (~77 rates)")
    print("  ✓ Charge penalty: 200x (strong constraint)")
    print("  ✓ Ionization rates FIXED (per user guidance)")
    print("\nV2 result: Ion/e = 4.58x ✓, but CH = 10.42x (still high)")
    print("V3 goal: More degrees of freedom → satisfy BOTH constraints!")

    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()
    print(f"  Selected {len(tunable_rates)} key rates")

    db = get_complete_rate_database()

    print("\nSelected rates:")
    print("-" * 80)
    for i, name in enumerate(tunable_rates[:10], 1):
        if name in db:
            rate = db[name]
            range_factor = rate.max / rate.min if rate.min > 0 else 1.0
            print(f"{i:2d}. {name:45s} [{rate.min:.2e}, {rate.max:.2e}] ({range_factor:.1f}x)")
    print(f"    ... and {len(tunable_rates) - 10} more")

    params_base = {
        'E_field': 50,
        'L_discharge': 0.45,
        'ne': 3.3e9,
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

    bounds.append((5.0, 300.0))  # E field (lower min to reduce ion loss)
    param_names.append('E_field')

    bounds.append((1.0e9, 5.0e9))  # Ne
    param_names.append('ne')

    print(f"\n Optimization parameters:")
    print(f"  Tunable rates: {len(param_names) - 2} (including ion reactions)")
    print(f"  E field: [5, 300] V/cm (lower min to retain ions)")
    print(f"  Ne: [1.0e9, 5.0e9] cm^-3")
    print(f"  Total parameters: {len(param_names)}")

    print(f"\n Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm^-3")
    print(f"  CH: {TARGETS['CH']:.2e} cm^-3  (heavily weighted)")
    print(f"  C2: {TARGETS['C2']:.2e} cm^-3")
    print(f"  Ion/e ratio: {TARGET_ION_E_RATIO:.1f}x (range: 3-6x for sheath region)")

    print("\n" + "=" * 80)
    print(" Running Charge-Balanced Optimization (15 iterations, pop=6)")
    print("=" * 80)
    print("\nDetailed logs will be saved to optimization_results_charge_balanced_v3/")
    print("Best results will be logged with full chemistry breakdown.\n")

    # Create warm-start population
    init_pop = create_warm_start_population(param_names, bounds, db, popsize=6)

    start_time = time.time()

    result = differential_evolution(
        objective_function,
        bounds,
        args=(param_names, params_base),
        strategy='best1bin',
        maxiter=15,   # Reasonable number to complete
        popsize=6,    # Smaller for speed
        tol=0.01,
        atol=0.0,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=42,
        init=init_pop,  # WARM START from previous best
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

    # Save final results
    final_params = {
        'Ne': result.x[-1],
        'E_field': result.x[-2],
        'rates': {name: result.x[i] for i, name in enumerate(param_names[:-2])}
    }

    with open('optimization_results_charge_balanced_v3/FINAL_RESULT.json', 'w') as f:
        json.dump({
            'objective': result.fun,
            'parameters': final_params,
            'success': result.success,
            'message': result.message,
            'iterations': result.nit,
            'evaluations': result.nfev
        }, f, indent=2)

    print(f"\n✓ Final parameters saved to optimization_results_charge_balanced_v3/FINAL_RESULT.json")
    print(f"✓ Best result logs saved to optimization_results_charge_balanced_v3/best_iteration_*.json")
    print("\nUse these logs to analyze optimal parameter combinations.")


if __name__ == '__main__':
    main()
