#!/usr/bin/env python3
"""
FINAL optimization to nail target densities within ±30%.

Starting from best solution (f=684):
- H: 2.52e13 (need 5.18e13, 0.49x) → Need +2x
- CH: 6.84e9 (need 1.00e9, 6.84x) → Need -7x
- C2: 5.11e10 (need 1.30e11, 0.39x) → Need +2.5x

Strategy:
- Reduce C2 removal to build up C2 (relax wall sticking/volume loss)
- Minimize C2+H→CH to prevent CH increase from more C2
- Reduce H wall sticking to build up H
- Keep CH removal high

Extended optimization (30 iterations, pop=8) for convergence.
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

# Create results directory
os.makedirs('optimization_results_final', exist_ok=True)

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
    """Select ~30 most impactful rates for optimization."""
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

    selected = target_rates | large_range_rates | flagged_rates
    selected = {name for name in selected if name in db}

    selected_with_range = [
        (name, db[name].max / db[name].min if db[name].min > 0 else 1.0)
        for name in selected
    ]
    selected_with_range.sort(key=lambda x: x[1], reverse=True)

    selected_names = [name for name, _ in selected_with_range[:30]]

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

        results = {
            'H': get_density('H'),
            'CH': get_density('CH'),
            'C2': get_density('C2')
        }

        # If logging requested, save detailed analysis
        if log_file:
            all_densities = {species[i]: float(y_final[i]) for i in range(ns)}
            chemistry = analyze_chemistry(y_final, species, params)

            log_data = {
                'Ne': ne,
                'E_field': E_field,
                'rate_values': {k: float(v) for k, v in rate_values.items()},
                'target_densities': results,
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

    weights = {
        'H': 1.0,
        'CH': 20.0,
        'C2': 3.0
    }

    error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        error += weights[species] * rel_error ** 2

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
        log_file = f'optimization_results_final/best_iteration_{iteration_counter:04d}_f{error:.1f}.json'
        run_simulation_with_logging(rate_values, E_field, ne, params_updated, log_file)

        print(f"\n  *** NEW BEST: f(x) = {error:.2f} at evaluation {objective_function.counter}")
        print(f"      H: {results['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {results['H']/TARGETS['H']:.2f}x)")
        print(f"      CH: {results['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {results['CH']/TARGETS['CH']:.2f}x)")
        print(f"      C2: {results['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {results['C2']/TARGETS['C2']:.2f}x)")
        print(f"      Saved to: {log_file}\n")

    return error


def create_warm_start_population(param_names, bounds, db, popsize=8):
    """
    Create initial population from previous best (f=684) with strategic adjustments.

    Previous best: Ne=1.28e9, E=300, f(x)=684
    - H: 0.49x (need +2x) → reduce H wall sticking
    - CH: 6.84x (need -7x) → minimize C2→CH, keep CH removal high
    - C2: 0.39x (need +2.5x) → reduce C2 wall sticking/volume loss
    """
    print("\n  Creating warm-start population from f=684 solution...")

    # Load previous best parameters
    with open('optimization_results/FINAL_RESULT.json', 'r') as f:
        prev_best = json.load(f)['parameters']

    population = np.zeros((popsize, len(bounds)))

    # Strategy-based rate adjustments for target balance
    reduce_to_min = {
        'C2_H_CH_C_cm3_5_4',       # C2+H→CH: minimize (prevent CH from C2)
        'stick_H_9_1',              # H sticking: reduce (build up H)
        'stick_C2_9_9',             # C2 sticking: reduce (build up C2, was at MAX)
        'loss_C2_11_3',             # C2 volume loss: reduce (build up C2, was at MAX)
        'C2H2_H_C2_H2_H_cm3_5_4',  # C2H2→C2: reduce slightly (less C2 feeding)
    }

    keep_high = {
        'stick_CH_9_3',             # CH sticking: keep high (remove CH)
        'loss_CH_11_9',             # CH volume loss: keep high (remove CH)
        'CH2_H_CH_H2_cm3_5_4',     # CH2→CH: keep moderate (major CH source)
    }

    for i in range(popsize):
        for j, (name, (low, high)) in enumerate(zip(param_names, bounds)):
            if name == 'E_field':
                # Keep E high (was at MAX 300)
                if i == 0:
                    population[i, j] = prev_best['E_field']
                else:
                    population[i, j] = np.random.uniform(270, 300)

            elif name == 'ne':
                # Keep Ne low but explore slightly higher for H balance
                if i == 0:
                    population[i, j] = prev_best['Ne']
                else:
                    population[i, j] = np.random.uniform(1.2e9, 1.8e9)

            elif name in reduce_to_min:
                # Bias toward minimum for these critical rates
                range_width = high - low
                if i == 0 and name in prev_best['rates']:
                    # First member: start from previous, push toward min
                    prev_val = prev_best['rates'][name]
                    population[i, j] = max(low, prev_val - 0.3 * range_width)
                else:
                    population[i, j] = np.random.uniform(low, low + 0.2 * range_width)

            elif name in keep_high:
                # Bias toward maximum
                range_width = high - low
                if i == 0 and name in prev_best['rates']:
                    population[i, j] = prev_best['rates'][name]
                else:
                    population[i, j] = np.random.uniform(high - 0.3 * range_width, high)

            elif i == 0 and name in prev_best['rates']:
                # First member: start from previous best with small perturbation
                prev_val = prev_best['rates'][name]
                noise = 0.1 * (high - low)
                population[i, j] = np.clip(prev_val + np.random.uniform(-noise, noise), low, high)

            else:
                # Others: random within bounds
                population[i, j] = np.random.uniform(low, high)

    print("  ✓ Strategy:")
    print("    - Reduce: C2 removal, H removal, C2→CH rate")
    print("    - Keep high: CH removal")
    print("    - Start from f=684 with adjustments for H/C2/CH balance")

    return population


def main():
    global iteration_counter

    print("=" * 80)
    print(" FINAL OPTIMIZATION: TARGET DENSITIES WITHIN ±30%")
    print("=" * 80)
    print("\nStarting from best solution: f(x) = 684")
    print("\nCurrent status:")
    print("  H:  0.49x target → Need +2x (reduce H loss)")
    print("  CH: 6.84x target → Need -7x (minimize C2→CH, keep CH loss high)")
    print("  C2: 0.39x target → Need +2.5x (reduce C2 loss)")
    print("\nStrategy: Reduce C2/H removal, minimize C2→CH, keep CH removal high")

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

    bounds.append((10.0, 300.0))  # E field
    param_names.append('E_field')

    bounds.append((1.0e9, 5.0e9))  # Ne
    param_names.append('ne')

    print(f"\n Optimization parameters:")
    print(f"  Tunable rates: {len(param_names) - 2}")
    print(f"  E field: [10, 300] V/cm")
    print(f"  Ne: [1.0e9, 5.0e9] cm^-3")
    print(f"  Total parameters: {len(param_names)}")

    print(f"\n Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm^-3")
    print(f"  CH: {TARGETS['CH']:.2e} cm^-3  (heavily weighted)")
    print(f"  C2: {TARGETS['C2']:.2e} cm^-3")

    print("\n" + "=" * 80)
    print(" Running Optimization (30 iterations, pop=8)")
    print("=" * 80)
    print("\nGoal: All targets within ±30%")
    print("  H:  [3.63e13, 6.73e13] cm^-3")
    print("  CH: [7.0e8, 1.3e9] cm^-3")
    print("  C2: [9.1e10, 1.69e11] cm^-3")
    print("\nDetailed logs will be saved to optimization_results_final/")
    print("Best results will be logged with full chemistry breakdown.\n")

    # Create warm-start population
    init_pop = create_warm_start_population(param_names, bounds, db, popsize=8)

    start_time = time.time()

    result = differential_evolution(
        objective_function,
        bounds,
        args=(param_names, params_base),
        strategy='best1bin',
        maxiter=30,   # Extended for convergence
        popsize=8,    # Larger for exploration
        tol=0.01,
        atol=0.0,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=43,      # Different seed for new optimization
        init=init_pop,  # WARM START from f=684
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

    with open('optimization_results_final/FINAL_RESULT.json', 'w') as f:
        json.dump({
            'objective': result.fun,
            'parameters': final_params,
            'success': result.success,
            'message': result.message,
            'iterations': result.nit,
            'evaluations': result.nfev
        }, f, indent=2)

    print(f"\n✓ Final parameters saved to optimization_results_final/FINAL_RESULT.json")
    print(f"✓ Best result logs saved to optimization_results_final/best_iteration_*.json")
    print("\nUse these logs to analyze optimal parameter combinations.")

    # Check if targets achieved
    final_results = run_simulation_with_logging(final_params['rates'], final_params['E_field'],
                                                 final_params['Ne'], params_base)
    if final_results:
        H_ratio = final_results['H'] / TARGETS['H']
        CH_ratio = final_results['CH'] / TARGETS['CH']
        C2_ratio = final_results['C2'] / TARGETS['C2']

        within_30pct = all([
            0.7 <= H_ratio <= 1.3,
            0.7 <= CH_ratio <= 1.3,
            0.7 <= C2_ratio <= 1.3
        ])

        if within_30pct:
            print("\n" + "=" * 80)
            print("  ✓✓✓ SUCCESS: ALL TARGETS WITHIN ±30%! ✓✓✓")
            print("=" * 80)
        else:
            print("\n  Closest to targets achieved - see logs for details")


if __name__ == '__main__':
    main()
