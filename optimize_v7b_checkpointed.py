#!/usr/bin/env python3
"""
V7b: BALANCED OPTIMIZATION WITH CHECKPOINTING
Balances species accuracy with charge balance (weight=150).
Survives system restarts by saving progress after each iteration batch.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import time
import signal
import os
import json
import pickle
from contextlib import contextmanager
from pathlib import Path

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

TARGET_ION_E_RATIO = 4.5
CHARGE_BALANCE_WEIGHT = 150.0  # NO CHARGE PENALTY

CHECKPOINT_FILE = 'v7b_checkpoint.pkl'
RESULTS_DIR = 'optimization_results_v7b_pure_species'
ITERATIONS_PER_BATCH = 1  # Save checkpoint after EVERY iteration (survive fast reboots)
TOTAL_ITERATIONS = 25

os.makedirs(RESULTS_DIR, exist_ok=True)

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

    selected |= {name for name, rate in db.items() if (rate.max / rate.min > 1.5) and (rate.min > 0)}
    selected |= {name for name in db.keys() if name.startswith('stick_')}
    selected |= {name for name in db.keys() if name.startswith('loss_')}

    ionization_keywords = ['CH4Plus', 'CH3Plus', 'ArPlus', 'CH5Plus', 'CHPlus', 'CH2Plus',
                          'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus']
    selected |= {name for name in db.keys() if name.startswith('e_') and
                 not any(ion in name for ion in ionization_keywords)}

    dissoc_recomb = [
        'e_CH4Plus_CH3_H_cm3_6_4', 'e_CH4Plus_CH2_H2_cm3_6_9',
        'e_CH4Plus_CH_H2_H_cm3_6_11', 'e_CH4Plus_C_2H2_cm3_6_13',
        'CH5Plus_e_CH4_H_cm3_6_3', 'CH5Plus_e_CH3_H2_cm3_6_8',
        'CH5Plus_e_CH2_H2_H_cm3_6_10', 'CH5Plus_e_CH3_2H_cm3_6_12',
        'ArPlus_e_Ar_cm3_6_1', 'CH3Plus_e_CH3_cm3_6_2',
        'C2H5Plus_e_C2H4_H_cm3_6_14', 'C2H4Plus_e_C2H2_H2_cm3_6_15',
        'C2H3Plus_e_C2H2_H_cm3_6_16', 'C2HPlus_e_C2_H_cm3_6_18',
        'C2H5Plus_e_C2H4_H_cm3_6_28',
        'ArHPlus_e_Ar_H_cm3_6_29',
        'e_CH4_CH3_HMinus_cm3_8_1',
    ]
    selected |= set(dissoc_recomb)

    secondary_tertiary = [
        'H_CH4_CH3_H2_cm3_7_25', 'CH3_CH3_C2H2_H2_H2_cm3_7_49',
        'CH2_CH3_C2H2_H_H2_cm3_7_62', 'C2H2_H_C2_H2_H_cm3_7_50',
        'CH2_H_CH_H2_cm3_7_1', 'CH_CH4_C2H4_H_cm3_7_20',
        'C2_H_CH_C_cm3_7_6', 'CH_H_C_H2_cm3_7_3',
    ]
    selected |= set(secondary_tertiary)

    selected |= {name for name, rate in db.items() if rate.flag}
    selected = {name for name in selected if name in db}

    return sorted(selected)

def run_simulation(rate_values, E_field, ne, params_base):
    try:
        params = params_base.copy()
        params['E_field'] = E_field
        params['ne'] = ne
        k = define_rates(params)
        db = get_complete_rate_database()

        for name, val in rate_values.items():
            if name in k and name in db:
                k[name] = np.clip(val, db[name].min, db[name].max)

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

        with time_limit(30):
            sol = solve_ivp(ode_func, (0, 100), y0, method='BDF',
                           rtol=1e-5, atol=1e-6, max_step=10.0)

        if not sol.success:
            return None

        y_final = sol.y[:, -1]
        all_densities = {species[i]: float(y_final[i]) for i in range(ns)}

        return {
            'H': all_densities.get('H', 0),
            'CH': all_densities.get('CH', 0),
            'C2': all_densities.get('C2', 0),
            'all_densities': all_densities
        }
    except (TimeoutException, Exception):
        return None

def objective_function(x, param_names, params_base, eval_counter):
    eval_counter['count'] += 1

    if eval_counter['count'] % 10 == 0:
        elapsed = time.time() - eval_counter['start_time']
        print(f"  [{eval_counter['count']} evaluations, {elapsed/60:.1f} min elapsed]")

    rate_values = {name: val for name, val in zip(param_names[:-2], x[:-2])}
    E_field = x[-2]
    ne = x[-1]

    results = run_simulation(rate_values, E_field, ne, params_base)

    if results is None:
        return 1e10

    weights = {'H': 1.0, 'CH': 30.0, 'C2': 3.0}
    error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        error += weights[species] * rel_error ** 2

    # ADD CHARGE BALANCE PENALTY
    all_dens = results['all_densities']
    positive_ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                     'H3Plus', 'CHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                     'C2H3Plus', 'C2HPlus']
    total_positive = sum(all_dens.get(ion, 0) for ion in positive_ions)
    ion_e_ratio = total_positive / ne if ne > 0 else 0
    charge_error = abs(ion_e_ratio - TARGET_ION_E_RATIO) / TARGET_ION_E_RATIO
    error += CHARGE_BALANCE_WEIGHT * charge_error

    if error < eval_counter['best_obj']:
        eval_counter['best_obj'] = error
        eval_counter['best_x'] = x.copy()
        eval_counter['best_result'] = results

        print(f"\n  *** NEW BEST: f(x) = {error:.2f} at evaluation {eval_counter['count']}")
        print(f"      H: {results['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {results['H']/TARGETS['H']:.2f}x)")
        print(f"      CH: {results['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {results['CH']/TARGETS['CH']:.2f}x)")
        print(f"      C2: {results['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {results['C2']/TARGETS['C2']:.2f}x)")
        print(f"      Ion/e: {ion_e_ratio:.2f}x (target: {TARGET_ION_E_RATIO}x, weight: {CHARGE_BALANCE_WEIGHT})\n")

    return error

def load_checkpoint():
    """Load checkpoint if it exists"""
    if Path(CHECKPOINT_FILE).exists():
        print(f"\n{'='*80}")
        print(f" CHECKPOINT FOUND - Resuming from {CHECKPOINT_FILE}")
        print(f"{'='*80}\n")
        with open(CHECKPOINT_FILE, 'rb') as f:
            checkpoint = pickle.load(f)
        return checkpoint
    return None

def save_checkpoint(checkpoint):
    """Save checkpoint to disk"""
    with open(CHECKPOINT_FILE, 'wb') as f:
        pickle.dump(checkpoint, f)
    print(f"\n  → Checkpoint saved to {CHECKPOINT_FILE}\n")

def main():
    print("=" * 80)
    print(" V7b OPTIMIZATION WITH CHECKPOINTING - BALANCED APPROACH")
    print("=" * 80)
    print("\nSTRATEGY:")
    print("  ✓ BALANCED: Species targets + charge balance (WEIGHT = 150)")
    print("  ✓ Target Ion/e ratio: 4.5x")
    print("  ✓ Checkpointing: survives system restarts")
    print(f"  ✓ {TOTAL_ITERATIONS} iterations, {ITERATIONS_PER_BATCH} per batch")

    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()
    print(f"  V7b selection: {len(tunable_rates)} rates")

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

    # Load or create checkpoint
    checkpoint = load_checkpoint()
    if checkpoint:
        completed_iters = checkpoint['completed_iterations']
        best_x = checkpoint['best_x']
        best_obj = checkpoint['best_obj']
        total_evals = checkpoint['total_evaluations']
        print(f"Resuming from iteration {completed_iters}/{TOTAL_ITERATIONS}")
        print(f"Best objective so far: {best_obj:.2f}")
        print(f"Total evaluations: {total_evals}")
    else:
        completed_iters = 0
        best_x = None
        best_obj = 1e10
        total_evals = 0
        print(f"\nStarting fresh optimization: {TOTAL_ITERATIONS} iterations")

    overall_start = time.time()

    # Run optimization in batches
    while completed_iters < TOTAL_ITERATIONS:
        remaining_iters = TOTAL_ITERATIONS - completed_iters
        batch_iters = min(ITERATIONS_PER_BATCH, remaining_iters)

        print("\n" + "=" * 80)
        print(f" BATCH: Iterations {completed_iters+1}-{completed_iters+batch_iters} of {TOTAL_ITERATIONS}")
        print("=" * 80)

        eval_counter = {
            'count': 0,
            'start_time': time.time(),
            'best_obj': best_obj,
            'best_x': best_x,
            'best_result': None
        }

        # Run differential evolution for this batch
        result = differential_evolution(
            objective_function,
            bounds,
            args=(param_names, params_base, eval_counter),
            strategy='best1bin',
            maxiter=batch_iters,
            popsize=6,
            tol=0.01,
            atol=0.0,
            mutation=(0.5, 1.0),
            recombination=0.7,
            seed=42 + completed_iters,  # Different seed for each batch
            disp=True,
            workers=1,
            updating='deferred',
            polish=(batch_iters == remaining_iters),  # Polish only on last batch
            x0=best_x  # Warm start from best solution
        )

        # Update checkpoint
        completed_iters += batch_iters
        total_evals += eval_counter['count']
        best_x = eval_counter['best_x'] if eval_counter['best_x'] is not None else result.x
        best_obj = eval_counter['best_obj']

        checkpoint = {
            'completed_iterations': completed_iters,
            'total_iterations': TOTAL_ITERATIONS,
            'best_x': best_x,
            'best_obj': best_obj,
            'total_evaluations': total_evals,
            'param_names': param_names,
            'timestamp': time.time()
        }
        save_checkpoint(checkpoint)

        elapsed = time.time() - overall_start
        print(f"\n→ Progress: {completed_iters}/{TOTAL_ITERATIONS} iterations complete")
        print(f"→ Best objective: {best_obj:.2f}")
        print(f"→ Total time: {elapsed/60:.1f} minutes")

    # Final results
    print("\n" + "=" * 80)
    print(" V7b OPTIMIZATION COMPLETE!")
    print("=" * 80)

    elapsed_total = time.time() - overall_start
    print(f"\nTotal runtime: {elapsed_total/60:.1f} minutes")
    print(f"Total evaluations: {total_evals}")
    print(f"Final objective: {best_obj:.6f}")

    # Save final result
    with open(f'{RESULTS_DIR}/FINAL_RESULT.json', 'w') as f:
        json.dump({
            'objective': float(best_obj),
            'parameters': {
                'Ne': float(best_x[-1]),
                'E_field': float(best_x[-2]),
                'rates': {name: float(best_x[i]) for i, name in enumerate(param_names[:-2])}
            },
            'iterations': TOTAL_ITERATIONS,
            'evaluations': total_evals,
            'runtime_minutes': elapsed_total/60
        }, f, indent=2)

    # Clean up checkpoint file
    if Path(CHECKPOINT_FILE).exists():
        os.remove(CHECKPOINT_FILE)
        print(f"\n→ Checkpoint file removed (optimization complete)")

if __name__ == '__main__':
    main()
