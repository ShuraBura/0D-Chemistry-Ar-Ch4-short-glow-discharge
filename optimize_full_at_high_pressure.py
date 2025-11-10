#!/usr/bin/env python3
"""
COMPREHENSIVE OPTIMIZATION AT HIGH PRESSURE (500 mTorr)

Combines the best insights from all previous analyses:
1. High pressure (500 mTorr) - gives 51% improvement baseline
2. Temperature-dependent rates (Te optimization)
3. ~40 key tunable rates within literature bounds

Focus areas based on chemistry analysis:
- C2H2 production (CH3+CH3→C2H2 dominates at 64.9%)
- C2H2→C2 pathway (83.8% of C2 production)
- C2+H→CH pathway (72.6% of CH production)
- CH loss pathways

Parameters: ~43 total
- Te: 0.3-7 eV
- Ne: 1e8-5e9 cm⁻³
- E-field: 10-500 V/cm
- ~40 reaction rates within literature bounds
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

# High pressure setting (from pressure scan finding)
PRESSURE_MTORR = 500.0

# Create results directory
os.makedirs('optimization_results_full_high_P', exist_ok=True)

# Track best result
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


def pressure_to_density(pressure_mTorr, T_K=400):
    """Convert pressure to number density."""
    kB = 1.38064852e-23  # J/K
    Torr_to_Pa = 133.322
    P_Pa = pressure_mTorr * 1e-3 * Torr_to_Pa
    n_m3 = P_Pa / (kB * T_K)
    n_cm3 = n_m3 * 1e-6
    return n_cm3


def select_tunable_rates():
    """
    Select ~40 most important tunable rates based on chemistry analysis.

    Priority:
    1. C2H2 production (affects entire cascade)
    2. C2H2 → C2 conversion
    3. C2 + H → CH (dominant CH source)
    4. CH loss pathways
    5. Rates with large uncertainty ranges
    """
    db = get_complete_rate_database()

    # Critical reactions identified in chemistry analysis
    critical_reactions = [
        # C2H2 production (top sources)
        'CH3_CH3_C2H2_H2_H2_cm3_7_49',      # 64.9% of C2H2 production!
        'C_CH3_C2_H2_H_cm3_7_8',            # 9.0%
        'CH3_C2H5_C2H2_CH3_H2_cm3_7_61',   # 6.6%
        'CH2_CH3_C2H2_H_H2_cm3_7_62',      # 6.5%
        'CH2_C2H3_C2H2_CH3_cm3_7_53',      # 5.4%
        'CH2_CH2_C2H2_H2_cm3_7_15',        # 4.5%

        # C2 production from C2H2
        'C2H2_H_C2_H2_H_cm3_7_50',         # 83.8% of C2 production!
        'C2H2_C_C2_CH2_cm3_7_19',          # 12.0%

        # CH production (C2 is main source)
        'C2_H_CH_C_cm3_7_6',               # 72.6% of CH production!
        'C_H_CH_cm3_7_12',                 # 13.7%
        'CH2_H_CH_H2_cm3_7_1',             # 12.4%

        # CH loss (need to increase these)
        'CH_CH4_C2H4_H_cm3_7_20',          # 71.6% of CH loss
        'CH_CH4_CH2_CH3_cm3_7_39',         # 4.8%
        'CH_H_C_H2_cm3_7_3',               # 1.4%
        'CH_CH3_C2H4_cm3_7_5',             # 2.1%

        # CH3 production/loss (affects C2H2)
        'e_CH4_CH3_H_cm3_1_1',             # CH3 from CH4
        'CH3_CH3_C2H2_H2_H2_cm3_7_49',    # CH3 sink (already listed)
        'CH3_H_CH4_cm3_7_35',              # CH3 sink

        # CH2 chemistry (CH precursor)
        'CH2_H_CH_H2_cm3_7_1',             # Already listed
        'CH2_CH2_C2H2_H2_cm3_7_15',       # Already listed
        'e_CH4_CH2_H2_cm3_1_2',           # CH2 from CH4
    ]

    # Get rates that affect our targets
    h_rates = set(get_tunable_rates_for_target('H').keys())
    ch_rates = set(get_tunable_rates_for_target('CH').keys())
    c2_rates = set(get_tunable_rates_for_target('C2').keys())

    target_rates = h_rates | ch_rates | c2_rates

    # Rates with large uncertainty (>3× range)
    large_range_rates = {
        name for name, rate in db.items()
        if (rate.max / rate.min > 3.0) and (rate.min > 0)
    }

    # Flagged rates from database
    flagged_rates = {name for name, rate in db.items() if rate.flag}

    # Combine all selections
    selected = set(critical_reactions) | target_rates | large_range_rates | flagged_rates
    selected = {name for name in selected if name in db}

    # Sort by range (larger range = more tunable)
    selected_with_range = [
        (name, db[name].max / db[name].min if db[name].min > 0 else 1.0)
        for name in selected
    ]
    selected_with_range.sort(key=lambda x: x[1], reverse=True)

    # Take top 40 rates
    selected_names = [name for name, _ in selected_with_range[:40]]

    return selected_names


def run_simulation(rate_values, Te, ne, E_field, params_base, log_file=None):
    """Run simulation with given parameters at high pressure."""

    try:
        # Calculate neutral density at high pressure
        n_total = pressure_to_density(PRESSURE_MTORR)

        params = params_base.copy()
        params['E_field'] = E_field
        params['ne'] = ne
        params['Te'] = Te

        # Get temperature-dependent rates
        k = define_rates(params)
        db = get_complete_rate_database()

        # Apply tuned rate values
        for name, val in rate_values.items():
            if name in k and name in db:
                val = np.clip(val, db[name].min, db[name].max)
                k[name] = val

        # Ensure all rates within bounds
        for name, rate_db in db.items():
            if name in k:
                if k[name] < rate_db.min:
                    k[name] = rate_db.min
                elif k[name] > rate_db.max:
                    k[name] = rate_db.max

        params['k'] = k
        params['R'], params['tags'] = build_reactions(params)

        # Initial conditions
        species = params['species']
        ns = len(species)
        y0 = np.ones(ns) * 1e3

        def set_density(name, value):
            try:
                idx = species.index(name)
                y0[idx] = value
            except ValueError:
                pass

        # Set neutral densities based on high pressure
        n_Ar = 0.85 * n_total
        n_CH4 = 0.15 * n_total

        set_density('e', ne)
        set_density('Ar', n_Ar)
        set_density('CH4', n_CH4)
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

        # Run simulation with timeout
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
            'C2': get_density('C2'),
            'C2H2': get_density('C2H2'),
        }

        # Save detailed log if requested
        if log_file:
            all_densities = {species[i]: float(y_final[i]) for i in range(ns)}
            log_data = {
                'pressure_mTorr': PRESSURE_MTORR,
                'n_total': float(n_total),
                'Te': Te,
                'Ne': ne,
                'E_field': E_field,
                'rate_values': {k: float(v) for k, v in rate_values.items()},
                'target_densities': {k: float(results[k]) for k in ['H', 'CH', 'C2']},
                'all_densities': all_densities,
            }

            with open(log_file, 'w') as f:
                json.dump(log_data, f, indent=2)

        return results

    except Exception as e:
        return None


def objective_function(x, param_names, params_base):
    """Objective function with emphasis on CH (biggest problem)."""
    global best_result

    if not hasattr(objective_function, 'counter'):
        objective_function.counter = 0
        objective_function.start_time = time.time()

    objective_function.counter += 1

    if objective_function.counter % 20 == 0:
        elapsed = time.time() - objective_function.start_time
        print(f"  [{objective_function.counter} evaluations, {elapsed/60:.1f} min elapsed]")

    # Extract parameters: rates, then Te, Ne, E
    rate_values = {name: val for name, val in zip(param_names[:-3], x[:-3])}
    Te = x[-3]
    ne = x[-2]
    E_field = x[-1]

    # Run simulation
    results = run_simulation(rate_values, Te, ne, E_field, params_base)

    if results is None:
        return 1e10

    # Weighted error - CH still heavily weighted
    weights = {
        'H': 1.0,
        'CH': 20.0,  # Still the biggest problem
        'C2': 3.0
    }

    error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        error += weights[species] * rel_error ** 2

    # Track and save best result
    if error < best_result['objective']:
        best_result['objective'] = error
        best_result['params'] = {
            'Te': Te,
            'Ne': ne,
            'E_field': E_field,
            'pressure_mTorr': PRESSURE_MTORR,
            'rates': dict(rate_values)
        }
        best_result['densities'] = results

        # Save detailed log
        log_file = f'optimization_results_full_high_P/best_f{error:.1f}_P{PRESSURE_MTORR:.0f}.json'
        run_simulation(rate_values, Te, ne, E_field, params_base, log_file)

        print(f"\n  *** NEW BEST: f(x) = {error:.2f} at evaluation {objective_function.counter}")
        print(f"      P: {PRESSURE_MTORR:.0f} mTorr, Te: {Te:.2f} eV, Ne: {ne:.2e}, E: {E_field:.1f} V/cm")
        print(f"      H: {results['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {results['H']/TARGETS['H']:.2f}x)")
        print(f"      CH: {results['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {results['CH']/TARGETS['CH']:.2f}x)")
        print(f"      C2: {results['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {results['C2']/TARGETS['C2']:.2f}x)")
        print(f"      Saved to: {log_file}\n")

    return error


def main():
    print("=" * 80)
    print(" COMPREHENSIVE OPTIMIZATION AT HIGH PRESSURE")
    print("=" * 80)
    print(f"\nPressure: {PRESSURE_MTORR} mTorr (51% improvement baseline)")
    print("Optimizing: ~43 parameters total")
    print("  - Te: 0.3-7 eV")
    print("  - Ne: 1e8-5e9 cm⁻³")
    print("  - E-field: 10-500 V/cm")
    print("  - ~40 key reaction rates (within literature bounds)")

    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()
    print(f"  Selected {len(tunable_rates)} key rates")

    db = get_complete_rate_database()

    print("\nTop 10 tunable rates:")
    print("-" * 80)
    for i, name in enumerate(tunable_rates[:10], 1):
        if name in db:
            rate = db[name]
            range_factor = rate.max / rate.min if rate.min > 0 else 1.0
            print(f"{i:2d}. {name:45s} range: {range_factor:.1f}×")
    print(f"    ... and {len(tunable_rates) - 10} more")

    # Base parameters
    params_base = {
        'E_field': 50,
        'L_discharge': 0.45,
        'ne': 3.3e9,
        'Te': 1.0,
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

    # Setup bounds
    bounds = []
    param_names = []

    # Reaction rate bounds
    for name in tunable_rates:
        if name in db:
            bounds.append((db[name].min, db[name].max))
            param_names.append(name)

    # Physical parameter bounds
    bounds.append((0.3, 7.0))      # Te (eV)
    param_names.append('Te')

    bounds.append((1e8, 5e9))      # Ne (cm⁻³)
    param_names.append('Ne')

    bounds.append((10.0, 500.0))   # E-field (V/cm)
    param_names.append('E_field')

    print(f"\n Optimization setup:")
    print(f"  Total parameters: {len(param_names)}")
    print(f"  Tunable rates: {len(tunable_rates)}")
    print(f"  Physical params: 3 (Te, Ne, E)")

    print(f"\n Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm⁻³")
    print(f"  CH: {TARGETS['CH']:.2e} cm⁻³  (currently 43× too high at 500 mTorr)")
    print(f"  C2: {TARGETS['C2']:.2e} cm⁻³  (currently 7× too high)")

    print("\n" + "=" * 80)
    print(" Running Optimization (20 iterations, pop=12)")
    print("=" * 80)
    print("\nThis combines:")
    print("  ✓ High pressure (500 mTorr)")
    print("  ✓ Temperature-dependent rates")
    print("  ✓ 40 tunable rates")
    print("  ✓ Te/Ne/E optimization")
    print("\nExpected: significant improvement over all previous optimizations!\n")

    start_time = time.time()

    result = differential_evolution(
        objective_function,
        bounds,
        args=(param_names, params_base),
        strategy='best1bin',
        maxiter=20,
        popsize=12,
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

    print("\n" + "=" * 80)
    print(f" OPTIMIZATION COMPLETE ({elapsed/60:.1f} minutes)")
    print("=" * 80)

    print(f"\nFinal objective: {result.fun:.6f}")
    print(f"Iterations: {result.nit}")
    print(f"Evaluations: {result.nfev}")

    # Extract optimized parameters
    Te_opt = result.x[-3]
    Ne_opt = result.x[-2]
    E_opt = result.x[-1]

    print(f"\nOptimized physical parameters:")
    print(f"  Pressure: {PRESSURE_MTORR} mTorr (fixed)")
    print(f"  Te: {Te_opt:.2f} eV")
    print(f"  Ne: {Ne_opt:.2e} cm⁻³")
    print(f"  E-field: {E_opt:.1f} V/cm")

    # Save final results
    final_data = {
        'objective': float(result.fun),
        'parameters': {
            'pressure_mTorr': float(PRESSURE_MTORR),
            'Te': float(Te_opt),
            'Ne': float(Ne_opt),
            'E_field': float(E_opt),
            'rates': {name: float(result.x[i]) for i, name in enumerate(param_names[:-3])}
        },
        'success': result.success,
        'message': result.message,
        'iterations': int(result.nit),
        'evaluations': int(result.nfev)
    }

    with open('optimization_results_full_high_P/FINAL_RESULT.json', 'w') as f:
        json.dump(final_data, f, indent=2)

    print(f"\n✓ Results saved to optimization_results_full_high_P/")
    print("\nRe-run simulation with these parameters to validate and analyze chemistry!")


if __name__ == '__main__':
    main()
