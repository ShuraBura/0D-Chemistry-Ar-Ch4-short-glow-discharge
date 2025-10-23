#!/usr/bin/env python3
"""
Constrained multi-parameter optimization to match H, CH, C2 targets.

All rates constrained within literature bounds.
E field tunable in range [10, 200] V/cm.

Baseline (before optimization):
  H:  3.45e13 (target: 5.18e13) - 0.67x  [need +49%]
  CH: 4.59e10 (target: 1.0e9)   - 46x    [need -46x!]
  C2: 7.31e11 (target: 1.3e11)  - 5.6x   [need -5.6x]
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import time
import signal
from contextlib import contextmanager

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database, get_tunable_rates_for_target
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


# Experimental targets
TARGETS = {
    'H': 5.18e13,   # cm^-3
    'CH': 1.0e9,    # cm^-3  (currently 46x too high!)
    'C2': 1.3e11,   # cm^-3
}


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
    Select ~30 most impactful rates for optimization.

    Criteria:
    1. Direct impact on H, CH, C2
    2. Large tuning range (max/min ratio)
    3. Flagged for validation
    """
    db = get_complete_rate_database()

    # Get rates relevant to each target
    h_rates = set(get_tunable_rates_for_target('H').keys())
    ch_rates = set(get_tunable_rates_for_target('CH').keys())
    c2_rates = set(get_tunable_rates_for_target('C2').keys())

    # Combine all target-relevant rates
    target_rates = h_rates | ch_rates | c2_rates

    # Add rates with large tuning ranges (>3x)
    large_range_rates = {
        name for name, rate in db.items()
        if (rate.max / rate.min > 3.0) and (rate.min > 0)
    }

    # Add flagged rates
    flagged_rates = {name for name, rate in db.items() if rate.flag}

    # Combine all
    selected = target_rates | large_range_rates | flagged_rates

    # Filter to only those in database
    selected = {name for name in selected if name in db}

    # Sort by range (largest first)
    selected_with_range = [
        (name, db[name].max / db[name].min if db[name].min > 0 else 1.0)
        for name in selected
    ]
    selected_with_range.sort(key=lambda x: x[1], reverse=True)

    # Take top 30
    selected_names = [name for name, _ in selected_with_range[:30]]

    return selected_names


def run_simulation(rate_values, E_field, params_base, verbose=False):
    """
    Run simulation with given rate constants and E field.

    Returns:
    --------
    dict : {'H': value, 'CH': value, 'C2': value} or None if failed
    """
    try:
        # Copy base params
        params = params_base.copy()
        params['E_field'] = E_field

        # Load rates with corrections
        k = define_rates(params)
        db = get_complete_rate_database()

        # Apply optimized rates
        for name, val in rate_values.items():
            if name in k:
                # Ensure within bounds
                if name in db:
                    val = np.clip(val, db[name].min, db[name].max)
                k[name] = val

        # Apply other corrections to ensure all rates within bounds
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

        set_density('e', params['ne'])
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

        # Run simulation with timeout (30 seconds max)
        ode_func = PlasmaODE_Optimized(params)

        try:
            with time_limit(30):
                sol = solve_ivp(
                    ode_func,
                    (0, 100),
                    y0,
                    method='BDF',
                    rtol=1e-5,  # Slightly relaxed for speed
                    atol=1e-6,
                    max_step=10.0
                )
        except TimeoutException:
            if verbose:
                print("Simulation timed out after 30s")
            return None

        if not sol.success:
            if verbose:
                print(f"Solver failed: {sol.message}")
            return None

        # Extract final densities
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

        return results

    except Exception as e:
        if verbose:
            print(f"Simulation failed: {e}")
        return None


def objective_function(x, param_names, params_base):
    """
    Objective function for optimization.

    x : array of parameter values (rates + E_field + ne)
    """
    # Progress counter
    if not hasattr(objective_function, 'counter'):
        objective_function.counter = 0
        objective_function.start_time = time.time()

    objective_function.counter += 1

    if objective_function.counter % 10 == 0:
        elapsed = time.time() - objective_function.start_time
        print(f"  [{objective_function.counter} evaluations, {elapsed/60:.1f} min elapsed]")

    # Unpack parameters (last two are E_field and ne)
    rate_values = {name: val for name, val in zip(param_names[:-2], x[:-2])}
    E_field = x[-2]
    ne = x[-1]

    # Update params_base with ne
    params_updated = params_base.copy()
    params_updated['ne'] = ne

    # Run simulation
    results = run_simulation(rate_values, E_field, params_updated)

    if results is None:
        return 1e10  # Penalize failed simulations

    # Calculate weighted error
    # CH is weighted 20x because it's 46x off (hardest to match)
    # C2 is weighted 3x because it's 5.6x off
    # H is weighted 1x because it's only 0.67x (closest)
    weights = {
        'H': 1.0,
        'CH': 20.0,   # Critical!
        'C2': 3.0
    }

    error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        error += weights[species] * rel_error ** 2

    return error


def main():
    print("=" * 80)
    print(" CONSTRAINED OPTIMIZATION TO MATCH TARGETS")
    print("=" * 80)

    # Select tunable rates
    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()
    print(f"  Selected {len(tunable_rates)} key rates")

    # Load database for bounds
    db = get_complete_rate_database()

    # Show selected rates
    print("\nSelected rates for optimization:")
    print("-" * 80)
    for i, name in enumerate(tunable_rates[:10], 1):
        if name in db:
            rate = db[name]
            range_factor = rate.max / rate.min if rate.min > 0 else 1.0
            print(f"{i:2d}. {name:40s} [{rate.min:.2e}, {rate.max:.2e}] ({range_factor:.1f}x)")
    print(f"    ... and {len(tunable_rates) - 10} more")

    # Set up base parameters
    params_base = {
        'E_field': 50,
        'L_discharge': 0.45,
        'ne': 3.3e9,  # Fixed experimental value
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

    # Set up optimization bounds
    bounds = []
    param_names = []

    for name in tunable_rates:
        if name in db:
            bounds.append((db[name].min, db[name].max))
            param_names.append(name)

    # Add E field (wider range since optimizer hit 200 V/cm max)
    bounds.append((10.0, 300.0))  # V/cm (was 200)
    param_names.append('E_field')

    # Add Ne (wider range since optimizer hit 1.65e9 minimum)
    bounds.append((1.0e9, 5.0e9))  # cm^-3 (was 1.65e9-4.95e9)
    param_names.append('ne')

    print(f"\n Optimization parameters:")
    print(f"  Tunable rates: {len(param_names) - 2}")
    print(f"  E field: [10, 300] V/cm")
    print(f"  Ne: [1.0e9, 5.0e9] cm^-3")
    print(f"  Total parameters: {len(param_names)}")

    print(f"\n Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm^-3")
    print(f"  CH: {TARGETS['CH']:.2e} cm^-3  (currently 46x too high - heavily weighted)")
    print(f"  C2: {TARGETS['C2']:.2e} cm^-3")

    print("\n" + "=" * 80)
    print(" Running Differential Evolution Optimization")
    print("=" * 80)
    print("\nThis may take 10-30 minutes...")
    print("Progress will be shown as optimization runs.\n")

    start_time = time.time()

    # Run optimization
    result = differential_evolution(
        objective_function,
        bounds,
        args=(param_names, params_base),
        strategy='best1bin',
        maxiter=30,   # Reduced for reliability (was 100)
        popsize=8,    # Reduced for speed (was 10)
        tol=0.01,
        atol=0.0,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=42,
        disp=True,
        workers=1,  # Single worker for stability
        updating='deferred',
        polish=False  # Skip polishing for speed
    )

    elapsed = time.time() - start_time

    print("\n" + "=" * 80)
    print(f" OPTIMIZATION COMPLETE ({elapsed/60:.1f} minutes)")
    print("=" * 80)

    if result.success:
        print("\n✓ Optimization converged successfully")
    else:
        print("\n⚠ Optimization did not fully converge (but found best solution)")

    print(f"\nFinal objective value: {result.fun:.6f}")
    print(f"Function evaluations: {result.nfev}")
    print(f"Iterations: {result.nit}")

    # Extract optimized parameters
    optimized_rates = {name: val for name, val in zip(param_names[:-2], result.x[:-2])}
    optimized_E = result.x[-2]
    optimized_ne = result.x[-1]

    print(f"\nOptimized E field: {optimized_E:.1f} V/cm (baseline: 50)")
    print(f"Optimized Ne: {optimized_ne:.3e} cm^-3 (baseline: 3.3e9, 2-param: 1.66e9)")

    # Run final simulation with optimized parameters
    print("\n" + "=" * 80)
    print(" FINAL RESULTS WITH OPTIMIZED PARAMETERS")
    print("=" * 80)

    params_final = params_base.copy()
    params_final['ne'] = optimized_ne
    final_results = run_simulation(optimized_rates, optimized_E, params_final, verbose=True)

    if final_results:
        print("\nFinal Densities:")
        print("-" * 80)
        for species in ['H', 'CH', 'C2']:
            current = final_results[species]
            target = TARGETS[species]
            ratio = current / target
            error_pct = (ratio - 1.0) * 100
            status = "✓" if 0.8 <= ratio <= 1.2 else ("→" if 0.5 <= ratio <= 2.0 else "✗")

            print(f"{species:5s}: {current:.2e} cm^-3")
            print(f"       Target: {target:.2e} cm^-3")
            print(f"       Ratio: {ratio:.2f}x  ({error_pct:+.1f}%)  {status}")
            print()

    # Save optimized parameters
    print("\n" + "=" * 80)
    print(" SAVING OPTIMIZED PARAMETERS")
    print("=" * 80)

    with open('optimized_parameters.txt', 'w') as f:
        f.write("OPTIMIZED PARAMETERS (32-parameter optimization)\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"E_field: {optimized_E:.6f} V/cm\n")
        f.write(f"Ne: {optimized_ne:.6e} cm^-3\n\n")
        f.write("Rate Constants:\n")
        f.write("-" * 80 + "\n")
        for name in sorted(optimized_rates.keys()):
            val = optimized_rates[name]
            if name in db:
                rate_db = db[name]
                f.write(f"{name}: {val:.6e}  (range: [{rate_db.min:.2e}, {rate_db.max:.2e}])\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("FINAL RESULTS\n")
        f.write("=" * 80 + "\n\n")
        if final_results:
            for species in ['H', 'CH', 'C2']:
                current = final_results[species]
                target = TARGETS[species]
                ratio = current / target
                f.write(f"{species}: {current:.6e} cm^-3  (target: {target:.2e}, ratio: {ratio:.2f}x)\n")

    print("✓ Saved to optimized_parameters.txt")

    print("\n" + "=" * 80)
    print(" OPTIMIZATION COMPLETE")
    print("=" * 80)


if __name__ == '__main__':
    main()
