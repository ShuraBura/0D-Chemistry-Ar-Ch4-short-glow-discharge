#!/usr/bin/env python3
"""
Simple 2-parameter optimization: Ne and E field only.
Keep all reaction rates at literature values.

This tests if targets are achievable by tuning just the physical parameters
before resorting to adjusting chemistry rates.

Ne range: 1.65e9 to 4.95e9 (3.3e9 ± 50%)
E field range: 10 to 200 V/cm
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import time

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


# Experimental targets
TARGETS = {
    'H': 5.18e13,   # cm^-3
    'CH': 1.0e9,    # cm^-3
    'C2': 1.3e11,   # cm^-3
}


def run_simulation(ne, E_field, verbose=False):
    """
    Run simulation with given Ne and E field.
    All rates kept at literature values.

    Returns:
    --------
    dict : {'H': value, 'CH': value, 'C2': value} or None if failed
    """
    try:
        # Setup parameters
        params = {
            'E_field': E_field,
            'L_discharge': 0.45,
            'ne': ne,
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

        # Load rates with literature corrections
        k = define_rates(params)
        db = get_complete_rate_database()

        # Ensure all rates within literature bounds
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

        # Run simulation
        ode_func = PlasmaODE_Optimized(params)
        sol = solve_ivp(
            ode_func,
            (0, 100),
            y0,
            method='BDF',
            rtol=1e-5,
            atol=1e-6,
            max_step=10.0
        )

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


def objective_function(x, verbose=False):
    """
    Objective function for 2-parameter optimization.

    x[0] = Ne (electron density)
    x[1] = E field
    """
    ne = x[0]
    E_field = x[1]

    # Run simulation
    results = run_simulation(ne, E_field, verbose=verbose)

    if results is None:
        return 1e10  # Penalize failed simulations

    # Calculate weighted error
    # CH weighted 20x (hardest to match, 46x off)
    # C2 weighted 3x (5.6x off)
    # H weighted 1x (0.67x, closest)
    weights = {
        'H': 1.0,
        'CH': 20.0,
        'C2': 3.0
    }

    error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        error += weights[species] * rel_error ** 2

    return error


def main():
    print("=" * 80)
    print(" 2-PARAMETER OPTIMIZATION: Ne and E field")
    print("=" * 80)
    print("\nStrategy: Tune only physical parameters, keep all rates at literature values")
    print("\nParameters:")
    print("  Ne (electron density): [1.65e9, 4.95e9] cm^-3 (3.3e9 ± 50%)")
    print("  E field: [10, 200] V/cm")
    print("\nAll reaction rates constrained to literature values.")

    print("\n Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm^-3")
    print(f"  CH: {TARGETS['CH']:.2e} cm^-3  (currently 46x too high)")
    print(f"  C2: {TARGETS['C2']:.2e} cm^-3")

    # Run baseline with original Ne=3.3e9, E=50
    print("\n" + "=" * 80)
    print(" BASELINE (Ne=3.3e9, E=50)")
    print("=" * 80)

    baseline = run_simulation(3.3e9, 50.0, verbose=True)
    if baseline:
        print("\nBaseline Densities:")
        for species in ['H', 'CH', 'C2']:
            current = baseline[species]
            target = TARGETS[species]
            ratio = current / target
            print(f"  {species:5s}: {current:.2e} cm^-3  (target: {target:.2e}, ratio: {ratio:.2f}x)")

    baseline_error = objective_function([3.3e9, 50.0])
    print(f"\nBaseline objective: f(x) = {baseline_error:.2f}")

    # Set up optimization
    print("\n" + "=" * 80)
    print(" Running Optimization (2 parameters)")
    print("=" * 80)

    bounds = [
        (1.65e9, 4.95e9),  # Ne: 3.3e9 ± 50%
        (10.0, 200.0)       # E field
    ]

    print("\nThis should complete quickly (~2-5 minutes)...")
    print("Progress:\n")

    start_time = time.time()

    result = differential_evolution(
        objective_function,
        bounds,
        strategy='best1bin',
        maxiter=50,
        popsize=10,
        tol=0.01,
        atol=0.0,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=42,
        disp=True,
        workers=1,
        updating='deferred',
        polish=True  # Enable polishing for 2-param case
    )

    elapsed = time.time() - start_time

    print("\n" + "=" * 80)
    print(f" OPTIMIZATION COMPLETE ({elapsed/60:.1f} minutes)")
    print("=" * 80)

    if result.success:
        print("\n✓ Optimization converged successfully")
    else:
        print("\n⚠ Optimization stopped but found best solution")

    print(f"\nFinal objective value: {result.fun:.6f}")
    print(f"Improvement: {(baseline_error - result.fun) / baseline_error * 100:.1f}%")
    print(f"Function evaluations: {result.nfev}")
    print(f"Iterations: {result.nit}")

    # Extract optimized parameters
    opt_ne = result.x[0]
    opt_E = result.x[1]

    print(f"\nOptimized Parameters:")
    print(f"  Ne: {opt_ne:.3e} cm^-3 (baseline: 3.3e9, change: {(opt_ne/3.3e9 - 1)*100:+.1f}%)")
    print(f"  E:  {opt_E:.1f} V/cm (baseline: 50, change: {(opt_E/50 - 1)*100:+.1f}%)")

    # Run final simulation with optimized parameters
    print("\n" + "=" * 80)
    print(" FINAL RESULTS WITH OPTIMIZED Ne AND E")
    print("=" * 80)

    final_results = run_simulation(opt_ne, opt_E, verbose=True)

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

        # Compare to baseline
        print("\n" + "=" * 80)
        print(" COMPARISON: Baseline vs Optimized")
        print("=" * 80)
        print("\n                Baseline      Optimized     Change")
        print("-" * 60)
        for species in ['H', 'CH', 'C2']:
            base = baseline[species]
            opt = final_results[species]
            change = (opt - base) / base * 100
            arrow = "↑" if change > 0 else "↓"
            print(f"{species:5s}:      {base:.2e}    {opt:.2e}    {change:+6.1f}% {arrow}")

    # Save results
    with open('optimization_ne_E_results.txt', 'w') as f:
        f.write("2-PARAMETER OPTIMIZATION RESULTS (Ne and E only)\n")
        f.write("=" * 80 + "\n\n")
        f.write("All reaction rates kept at literature values.\n\n")
        f.write(f"Optimized Ne: {opt_ne:.6e} cm^-3\n")
        f.write(f"Optimized E:  {opt_E:.6f} V/cm\n\n")
        f.write(f"Objective value: {result.fun:.6f}\n")
        f.write(f"Improvement over baseline: {(baseline_error - result.fun) / baseline_error * 100:.1f}%\n\n")
        f.write("=" * 80 + "\n")
        f.write("FINAL RESULTS\n")
        f.write("=" * 80 + "\n\n")
        if final_results:
            for species in ['H', 'CH', 'C2']:
                current = final_results[species]
                target = TARGETS[species]
                ratio = current / target
                f.write(f"{species}: {current:.6e} cm^-3  (target: {target:.2e}, ratio: {ratio:.2f}x)\n")

    print("\n✓ Results saved to optimization_ne_E_results.txt")
    print("\n" + "=" * 80)
    print(" NEXT STEPS")
    print("=" * 80)
    print("\nIf targets are close:")
    print("  → Success! Physical parameters alone can match experiment")
    print("\nIf targets still far:")
    print("  → Need to also tune reaction rates within literature bounds")
    print("  → Run full 31-parameter optimization with optimized Ne and E as starting point")


if __name__ == '__main__':
    main()
