#!/usr/bin/env python3
"""
3-parameter optimization: Te, Ne, and E field.
Keep all reaction rates at literature values.

This tests if targets are achievable by tuning physical parameters
with our new temperature-dependent rates.

Te range: 0.3 to 7 eV
Ne range: 1e8 to 5e9 cm⁻³
E field range: 10 to 500 V/cm
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import time
import json

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

# Track best result
best_result = {'objective': 1e10, 'params': None, 'densities': None}


def run_simulation(Te, ne, E_field, verbose=False):
    """
    Run simulation with given Te, Ne and E field.
    All rates use temperature-dependent values from define_rates.py.

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
            'Te': Te,  # ← NEW! Electron temperature for rate calculation
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

        # Load rates with temperature dependence
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
    Objective function for 3-parameter optimization.

    x[0] = Te (electron temperature)
    x[1] = Ne (electron density)
    x[2] = E field
    """
    Te = x[0]
    ne = x[1]
    E_field = x[2]

    # Run simulation
    results = run_simulation(Te, ne, E_field, verbose=verbose)

    if results is None:
        return 1e10  # Penalize failed simulations

    # Calculate weighted error
    # CH weighted 20x (hardest to match, 59x off)
    # C2 weighted 3x (7x off)
    # H weighted 1x (0.66x, closest)
    weights = {
        'H': 1.0,
        'CH': 20.0,
        'C2': 3.0
    }

    error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        error += weights[species] * rel_error ** 2

    # Track best result
    global best_result
    if error < best_result['objective']:
        best_result['objective'] = error
        best_result['params'] = {'Te': Te, 'Ne': ne, 'E_field': E_field}
        best_result['densities'] = results

        print(f"\n  *** NEW BEST: f(x) = {error:.2f}")
        print(f"      Te: {Te:.2f} eV, Ne: {ne:.2e} cm⁻³, E: {E_field:.1f} V/cm")
        print(f"      H: {results['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {results['H']/TARGETS['H']:.2f}x)")
        print(f"      CH: {results['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {results['CH']/TARGETS['CH']:.2f}x)")
        print(f"      C2: {results['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {results['C2']/TARGETS['C2']:.2f}x)\n")

    return error


def main():
    print("=" * 80)
    print(" 3-PARAMETER OPTIMIZATION: Te, Ne, and E field")
    print("=" * 80)
    print("\nStrategy: Tune physical parameters with TEMPERATURE-DEPENDENT rates")
    print("          Testing our newly implemented Te-dependent define_rates.py!")
    print("\nParameters:")
    print("  Te (electron temperature): [0.3, 7.0] eV")
    print("  Ne (electron density): [1e8, 5e9] cm^-3")
    print("  E field: [10, 500] V/cm")
    print("\nAll reaction rates temperature-dependent via define_rates.py")
    print("Electron-impact, ionization, and recombination rates scale with Te")

    print("\n Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm^-3")
    print(f"  CH: {TARGETS['CH']:.2e} cm^-3  (currently 59x too high)")
    print(f"  C2: {TARGETS['C2']:.2e} cm^-3  (currently 7x too high)")

    # Run baseline with Te=1.0, Ne=3.3e9, E=50
    print("\n" + "=" * 80)
    print(" BASELINE (Te=1.0 eV, Ne=3.3e9, E=50)")
    print("=" * 80)

    baseline = run_simulation(1.0, 3.3e9, 50.0, verbose=True)
    if baseline:
        print("\nBaseline Densities:")
        for species in ['H', 'CH', 'C2']:
            current = baseline[species]
            target = TARGETS[species]
            ratio = current / target
            print(f"  {species:5s}: {current:.2e} cm^-3  (target: {target:.2e}, ratio: {ratio:.2f}x)")

    baseline_error = objective_function([1.0, 3.3e9, 50.0])
    print(f"\nBaseline objective: f(x) = {baseline_error:.2f}")

    # Set up optimization
    print("\n" + "=" * 80)
    print(" Running Optimization (3 parameters)")
    print("=" * 80)

    bounds = [
        (0.3, 7.0),      # Te: 0.3-7 eV
        (1e8, 5e9),      # Ne: 1e8 to 5e9 cm⁻³
        (10.0, 500.0)    # E field: 10-500 V/cm
    ]

    print("\nThis will test if Te variation can help match targets...")
    print("Higher Te → stronger electron-impact/ionization, weaker recombination")
    print("This may change the CH/C2 balance!\n")
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
        polish=True
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
    opt_Te = result.x[0]
    opt_ne = result.x[1]
    opt_E = result.x[2]

    print(f"\nOptimized Parameters:")
    print(f"  Te: {opt_Te:.2f} eV (baseline: 1.0, change: {(opt_Te/1.0 - 1)*100:+.1f}%)")
    print(f"  Ne: {opt_ne:.3e} cm^-3 (baseline: 3.3e9, change: {(opt_ne/3.3e9 - 1)*100:+.1f}%)")
    print(f"  E:  {opt_E:.1f} V/cm (baseline: 50, change: {(opt_E/50 - 1)*100:+.1f}%)")

    # Calculate E/N for validation
    n_total = 9.66e15  # cm⁻³
    E_over_N_Td = (opt_E / n_total) * 1e17
    print(f"\n  E/N: {E_over_N_Td:.1f} Td")
    print(f"  Expected Te from E/N ~ 40-100 Td: 2-5 eV")
    if 1.0 <= opt_Te <= 6.0:
        print(f"  → Te={opt_Te:.2f} eV is physically reasonable ✓")
    else:
        print(f"  → Te={opt_Te:.2f} eV may be outside typical range")

    # Run final simulation with optimized parameters
    print("\n" + "=" * 80)
    print(" FINAL RESULTS WITH OPTIMIZED Te, Ne, AND E")
    print("=" * 80)

    final_results = run_simulation(opt_Te, opt_ne, opt_E, verbose=True)

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
    result_data = {
        'optimized_parameters': {
            'Te': float(opt_Te),
            'Ne': float(opt_ne),
            'E_field': float(opt_E),
            'E_over_N_Td': float(E_over_N_Td)
        },
        'objective_value': float(result.fun),
        'baseline_objective': float(baseline_error),
        'improvement_percent': float((baseline_error - result.fun) / baseline_error * 100),
        'function_evaluations': int(result.nfev),
        'iterations': int(result.nit),
        'success': result.success
    }

    if final_results:
        result_data['final_densities'] = {
            'H': float(final_results['H']),
            'CH': float(final_results['CH']),
            'C2': float(final_results['C2'])
        }
        result_data['ratios_vs_target'] = {
            'H': float(final_results['H'] / TARGETS['H']),
            'CH': float(final_results['CH'] / TARGETS['CH']),
            'C2': float(final_results['C2'] / TARGETS['C2'])
        }

    with open('optimization_Te_Ne_E_results.json', 'w') as f:
        json.dump(result_data, f, indent=2)

    print("\n✓ Results saved to optimization_Te_Ne_E_results.json")

    print("\n" + "=" * 80)
    print(" TEMPERATURE DEPENDENCE IMPACT")
    print("=" * 80)
    print("\nThis optimization used our newly implemented temperature-dependent rates!")
    print("Higher Te increases: electron-impact (99×), ionization (770×)")
    print("Higher Te decreases: recombination (0.62×)")
    print("\nIf Te varied significantly from 1.0 eV, temperature dependence was crucial.")

    print("\n" + "=" * 80)
    print(" NEXT STEPS")
    print("=" * 80)
    print("\nIf targets are close:")
    print("  → Success! Physical parameters + Te dependence can match experiment")
    print("\nIf targets still far:")
    print("  → Need to also tune reaction rates within literature bounds")
    print("  → Run full multi-parameter optimization with rates + Te/Ne/E")


if __name__ == '__main__':
    main()
