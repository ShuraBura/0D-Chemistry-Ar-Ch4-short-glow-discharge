#!/usr/bin/env python3
"""
optimize_to_targets.py - Tune reaction rates to match experimental densities
Constraints: All rates must stay within literature ranges
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.integrate import solve_ivp
import time

from rate_database import get_rate_database, get_tunable_parameters
from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


# Experimental target densities (from MATLAB comments, line 10)
TARGETS = {
    'H': 5.18e13,   # cm^-3, measured at 4.5mm from cathode
    'CH': 1e9,       # cm^-3
    'C2': 1.3e11,    # cm^-3
}


def define_rates_from_dict(params, rate_values):
    """Create rate dict from parameter values."""
    # Start with ALL rates from original define_rates
    k = define_rates(params)

    # Override with tuned values
    for key, val in rate_values.items():
        k[key] = val

    return k


def define_rates_from_dict_OLD(params, rate_values):
    """OLD VERSION - Create rate dict from parameter values."""
    k = {}
    rate_db = get_rate_database()

    # Start with database values
    for key, rate_obj in rate_db.items():
        k[key] = rate_obj.value

    # Override with tuned values
    for key, val in rate_values.items():
        k[key] = val

    # Add fixed rates that don't have ranges (mobilities, etc.)


def calculate_initial_densities(params):
    """Initial conditions."""
    species = params['species']
    y0 = np.ones(len(species)) * 1e3

    def set_density(name, value):
        try:
            y0[species.index(name)] = value
        except ValueError:
            pass

    set_density('e', params['ne'])
    set_density('Ar', 0.85 * 9.66e15)
    set_density('CH4', 0.15 * 9.66e15)
    set_density('ArPlus', 1e7)
    set_density('CH4Plus', 1e5)
    set_density('CH3Plus', 1e5)
    set_density('H2', 1e12)
    set_density('ArStar', 5e6)
    set_density('H', 1e11)
    set_density('C2', 5e7)
    set_density('CH', 5e4)

    return y0


def run_simulation(rate_values, params_base, verbose=False):
    """
    Run simulation with given rate values.

    Returns final densities for H, CH, C2
    """
    params = params_base.copy()
    params['k'] = define_rates_from_dict(params, rate_values)
    params['R'], params['tags'] = build_reactions(params)

    y0 = calculate_initial_densities(params)
    ode_func = PlasmaODE_Optimized(params)

    try:
        sol = solve_ivp(
            ode_func,
            (0, 100),
            y0,
            method='BDF',
            t_eval=[100],
            rtol=1e-6,
            atol=1e-7
        )

        if not sol.success:
            if verbose:
                print(f"  Solver failed: {sol.message}")
            return None

        # Extract final densities
        results = {}
        for sp in ['H', 'CH', 'C2']:
            idx = params['species'].index(sp)
            results[sp] = sol.y[idx, -1]

        return results

    except Exception as e:
        if verbose:
            print(f"  Simulation error: {e}")
        return None


def objective_function(x, param_names, params_base, verbose=False):
    """
    Objective function to minimize: weighted error vs targets.

    Parameters:
    -----------
    x : array
        Parameter values to optimize
    param_names : list
        Names of parameters being optimized
    params_base : dict
        Base parameters
    verbose : bool
        Print diagnostics

    Returns:
    --------
    error : float
        Weighted sum of squared relative errors
    """
    # Create rate dict from optimization variables
    rate_values = {name: val for name, val in zip(param_names, x)}

    # Run simulation
    results = run_simulation(rate_values, params_base, verbose)

    if results is None:
        return 1e10  # Penalty for failed simulation

    # Calculate weighted error
    error = 0
    weights = {'H': 1.0, 'CH': 10.0, 'C2': 2.0}  # CH is hardest to match

    for sp in ['H', 'CH', 'C2']:
        target = TARGETS[sp]
        actual = results[sp]
        rel_error = (actual - target) / target
        error += weights[sp] * rel_error**2

    if verbose:
        print(f"  Error: {error:.3f}")
        for sp in ['H', 'CH', 'C2']:
            print(f"    {sp}: {results[sp]:.3e} (target: {TARGETS[sp]:.3e})")

    return error


def optimize_parameters():
    """Main optimization routine."""
    print("=" * 70)
    print(" Parameter Optimization to Match Experimental Targets")
    print("=" * 70)

    # Setup base parameters
    params_base = {
        'E_field': 50,
        'L_discharge': 0.45,
        'ne': 1e10,
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

    # Get tunable parameters
    tunable = get_tunable_parameters()
    rate_db = get_rate_database()

    print("\nTunable Parameters (within literature ranges):")
    param_names = []
    x0 = []
    bounds = []

    for key in list(tunable.keys())[:5]:  # Start with top 5 most important
        if key == 'loss_H_drift_gain':
            continue  # Handle separately
        param_names.append(key)
        rate_obj = rate_db[key]
        x0.append(rate_obj.value)
        bounds.append((rate_obj.min, rate_obj.max))
        print(f"  {key}:")
        print(f"    Range: [{rate_obj.min:.2e}, {rate_obj.max:.2e}]")
        print(f"    Current: {rate_obj.value:.2e}")
        print(f"    Note: {tunable[key]}")

    print(f"\nOptimizing {len(param_names)} parameters...")
    print(f"\nTargets:")
    for sp, target in TARGETS.items():
        print(f"  {sp}: {target:.3e} cm^-3")

    # Run baseline
    print("\nBaseline (before optimization):")
    baseline_rates = {name: val for name, val in zip(param_names, x0)}
    baseline_results = run_simulation(baseline_rates, params_base, verbose=True)

    # Run optimization
    print("\nRunning optimization...")
    print("  Method: Differential Evolution (global optimizer)")
    print("  Constraints: All rates within literature bounds")

    start = time.time()

    result = differential_evolution(
        objective_function,
        bounds,
        args=(param_names, params_base, False),
        strategy='best1bin',
        maxiter=50,
        popsize=10,
        tol=0.01,
        disp=True,
        workers=1
    )

    elapsed = time.time() - start

    print(f"\nOptimization completed in {elapsed:.1f} seconds")
    print(f"  Success: {result.success}")
    print(f"  Final error: {result.fun:.6f}")

    # Display optimized parameters
    print("\n" + "=" * 70)
    print(" Optimized Parameters")
    print("=" * 70)

    optimized_rates = {name: val for name, val in zip(param_names, result.x)}

    for i, name in enumerate(param_names):
        rate_obj = rate_db[name]
        old_val = x0[i]
        new_val = result.x[i]
        change = (new_val / old_val - 1) * 100

        print(f"\n{name}:")
        print(f"  Old: {old_val:.3e}")
        print(f"  New: {new_val:.3e} ({change:+.1f}%)")
        print(f"  Range: [{rate_obj.min:.2e}, {rate_obj.max:.2e}]")
        print(f"  Within bounds: {rate_obj.is_within_range(new_val)} ✓")

    # Final simulation with optimized parameters
    print("\n" + "=" * 70)
    print(" Final Results (optimized)")
    print("=" * 70)

    final_results = run_simulation(optimized_rates, params_base, verbose=False)

    print("\nDensity Comparison:")
    print(f"{'Species':<10} {'Target':<15} {'Baseline':<15} {'Optimized':<15} {'Error'}")
    print("-" * 70)

    for sp in ['H', 'CH', 'C2']:
        target = TARGETS[sp]
        baseline = baseline_results[sp]
        optimized = final_results[sp]
        error = abs(optimized - target) / target * 100

        print(f"{sp:<10} {target:.3e}  {baseline:.3e}  {optimized:.3e}  {error:.1f}%")

    print("=" * 70)

    return result, param_names


if __name__ == '__main__':
    result, param_names = optimize_parameters()
