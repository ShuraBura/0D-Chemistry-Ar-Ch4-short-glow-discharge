#!/usr/bin/env python3
"""
Run the BEST optimized model from optimization_results_targeted
"""

import numpy as np
from scipy.integrate import solve_ivp
import time
import json
import os

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


def load_best_model():
    """Load the best model parameters from optimization results."""

    # Find the best model file (lowest objective value)
    results_dir = 'optimization_results_targeted'
    best_file = None
    best_objective = float('inf')

    for filename in os.listdir(results_dir):
        if filename.startswith('best_iteration') and filename.endswith('.json'):
            # Extract objective value from filename
            try:
                obj_str = filename.split('_f')[1].split('.json')[0]
                obj_val = float(obj_str)
                if obj_val < best_objective:
                    best_objective = obj_val
                    best_file = filename
            except:
                continue

    if best_file is None:
        raise FileNotFoundError("No best model found in optimization_results_targeted/")

    print(f"\nLoading best model: {best_file}")
    print(f"  Objective value: {best_objective:.2f}")

    with open(os.path.join(results_dir, best_file), 'r') as f:
        model_data = json.load(f)

    return model_data


def run_simulation_with_params(model_data):
    """Run simulation using the optimized parameters."""

    # Extract parameters
    ne = model_data['Ne']
    E_field = model_data['E_field']
    rate_values = model_data['rate_values']

    print(f"\nModel parameters:")
    print(f"  Electron density (Ne): {ne:.3e} cm^-3")
    print(f"  Electric field (E): {E_field:.1f} V/cm")
    print(f"  Number of tuned rates: {len(rate_values)}")

    # Setup base parameters
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

    # Define base rates
    k = define_rates(params)

    # Apply optimized rate values
    db = get_complete_rate_database()
    for name, val in rate_values.items():
        if name in k:
            if name in db:
                val = np.clip(val, db[name].min, db[name].max)
            k[name] = val

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    # Setup initial conditions
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

    print("\n" + "=" * 80)
    print(" Running Optimized Simulation: 0 → 100 seconds")
    print("=" * 80)

    # Create ODE function
    ode_func = PlasmaODE_Optimized(params)

    # Run simulation
    t_span = (0, 100)
    t_eval = np.logspace(-3, 2, 50)  # Log-spaced time points for better resolution

    print(f"\nSolver: BDF (stiff)")
    print(f"Tolerances: rtol=1e-5, atol=1e-6")
    print("\nStarting...")

    start_time = time.time()

    sol = solve_ivp(
        ode_func,
        t_span,
        y0,
        method='BDF',
        t_eval=t_eval,
        rtol=1e-5,
        atol=1e-6,
        max_step=10.0
    )

    elapsed = time.time() - start_time

    if not sol.success:
        print(f"\n✗ Simulation FAILED: {sol.message}")
        return None

    print("\n" + "=" * 80)
    print(f" ✓ COMPLETED in {elapsed:.2f} seconds")
    print("=" * 80)
    print(f"  Status: {sol.message}")
    print(f"  Function evaluations: {sol.nfev}")
    print(f"  Jacobian evaluations: {sol.njev}")

    return sol, params


def display_results(sol, params, model_data):
    """Display final densities and compare to targets."""

    species = params['species']
    y_final = sol.y[:, -1]

    def get_density(name):
        try:
            idx = species.index(name)
            return y_final[idx]
        except ValueError:
            return 0.0

    # Experimental targets
    targets = {
        'H': 5.18e13,
        'CH': 1.0e9,
        'C2': 1.3e11,
    }

    print("\n" + "=" * 80)
    print(" Target Species Comparison")
    print("=" * 80)
    print(f"\n{'Species':<10} {'Simulation':<15} {'Target':<15} {'Ratio':<10}")
    print("-" * 60)

    for sp in ['H', 'CH', 'C2']:
        sim_val = get_density(sp)
        target_val = targets[sp]
        ratio = sim_val / target_val

        # Color coding based on how close to target
        if 0.5 <= ratio <= 2.0:
            status = "✓ GOOD"
        elif 0.2 <= ratio <= 5.0:
            status = "~ OK"
        else:
            status = "✗ OFF"

        print(f"{sp:<10} {sim_val:<15.3e} {target_val:<15.3e} {ratio:<10.2f}x {status}")

    # Additional important species
    print("\n" + "=" * 80)
    print(" Other Key Species")
    print("=" * 80)

    key_species = {
        'Ions': ['e', 'ArPlus', 'CH3Plus', 'CH5Plus', 'ArHPlus'],
        'Radicals': ['CH3', 'CH2', 'C'],
        'Molecules': ['H2', 'CH4', 'C2H2', 'C2H4', 'C2H6'],
        'Excited': ['ArStar']
    }

    for category, species_list in key_species.items():
        print(f"\n{category}:")
        for sp in species_list:
            dens = get_density(sp)
            if dens > 0:
                print(f"  {sp:<10s}: {dens:12.4e} cm^-3")

    print("\n" + "=" * 80)
    print(" Optimization Quality")
    print("=" * 80)
    print(f"  Objective value: {model_data.get('objective', 'N/A')}")

    # Calculate current objective for verification
    weights = {'H': 1.0, 'CH': 20.0, 'C2': 3.0}
    error = 0.0
    for sp in ['H', 'CH', 'C2']:
        rel_error = (get_density(sp) - targets[sp]) / targets[sp]
        error += weights[sp] * rel_error ** 2

    print(f"  Computed error: {error:.2f}")
    print(f"  Weighted by: H=1.0, CH=20.0, C2=3.0")


def main():
    print("=" * 80)
    print(" Running BEST Optimized Ar/CH4 Plasma Model")
    print("=" * 80)

    # Load best model
    model_data = load_best_model()

    # Run simulation
    result = run_simulation_with_params(model_data)

    if result is None:
        print("\nSimulation failed!")
        return

    sol, params = result

    # Display results
    display_results(sol, params, model_data)

    print("\n" + "=" * 80)
    print(" Simulation Complete")
    print("=" * 80)


if __name__ == '__main__':
    main()
