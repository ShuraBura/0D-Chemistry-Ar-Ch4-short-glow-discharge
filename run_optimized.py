#!/usr/bin/env python3
"""
Run optimized plasma simulation with performance comparison
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


def calculate_initial_densities(params):
    """Calculate initial species densities."""
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

    return y0


print("=" * 70)
print(" OPTIMIZED Ar/CH4 Plasma Simulation")
print("=" * 70)

# Setup parameters
params = {
    'E_field': 50,
    'L_discharge': 0.45,
    'ne': 3.3e9,  # Experimental value from sheath analysis
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

print("\nBuilding reaction network...")
params['k'] = define_rates(params)
params['R'], params['tags'] = build_reactions(params)
print(f"  Species: {len(params['species'])}")
print(f"  Reactions: {len(params['R'])}")

y0 = calculate_initial_densities(params)

# Create optimized ODE function
print("\nInitializing optimized ODE solver...")
ode_func = PlasmaODE_Optimized(params)
print("  ✓ Sparse matrices built")
print("  ✓ Vectorization ready")

# Run simulation
print("\n" + "=" * 70)
print(" Running Simulation: 0 → 100 seconds")
print("=" * 70)

t_span = (0, 100)
t_eval = [0, 0.1, 1, 10, 100]

print(f"\nSolver: BDF (stiff)")
print(f"Tolerances: rtol=1e-6, atol=1e-7")
print("\nStarting...")

start_time = time.time()

sol = solve_ivp(
    ode_func,
    t_span,
    y0,
    method='BDF',
    t_eval=t_eval,
    rtol=1e-6,
    atol=1e-7
)

elapsed = time.time() - start_time

print("\n" + "=" * 70)
print(f" ✓ COMPLETED in {elapsed:.2f} seconds")
print("=" * 70)
print(f"  Status: {sol.message}")
print(f"  Function evaluations: {sol.nfev}")
print(f"  Jacobian evaluations: {sol.njev}")
print(f"  LU decompositions: {sol.nlu}")

# Display key results
print("\n" + "=" * 70)
print(" Final Densities (t = 100 s)")
print("=" * 70)

key_species = {
    'Radicals': ['H', 'CH', 'C2', 'CH3'],
    'Ions': ['e', 'ArPlus', 'CH3Plus', 'CH5Plus'],
    'Molecules': ['H2', 'CH4', 'C2H2', 'C2H4', 'C2H6'],
    'Excited': ['ArStar']
}

for category, species_list in key_species.items():
    print(f"\n{category}:")
    for sp in species_list:
        try:
            idx = params['species'].index(sp)
            print(f"  {sp:10s}: {sol.y[idx, -1]:12.4e} cm^-3")
        except ValueError:
            pass

print("\n" + "=" * 70)
print(" Simulation Summary")
print("=" * 70)
print(f"  Simulation time: 0 → 100 s")
print(f"  Wall clock time: {elapsed:.2f} s")
print(f"  Speedup factor: {100/elapsed:.1f}x faster than real-time")
print("=" * 70)
