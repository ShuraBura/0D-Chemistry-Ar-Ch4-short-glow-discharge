#!/usr/bin/env python3
"""
Quick test of plasma simulation - shorter time span
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from define_rates import define_rates
from build_reactions import build_reactions
from odefun import PlasmaODE


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
    set_density('H', 1e11)

    return y0


# Quick test parameters
params = {
    'E_field': 50,
    'L_discharge': 0.45,
    'ne': 1e10,
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

print("Quick Test - Building reactions...")
params['k'] = define_rates(params)
params['R'], params['tags'] = build_reactions(params)
print(f"  {len(params['R'])} reactions for {len(params['species'])} species")

y0 = calculate_initial_densities(params)
ode_func = PlasmaODE(params)

print("\nRunning short simulation (0 to 0.1 s)...")
start = time.time()

sol = solve_ivp(
    ode_func,
    (0, 0.1),  # Only 0.1 seconds!
    y0,
    method='BDF',
    t_eval=[0, 0.01, 0.1],
    rtol=1e-6,
    atol=1e-7
)

elapsed = time.time() - start

print(f"\n✅ Success! Solved in {elapsed:.2f} seconds")
print(f"   Function evaluations: {sol.nfev}")
print(f"   Final H density: {sol.y[params['species'].index('H'), -1]:.3e} cm^-3")
