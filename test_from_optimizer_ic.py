#!/usr/bin/env python3
"""
Test if we reach the same steady state starting from optimizer's initial condition (H=1e11)
vs starting from the saved transient state (H=1e13)
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from build_reactions import build_reactions

print("=" * 80)
print("STEADY STATE TEST: INITIAL CONDITION INDEPENDENCE")
print("=" * 80)
print()

# Load best result
with open('optimization_results_charge_balanced/best_f41.3.json', 'r') as f:
    result = json.load(f)

# Setup ODE with optimized rates
species_list = ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star']

mobilities = {
    'ArPlus': 1.54e3, 'CH4Plus': 1.54e3, 'CH3Plus': 1.54e3,
    'CH5Plus': 1.54e3, 'ArHPlus': 1.54e3, 'H3Plus': 1.54e3,
    'CH2Plus': 1.54e3, 'C2H5Plus': 1.54e3, 'C2H4Plus': 1.54e3,
    'C2H3Plus': 1.54e3, 'C2HPlus': 1.54e3, 'H2Plus': 1.54e3,
    'CHPlus': 1.54e3, 'CH3Minus': 1.54e3, 'HMinus': 1.54e3
}

params = {
    'E_field': result['E_field'],
    'L_discharge': 0.45,
    'ne': result['Ne'],
    'Te': result['Te'],
    'species': species_list,
    'T': 400.0,
    'Tgas': 400.0,
    'pressure': 500.0,
    'mobilities': mobilities
}

k = define_rates(params)
for name, val in result['rate_values'].items():
    if name in k:
        k[name] = val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)
ode = PlasmaODE_Optimized(params)

H_idx = species_list.index('H')
print(f"✓ ODE created with optimized rates from best_f41.3.json")
print(f"  H_drift_gain = {ode.H_drift_gain:.2e} cm⁻³/s")
print()

# Test 1: Start from saved state (H = 1.03e13)
print("=" * 80)
print("TEST 1: Start from saved state")
print("=" * 80)
y1 = np.array([result['all_densities'][sp] for sp in species_list])
print(f"Initial H = {y1[H_idx]:.2e} cm⁻³")

sol1 = solve_ivp(
    ode, (0, 100), y1,
    method='BDF',
    t_eval=np.logspace(-6, 2, 100),
    rtol=1e-6, atol=1e-8, max_step=1.0
)

if sol1.success:
    print(f"Final H = {sol1.y[H_idx, -1]:.2e} cm⁻³")
    print(f"Growth: {sol1.y[H_idx, -1] / y1[H_idx]:.1f}×")
    H_ss1 = sol1.y[H_idx, -1]
else:
    print(f"Integration failed: {sol1.message}")
    H_ss1 = None

print()

# Test 2: Start from optimizer initial condition (H = 1e11)
print("=" * 80)
print("TEST 2: Start from optimizer IC")
print("=" * 80)

# Use same initial condition as optimizer
y2 = y1.copy()
y2[H_idx] = 1e11

print(f"Initial H = {y2[H_idx]:.2e} cm⁻³ (optimizer default)")

sol2 = solve_ivp(
    ode, (0, 100), y2,
    method='BDF',
    t_eval=np.logspace(-6, 2, 100),
    rtol=1e-6, atol=1e-8, max_step=1.0
)

if sol2.success:
    print(f"Final H = {sol2.y[H_idx, -1]:.2e} cm⁻³")
    print(f"Growth: {sol2.y[H_idx, -1] / y2[H_idx]:.0f}×")
    H_ss2 = sol2.y[H_idx, -1]
else:
    print(f"Integration failed: {sol2.message}")
    H_ss2 = None

print()

# Compare
print("=" * 80)
print("COMPARISON")
print("=" * 80)
print()

if H_ss1 and H_ss2:
    print(f"Steady state from saved IC:     {H_ss1:.2e} cm⁻³")
    print(f"Steady state from optimizer IC: {H_ss2:.2e} cm⁻³")
    print(f"Ratio: {H_ss2/H_ss1:.3f}")
    print()

    if abs(H_ss2/H_ss1 - 1.0) < 0.05:
        print("✓ Same steady state (within 5%)")
        print("  → Steady state is independent of initial condition")
        print()
        print("CONCLUSION: The optimizer SHOULD find H ~ 1.8e14 with tight tolerances!")
        print("If it's still showing H ~ 1e13, either:")
        print("  1. Optimizer is finding DIFFERENT rate constants")
        print("  2. Integration still not fully converged")
        print("  3. Some other issue")
    else:
        print("✗ Different steady states!")
        print("  → This is unexpected - steady state should be unique")
        print("  → Possible causes:")
        print("    - Integration didn't fully converge in one case")
        print("    - Nonlinear effects / multiple steady states")
else:
    print("✗ One or both integrations failed")
