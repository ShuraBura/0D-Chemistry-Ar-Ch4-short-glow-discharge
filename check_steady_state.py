#!/usr/bin/env python3
"""
Check if the optimization result is at steady state by re-integrating
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

print("=" * 80)
print("STEADY STATE VERIFICATION")
print("=" * 80)
print()

# Load result
with open('optimization_results_charge_balanced/best_f27.0.json', 'r') as f:
    result = json.load(f)

print("Loaded result f=27.0")
print(f"  H: {result['all_densities']['H']:.2e} cm⁻³")
print(f"  Te: {result['Te']:.3f} eV")
print(f"  Ne: {result['Ne']:.2e} cm⁻³")
print()

# Import required modules
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from build_reactions import build_reactions

# Build params as optimizer does
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
    'CHPlus': 1.54e3,
    'CH3Minus': 1.54e3, 'HMinus': 1.54e3  # Negative ions
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

print("Building ODE...")
try:
    k = define_rates(params)

    # Override with result rates
    for name, val in result['rate_values'].items():
        if name in k:
            k[name] = val

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    ode = PlasmaODE_Optimized(params)
    print(f"✓ ODE created")
    print(f"  H_drift_gain = {ode.H_drift_gain:.2e} cm⁻³/s")
    print()

except Exception as e:
    print(f"✗ Failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# Build state vector
y = np.array([result['all_densities'][sp] for sp in species_list])
H_idx = species_list.index('H')

print("Testing ODE at final state...")
dydt = ode(0, y)

print(f"At final state:")
print(f"  H = {y[H_idx]:.2e} cm⁻³")
print(f"  dH/dt = {dydt[H_idx]:.2e} cm⁻³/s")
print()

if abs(dydt[H_idx]) < 1e12:
    print("✓ At steady state (dH/dt ≈ 0)")
else:
    print(f"✗ NOT at steady state (relative rate: {dydt[H_idx]/y[H_idx]:.2e} s⁻¹)")

print()
print("=" * 80)
print("DRIFT TERM TEST")
print("=" * 80)
print()

# Test with/without drift
drift_orig = ode.H_drift_gain

ode.H_drift_gain = 0.0
dydt_no_drift = ode(0, y)

ode.H_drift_gain = drift_orig
dydt_with_drift = ode(0, y)

diff = dydt_with_drift[H_idx] - dydt_no_drift[H_idx]

print(f"dH/dt WITHOUT drift: {dydt_no_drift[H_idx]:.2e} cm⁻³/s")
print(f"dH/dt WITH drift:    {dydt_with_drift[H_idx]:.2e} cm⁻³/s")
print(f"Difference:          {diff:.2e} cm⁻³/s")
print(f"Expected (drift):    {drift_orig:.2e} cm⁻³/s")
print(f"Match: {abs(diff - drift_orig) / drift_orig < 0.01}")
print()

if abs(diff - drift_orig) / drift_orig < 0.01:
    print("✓✓✓ DRIFT IS APPLIED!")
    print()
    print("Now let's re-integrate from scratch to see if H changes:")
    print()

    # Start from a different initial condition with higher H
    y0 = y.copy()
    y0[H_idx] = 2.52e14  # Target H

    print(f"Initial H: {y0[H_idx]:.2e}")
    print("Integrating for 10 seconds...")

    sol = solve_ivp(
        ode,
        (0, 10),
        y0,
        method='BDF',
        rtol=1e-5,
        atol=1e-6,
        max_step=1.0
    )

    if sol.success:
        H_final = sol.y[H_idx, -1]
        print(f"Final H: {H_final:.2e} cm⁻³")
        print(f"Change: {H_final/y0[H_idx]:.3f}×")

        if H_final < 0.5 * y0[H_idx]:
            print()
            print("H DECREASED! There are strong H sinks pulling it down.")
        elif H_final > 2 * y0[H_idx]:
            print()
            print("H INCREASED! Drift is overwhelming sinks.")
        else:
            print()
            print("H relatively stable.")
    else:
        print(f"Integration failed: {sol.message}")
else:
    print("✗✗✗ DRIFT NOT APPLIED!")
    print(f"Ratio: {diff/drift_orig:.6f}")
