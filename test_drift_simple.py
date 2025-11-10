#!/usr/bin/env python3
"""
Simpler test: Load ODE from a completed integration and check drift
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

print("=" * 80)
print("SIMPLE DRIFT TERM TEST")
print("=" * 80)
print()

# Load result
with open('optimization_results_charge_balanced/best_f27.0.json', 'r') as f:
    result = json.load(f)

# Import and setup
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from build_reactions import build_reactions

# Setup params as optimization does
params = {
    'E_field': result['E_field'],
    'L_discharge': 0.45,
    'ne': result['Ne'],
    'Te': result['Te'],
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star']
}

# Define rates
k = define_rates(params)

# Override with result rates
for name, val in result['rate_values'].items():
    if name in k:
        k[name] = val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Create ODE
ode = PlasmaODE_Optimized(params)
species = params['species']
H_idx = species.index('H')

print(f"✓ ODE created")
print(f"  H_drift_gain = {ode.H_drift_gain:.2e} cm⁻³/s")
print(f"  H_idx = {H_idx}")
print()

# Get state from result
y = np.array([result['all_densities'][sp] for sp in species])

print(f"State from result:")
print(f"  H = {y[H_idx]:.2e} cm⁻³")
print()

# Evaluate ODE
dydt = ode(0, y)

print(f"dH/dt at this state: {dydt[H_idx]:.2e} cm⁻³/s")
print()

# Test with and without drift
print("=" * 80)
print("DRIFT TERM PRESENCE TEST")
print("=" * 80)
print()

# Save and zero drift
original_drift = ode.H_drift_gain
ode.H_drift_gain = 0.0
dydt_no_drift = ode(0, y)

# Restore drift
ode.H_drift_gain = original_drift
dydt_with_drift = ode(0, y)

difference = dydt_with_drift[H_idx] - dydt_no_drift[H_idx]

print(f"dH/dt WITHOUT drift: {dydt_no_drift[H_idx]:.2e} cm⁻³/s")
print(f"dH/dt WITH drift:    {dydt_with_drift[H_idx]:.2e} cm⁻³/s")
print(f"Difference:          {difference:.2e} cm⁻³/s")
print(f"Expected (drift):    {original_drift:.2e} cm⁻³/s")
print()

match = abs(difference - original_drift) / original_drift < 0.01

if match:
    print("✓✓ DRIFT TERM IS APPLIED CORRECTLY!")
    print()
    print(f"Since dH/dt ≈ 0 at steady state, the drift ({original_drift:.2e})")
    print(f"must be balanced by H sinks.")
else:
    print("✗✗ DRIFT TERM NOT APPLIED!")
    print(f"   Expected difference: {original_drift:.2e}")
    print(f"   Actual difference:   {difference:.2e}")
    print(f"   Ratio: {difference/original_drift:.3f}")

print()
print("=" * 80)
print("STEADY STATE CHECK")
print("=" * 80)
print()

at_steady_state = abs(dydt[H_idx]) < 1e12

print(f"dH/dt = {dydt[H_idx]:.2e} cm⁻³/s")

if at_steady_state:
    print("✓ At steady state (dH/dt ≈ 0)")
    print()
    print("If drift IS applied but H is low, then:")
    print("  → There must be LARGE H sinks balancing the drift")
    print("  → Need to find these sinks in the chemistry")
else:
    print("✗ NOT at steady state")
    print(f"  Relative rate: {dydt[H_idx]/y[H_idx]:.2e} s⁻¹")
