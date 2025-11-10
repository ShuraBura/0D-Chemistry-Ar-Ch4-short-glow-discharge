#!/usr/bin/env python3
"""
Direct test: Is the drift term actually being applied?
"""

import numpy as np
import json

print("=" * 80)
print("DRIFT TERM APPLICATION TEST")
print("=" * 80)
print()

# Load best result to get params
with open('optimization_results_charge_balanced/best_f27.0.json', 'r') as f:
    result = json.load(f)

# Import ODE
from odefun_optimized import PlasmaODE_Optimized
from build_reactions import build_reactions

# Build params exactly as optimization does
params = {
    'E_field': result['E_field'],
    'L_discharge': 0.45,
    'ne': result['Ne'],
    'Te': result['Te'],
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
    'k': result['rate_values']
}

params['R'], params['tags'] = build_reactions(params)

# Create ODE
ode = PlasmaODE_Optimized(params)

print(f"ODE created successfully")
print(f"  H_drift_gain: {ode.H_drift_gain:.2e} cm⁻³/s")
print(f"  H_idx: {ode.H_idx}")
print()

# Get state vector from result
species = params['species']
y = np.array([result['all_densities'][sp] for sp in species])

H_idx = species.index('H')
print(f"Current state:")
print(f"  H density: {y[H_idx]:.2e} cm⁻³")
print()

# Call ODE to get dydt
dydt = ode(0, y)

print(f"ODE evaluation:")
print(f"  dH/dt: {dydt[H_idx]:.2e} cm⁻³/s")
print()

# Now let's manually calculate what dH/dt SHOULD be
# by removing the drift term temporarily and comparing

# Save drift value
drift_value = ode.H_drift_gain

# Temporarily zero out drift
ode.H_drift_gain = 0.0
dydt_no_drift = ode(0, y)
dH_dt_no_drift = dydt_no_drift[H_idx]

# Restore drift
ode.H_drift_gain = drift_value
dydt_with_drift = ode(0, y)
dH_dt_with_drift = dydt_with_drift[H_idx]

print("=" * 80)
print("COMPARISON TEST")
print("=" * 80)
print()

print(f"dH/dt WITHOUT drift: {dH_dt_no_drift:.2e} cm⁻³/s")
print(f"dH/dt WITH drift:    {dH_dt_with_drift:.2e} cm⁻³/s")
print(f"Difference:          {dH_dt_with_drift - dH_dt_no_drift:.2e} cm⁻³/s")
print()

print(f"Expected difference: {drift_value:.2e} cm⁻³/s")
print(f"Actual difference:   {dH_dt_with_drift - dH_dt_no_drift:.2e} cm⁻³/s")
print()

if abs((dH_dt_with_drift - dH_dt_no_drift) - drift_value) < 1e10:
    print("✓ DRIFT TERM IS BEING APPLIED CORRECTLY!")
    print()
    print("The drift term IS working, but H is still low because:")
    print("  1. System is at steady state (dH/dt ≈ 0)")
    print("  2. Chemical production + drift = Total loss")
    print("  3. Need to investigate what sinks are balancing the drift")
else:
    print("✗ DRIFT TERM IS NOT BEING APPLIED!")
    print()
    print("The difference should equal drift_value but doesn't.")
    print("There's a bug in the ODE implementation.")

print()
print("=" * 80)
print("STEADY STATE CHECK")
print("=" * 80)
print()

print(f"dH/dt at current state: {dydt[H_idx]:.2e} cm⁻³/s")
print(f"Relative rate: {dydt[H_idx]/y[H_idx]:.2e} s⁻¹")
print()

if abs(dydt[H_idx]) < 1e12:
    print("✓ System IS at steady state (dH/dt ≈ 0)")
    print()
    print("This means the ODE integration has converged.")
    print("The 'low' H is the actual steady-state value.")
    print()
    print("To understand WHY H is low, we need to identify:")
    print("  1. What H sinks are consuming the drift source?")
    print("  2. Are there implicit sinks not in the balance analysis?")
else:
    print("✗ System NOT at steady state")
    print("Integration hasn't converged yet.")
