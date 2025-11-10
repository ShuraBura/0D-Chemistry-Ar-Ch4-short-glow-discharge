#!/usr/bin/env python3
"""
Test that H_drift_gain is actually being applied in the ODE
"""

import numpy as np
from scipy.integrate import solve_ivp
import json

# Load a recent optimization result
with open('optimization_results_charge_balanced/best_f61.8.json', 'r') as f:
    result = json.load(f)

# Load the ODE function
from odefun_optimized import PlasmaODE_Optimized
from build_reactions import build_reactions

species, R, stoich_matrix = build_reactions()

# Build params
params = {
    'species': species,
    'R': R,
    'stoich': stoich_matrix,
    'k': result['rate_values'],
    'tags': {}
}

# Create ODE
ode = PlasmaODE_Optimized(params)

print("=" * 80)
print("H DRIFT TERM TEST")
print("=" * 80)
print()
print(f"H_drift_gain value: {ode.H_drift_gain:.2e} cm⁻³/s")
print()

# Get final densities from result
y_final = np.array([result['all_densities'][sp] for sp in species])
H_idx = species.index('H')

print(f"H density from result: {y_final[H_idx]:.2e} cm⁻³")
print()

# Evaluate dydt at the final state
dydt = ode(0, y_final)

print(f"dH/dt at final state: {dydt[H_idx]:.2e} cm⁻³/s")
print()

# Check if it's close to zero (steady state)
if abs(dydt[H_idx]) < 1e12:
    print("✓ System is at steady state (dH/dt ≈ 0)")
elif dydt[H_idx] > 0:
    print(f"✗ H is INCREASING at {dydt[H_idx]:.2e} cm⁻³/s (not at steady state!)")
else:
    print(f"✗ H is DECREASING at {dydt[H_idx]:.2e} cm⁻³/s (not at steady state!)")

print()

# Now let's manually calculate what dH/dt should be
print("=" * 80)
print("MANUAL H BALANCE CALCULATION")
print("=" * 80)
print()

# Wall sticking
k_wall = result['rate_values'].get('stick_H_9_1', 0)
wall_loss = k_wall * y_final[H_idx]
print(f"Wall sticking loss: {wall_loss:.2e} cm⁻³/s (k={k_wall:.1f} s⁻¹)")

# Drift source
drift_source = ode.H_drift_gain
print(f"Drift source: {drift_source:.2e} cm⁻³/s")

# Rough estimate of net rate (ignoring chemical reactions for now)
net_rate_simple = drift_source - wall_loss
print(f"Net rate (drift - wall only): {net_rate_simple:.2e} cm⁻³/s")
print()

if net_rate_simple > 1e15:
    print(f"→ H should be INCREASING rapidly!")
    print(f"→ With this net rate, H would increase by {net_rate_simple * 0.001:.2e} cm⁻³ in 1 ms")
    print(f"→ Current H = {y_final[H_idx]:.2e}")
    print(f"→ After 1 ms: H ≈ {y_final[H_idx] + net_rate_simple * 0.001:.2e}")
elif abs(net_rate_simple) < 1e12:
    print("→ H is approximately at steady state (drift ≈ wall loss)")
else:
    print(f"→ H should be DECREASING")

print()
print("=" * 80)
print("CONCLUSION")
print("=" * 80)
print()

if abs(dydt[H_idx]) > 1e15:
    print("The ODE shows H is NOT at steady state.")
    print("This suggests:")
    print("  1. Integration time (100s) may be insufficient")
    print("  2. Or there's numerical stiffness preventing convergence")
    print("  3. Or the solver terminates early for some reason")
elif abs(net_rate_simple) > 1e15:
    print("Drift >> wall loss, but ODE shows steady state.")
    print("This suggests:")
    print("  1. Other H sinks are very large")
    print("  2. Or drift term is not being applied correctly")
else:
    print("System appears to be at steady state as expected.")
