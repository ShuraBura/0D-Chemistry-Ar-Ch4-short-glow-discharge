#!/usr/bin/env python3
"""
MINIMAL TEST: What is dH/dt at the final state?
Load result, evaluate ODE, print dH/dt. That's it.
"""

import numpy as np
import json

print("=" * 80)
print("MINIMAL DRIFT VERIFICATION TEST")
print("=" * 80)
print()

# Load the best result
with open('optimization_results_charge_balanced/best_f27.0.json', 'r') as f:
    result = json.load(f)

print("Step 1: Load ODE module")
from odefun_optimized import PlasmaODE_Optimized

# Check what H_drift_gain is hardcoded in the class
print(f"  Checking PlasmaODE_Optimized source...")

# Read the source to verify
with open('odefun_optimized.py', 'r') as f:
    for i, line in enumerate(f, 1):
        if 'H_drift_gain' in line and '=' in line:
            print(f"  Line {i}: {line.strip()}")

print()

# The result JSON should have all we need
print("Step 2: Extract state from result JSON")
print(f"  H density in result: {result['all_densities']['H']:.2e} cm⁻³")
print(f"  Ne: {result['Ne']:.2e} cm⁻³")
print(f"  Te: {result['Te']:.3f} eV")
print()

# Try to reconstruct just enough to test
print("Step 3: Build minimal params to create ODE")

# Get species list from result
species = list(result['all_densities'].keys())
print(f"  Found {len(species)} species in result")

# Build minimal params
from build_reactions import build_reactions
from define_rates import define_rates

params = {
    'E_field': result['E_field'],
    'L_discharge': 0.45,  # Standard value
    'ne': result['Ne'],
    'Te': result['Te'],
    'species': species,
    'pressure': 500.0,  # mTorr
}

print("  Calling define_rates...")
try:
    # Need to add more params for define_rates
    params['T'] = 400.0  # K
    params['Tgas'] = 400.0  # K

    # Mobilities - get from result or use defaults
    params['mobilities'] = {
        'ArPlus': 1.54e3,
        'CH4Plus': 1.54e3,
        'CH3Plus': 1.54e3,
        'CH5Plus': 1.54e3,
        'ArHPlus': 1.54e3,
        'H3Plus': 1.54e3,
        'CH2Plus': 1.54e3,
        'C2H5Plus': 1.54e3,
        'C2H4Plus': 1.54e3,
        'C2H3Plus': 1.54e3,
        'C2HPlus': 1.54e3,
        'H2Plus': 1.54e3,
        'CHPlus': 1.54e3,
    }

    k = define_rates(params)

    # Override with optimized rates from result
    for name, val in result['rate_values'].items():
        if name in k:
            k[name] = val

    params['k'] = k

    print("  Calling build_reactions...")
    params['R'], params['tags'] = build_reactions(params)

    print("  Creating ODE...")
    ode = PlasmaODE_Optimized(params)

    print(f"✓ ODE created successfully")
    print()

except Exception as e:
    print(f"✗ Failed to create ODE: {e}")
    print()
    print("Falling back to direct inspection of result...")
    exit(1)

# Now test the ODE
print("=" * 80)
print("Step 4: Evaluate dH/dt at final state")
print("=" * 80)
print()

# Build state vector from result
y = np.array([result['all_densities'][sp] for sp in species])
H_idx = species.index('H')

print(f"State vector built:")
print(f"  H = {y[H_idx]:.2e} cm⁻³")
print(f"  H index = {H_idx}")
print()

print(f"ODE drift value:")
print(f"  H_drift_gain = {ode.H_drift_gain:.2e} cm⁻³/s")
print()

# Evaluate ODE at this state
dydt = ode(0, y)

print(f"ODE evaluation:")
print(f"  dH/dt = {dydt[H_idx]:.2e} cm⁻³/s")
print(f"  dH/dt / H = {dydt[H_idx] / y[H_idx]:.2e} s⁻¹")
print()

if abs(dydt[H_idx]) < 1e12:
    print("✓ System IS at steady state (dH/dt ≈ 0)")
else:
    print(f"✗ System NOT at steady state (dH/dt = {dydt[H_idx]:.2e})")

print()

# Now the critical test: remove drift and see the difference
print("=" * 80)
print("Step 5: TEST WITH AND WITHOUT DRIFT")
print("=" * 80)
print()

# Save original
drift_original = ode.H_drift_gain

# Test WITHOUT drift
ode.H_drift_gain = 0.0
dydt_no_drift = ode(0, y)
dH_dt_no_drift = dydt_no_drift[H_idx]

# Test WITH drift
ode.H_drift_gain = drift_original
dydt_with_drift = ode(0, y)
dH_dt_with_drift = dydt_with_drift[H_idx]

difference = dH_dt_with_drift - dH_dt_no_drift

print(f"WITHOUT drift (H_drift_gain = 0):")
print(f"  dH/dt = {dH_dt_no_drift:.2e} cm⁻³/s")
print()

print(f"WITH drift (H_drift_gain = {drift_original:.2e}):")
print(f"  dH/dt = {dH_dt_with_drift:.2e} cm⁻³/s")
print()

print(f"DIFFERENCE:")
print(f"  dH/dt(with) - dH/dt(without) = {difference:.2e} cm⁻³/s")
print(f"  Expected difference = {drift_original:.2e} cm⁻³/s")
print(f"  Ratio = {difference / drift_original:.6f}")
print()

# Check if they match
if abs(difference - drift_original) / drift_original < 0.01:
    print("✓✓✓ DRIFT TERM IS BEING APPLIED!")
    print()
    print("The drift term is working correctly in the code.")
    print()
    print("Since dH/dt ≈ 0 at steady state, this means:")
    print(f"  Drift source: {drift_original:.2e} cm⁻³/s")
    print(f"  Must be balanced by sinks of: {drift_original:.2e} cm⁻³/s")
    print()
    print("From balance analysis we saw:")
    print(f"  Chemical production: 2.07e15 cm⁻³/s")
    print(f"  Wall sticking: 1.91e15 cm⁻³/s")
    print(f"  Other chemical sinks: ~0.09e15 cm⁻³/s")
    print(f"  Total chemical: ~2.00e15 cm⁻³/s")
    print()
    print("So where is the drift sink?")
    print(f"  Drift: 7.74e16 cm⁻³/s")
    print(f"  Chemical: 2.07e15 cm⁻³/s")
    print(f"  Total production: 7.95e16 cm⁻³/s")
    print()
    print("This should give H = 7.95e16 / 138 s⁻¹ = 5.8e14 cm⁻³")
    print(f"But actual H = {y[H_idx]:.2e} cm⁻³ (40× lower!)")
    print()
    print("MYSTERY: Where is the missing H sink?")
else:
    print("✗✗✗ DRIFT TERM IS NOT BEING APPLIED!")
    print()
    print(f"Expected difference: {drift_original:.2e}")
    print(f"Actual difference:   {difference:.2e}")
    print(f"Ratio: {difference/drift_original:.3f}")
    print()
    print("There is a BUG in the drift implementation!")
