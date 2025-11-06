#!/usr/bin/env python3
"""
Check ion densities in our best optimized result to verify no negative densities
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from odefun import PlasmaODE
from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from rate_database_complete import get_complete_rate_database

# Load checkpoint
with open('checkpoint_f3407.json', 'r') as f:
    checkpoint = json.load(f)

params_base = checkpoint['params'].copy()

# Apply BEST parameters discovered
params = params_base.copy()
params['rate_values'] = params_base['rate_values'].copy()

params['rate_values']['stick_C2_9_9'] = 1.25e3      # lit min
params['rate_values']['loss_C2_11_3'] = 100         # near min
params['rate_values']['C2_H_CH_C_cm3_7_6'] = 8.0e-11  # lit min
params['rate_values']['stick_C2H2_9_11'] = 472      # 50% reduction

print("="*80)
print("ION DENSITY CHECK - Best optimized parameters")
print("="*80)
print()

# Get rate dictionary
db = get_complete_rate_database()
k_dict = define_rates_tunable(params)

# Apply rate_values overrides
for name, val in params.get('rate_values', {}).items():
    if name in k_dict and name in db:
        k_dict[name] = np.clip(val, db[name].min, db[name].max)

params['k'] = k_dict

# Build reactions
reactions, tags = build_reactions(params)
params['R'] = reactions
params['tags'] = tags

# Setup ODE system
ode = PlasmaODE(params)

# Initial conditions
y0 = np.zeros(ode.ns)
y0[ode.species.index('e')] = params['ne']
y0[ode.species.index('Ar')] = 1.29e16
y0[ode.species.index('CH4')] = 1.29e15

# Solve to steady-state
sol = solve_ivp(
    ode,
    t_span=[0, 1e-3],
    y0=y0,
    method='BDF',
    rtol=1e-6,
    atol=1e-10,
    max_step=1e-5
)

if not sol.success:
    print(f"ERROR: Simulation failed - {sol.message}")
    exit(1)

# Extract final densities
y_final = sol.y[:, -1]

# Identify ion species (end with Plus or Minus)
positive_ions = []
negative_ions = []
neutrals = []

for i, species_name in enumerate(ode.species):
    density = y_final[i]
    if species_name.endswith('Plus'):
        positive_ions.append((species_name, density))
    elif species_name.endswith('Minus'):
        negative_ions.append((species_name, density))
    elif species_name == 'e':
        pass  # Handle separately
    else:
        neutrals.append((species_name, density))

# Electron density
ne = y_final[ode.species.index('e')]

print("ELECTRON DENSITY:")
print(f"  e: {ne:.2e} cm⁻³")
print()

print("POSITIVE ION DENSITIES:")
total_positive = 0
has_negative = False
for name, density in sorted(positive_ions, key=lambda x: -x[1]):
    if density < 0:
        print(f"  {name:15s}: {density:.2e} cm⁻³  ⚠️  NEGATIVE!")
        has_negative = True
    elif density > 1e6:  # Only show significant ions
        print(f"  {name:15s}: {density:.2e} cm⁻³")
    total_positive += density

print(f"\n  Total positive ions: {total_positive:.2e} cm⁻³")
print()

print("NEGATIVE ION DENSITIES:")
total_negative = 0
for name, density in sorted(negative_ions, key=lambda x: -x[1]):
    if density < 0:
        print(f"  {name:15s}: {density:.2e} cm⁻³  ⚠️  NEGATIVE!")
        has_negative = True
    elif density > 1e6:  # Only show significant ions
        print(f"  {name:15s}: {density:.2e} cm⁻³")
    total_negative += density

print(f"\n  Total negative ions: {total_negative:.2e} cm⁻³")
print()

# Charge balance
total_ions = total_positive - total_negative
charge_imbalance = abs(ne - total_ions)
relative_imbalance = charge_imbalance / ne * 100

print("CHARGE BALANCE:")
print(f"  Electrons:         {ne:.2e} cm⁻³")
print(f"  Positive ions:     {total_positive:.2e} cm⁻³")
print(f"  Negative ions:     {total_negative:.2e} cm⁻³")
print(f"  Net positive ions: {total_ions:.2e} cm⁻³")
print(f"  Charge imbalance:  {charge_imbalance:.2e} cm⁻³ ({relative_imbalance:.2e}%)")
print()

print("="*80)
if has_negative:
    print("⚠️  WARNING: NEGATIVE ION DENSITIES DETECTED!")
    print("This indicates a problem with the ODE solver or chemistry.")
else:
    print("✓ ALL ION DENSITIES ARE POSITIVE - No negative density issues!")

print()
print("KEY NEUTRAL SPECIES (top 10):")
for name, density in sorted(neutrals, key=lambda x: -x[1])[:10]:
    print(f"  {name:10s}: {density:.2e} cm⁻³")
