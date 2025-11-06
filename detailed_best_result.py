#!/usr/bin/env python3
"""
Get detailed breakdown of ALL species for the best result
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
db = get_complete_rate_database()

# BEST conditions
TE_BEST = 1.20
E_BEST = 75
NE_BEST = 1.0e8

# Target densities
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

print("="*80)
print("DETAILED BREAKDOWN: BEST RESULT")
print("="*80)
print()
print("PLASMA CONDITIONS:")
print(f"  Te = {TE_BEST} eV")
print(f"  E_field = {E_BEST} V/cm")
print(f"  ne = {NE_BEST:.2e} cm⁻³")
print()

# Setup simulation with ALL optimizations + CH consumption maximized
params = params_base.copy()
params['rate_values'] = params_base['rate_values'].copy()

# All optimizations
params['rate_values']['stick_H_9_1'] = 3.89e2
params['rate_values']['stick_C2_9_9'] = 1.25e3
params['rate_values']['loss_C2_11_3'] = 100
params['rate_values']['C2_H_CH_C_cm3_7_6'] = 8.0e-11
params['rate_values']['stick_C2H2_9_11'] = 500
params['rate_values']['loss_C2H2_11_19'] = 1000
params['rate_values']['loss_C2H2Star_11_25'] = 100
params['rate_values']['CH_CH2_C2H2_H_cm3_7_7'] = 1.20e-10
params['rate_values']['CH2_CH2_C2H2_H2_cm3_7_15'] = 1.20e-10
params['rate_values']['CH3_CH_C2H2_H2_cm3_7_16'] = 1.20e-10
params['rate_values']['CH2_C_C2H2_cm3_7_17'] = 1.20e-10
params['rate_values']['CH_C2H2_C3H2_H_cm3_7_22'] = 8.0e-11
params['rate_values']['CH_C2H2_C3H_H2_cm3_7_27'] = 8.0e-11
params['rate_values']['CH_C2H2_C2H_CH2_cm3_7_29'] = 8.0e-11
params['rate_values']['e_CH4_CH_H_H2_cm3_1_11'] = 2.0e-11
params['rate_values']['e_CH4_CH_H2_H_vib_cm3_1_3'] = 2.0e-11

# CH consumption maximized
params['rate_values']['loss_CH_11_9'] = 1.00e4
params['rate_values']['stick_CH_9_3'] = 6.25e3
params['rate_values']['CH_CH4_CH2_CH3_cm3_7_39'] = 1.20e-11
params['rate_values']['CH_H_CH2_cm3_7_21'] = 1.20e-10
params['rate_values']['CH_CH3_C2H3_H_cm3_7_10'] = 1.20e-10
params['rate_values']['CH_CH3_C2H2_H2_cm3_7_23'] = 1.20e-10
params['rate_values']['CH_H2_CH2_H_cm3_7_30'] = 1.20e-11
params['rate_values']['CH_H2_CH2_H_cm3_7_51'] = 1.20e-11
params['rate_values']['C2H2_CH_C3_H2_H_cm3_7_56'] = 1.20e-10
params['rate_values']['CH_C2H2_C3H2_H_cm3_7_22'] = 1.20e-10

params['Te'] = TE_BEST
params['E_field'] = E_BEST
params['ne'] = NE_BEST

# Get rate dictionary
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
print("Running simulation...")
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
    print("SIMULATION FAILED!")
    exit(1)

print("Simulation complete!")
print()

# Extract final densities
y_final = sol.y[:, -1]

# Calculate total ion density
ni_total = 0
ni_minus = 0
for i_sp, species_name in enumerate(ode.species):
    if species_name.endswith('Plus'):
        ni_total += y_final[i_sp]
    elif species_name.endswith('Minus'):
        ni_minus += y_final[i_sp]
        ni_total -= y_final[i_sp]

ne = y_final[ode.species.index('e')]
ni_ne_ratio = ni_total / ne

print("="*80)
print("NEUTRALS (cm⁻³)")
print("="*80)
print()

neutral_species = []
for i_sp, species_name in enumerate(ode.species):
    if not species_name.endswith('Plus') and not species_name.endswith('Minus') and species_name != 'e':
        density = y_final[i_sp]
        if density > 1e5:  # Only show significant densities
            neutral_species.append((species_name, density))

# Sort by density
neutral_species.sort(key=lambda x: x[1], reverse=True)

for sp_name, density in neutral_species:
    # Check if it's a target species
    if sp_name in TARGETS:
        ratio = density / TARGETS[sp_name]
        status = "✓" if 0.6 <= ratio <= 1.4 else "✗"
        print(f"  {sp_name:12s} = {density:12.3e}  ({ratio:6.2f}× target) {status}")
    else:
        print(f"  {sp_name:12s} = {density:12.3e}")

print()
print("="*80)
print("POSITIVE IONS (cm⁻³)")
print("="*80)
print()

positive_ions = []
for i_sp, species_name in enumerate(ode.species):
    if species_name.endswith('Plus'):
        density = y_final[i_sp]
        if density > 1e3:  # Only show significant densities
            positive_ions.append((species_name, density))

# Sort by density
positive_ions.sort(key=lambda x: x[1], reverse=True)

for sp_name, density in positive_ions:
    fraction = density / (ni_total + ni_minus) * 100
    print(f"  {sp_name:12s} = {density:12.3e}  ({fraction:5.1f}% of ions)")

print()
print(f"  TOTAL (ni+) = {ni_total + ni_minus:12.3e}")
print()

print("="*80)
print("NEGATIVE IONS (cm⁻³)")
print("="*80)
print()

negative_ions = []
for i_sp, species_name in enumerate(ode.species):
    if species_name.endswith('Minus'):
        density = y_final[i_sp]
        if density > 1e3:  # Only show significant densities
            negative_ions.append((species_name, density))

# Sort by density
negative_ions.sort(key=lambda x: x[1], reverse=True)

if negative_ions:
    for sp_name, density in negative_ions:
        print(f"  {sp_name:12s} = {density:12.3e}")
    print()
    print(f"  TOTAL (ni-) = {ni_minus:12.3e}")
else:
    print("  None detected")

print()
print("="*80)
print("CHARGE BALANCE")
print("="*80)
print()
print(f"  ne          = {ne:12.3e} cm⁻³")
print(f"  ni (total)  = {ni_total:12.3e} cm⁻³")
print(f"  ni/ne ratio = {ni_ne_ratio:12.3f}  {'✓ IN RANGE' if 2.0 <= ni_ne_ratio <= 6.0 else '✗ OUT OF RANGE'}")
print()

print("="*80)
print("TARGET SPECIES SUMMARY")
print("="*80)
print()

for sp_name, target in TARGETS.items():
    idx = ode.species.index(sp_name)
    density = y_final[idx]
    ratio = density / target

    if ratio < 0.6:
        status = f"✗ TOO LOW (need {0.6*target:.2e})"
        err = (0.6 / ratio) ** 2
    elif ratio > 1.4:
        status = f"✗ TOO HIGH (need {1.4*target:.2e})"
        err = (ratio / 1.4) ** 2
    else:
        status = "✓ IN RANGE"
        err = 1.0

    print(f"  {sp_name:12s} = {density:12.3e} cm⁻³  ({ratio:6.2f}× target)")
    print(f"  {'':12s}   Target: {target:12.3e} cm⁻³  Range: [{0.6*target:.2e}, {1.4*target:.2e}]")
    print(f"  {'':12s}   {status}  (penalty: {err:.2f})")
    print()

# Calculate objective
errors = []
for sp, target in TARGETS.items():
    idx = ode.species.index(sp)
    density = y_final[idx]
    ratio = density / target
    if ratio < 0.6:
        err = (0.6 / ratio) ** 2
    elif ratio > 1.4:
        err = (ratio / 1.4) ** 2
    else:
        err = 1.0
    errors.append(err)

objective = np.prod(errors)

print("="*80)
print("OBJECTIVE FUNCTION")
print("="*80)
print()
print(f"  f(x) = {objective:.2f}")
print(f"  Individual penalties: {errors}")
print(f"  Status: {'✓ ALL IN RANGE' if objective == 1.0 else '✗ SOME OUT OF RANGE'}")
print()
print("="*80)

# Summary table
print()
print("="*80)
print("COMPLETE SUMMARY")
print("="*80)
print()
print("Plasma Conditions:")
print(f"  Te       = {TE_BEST} eV")
print(f"  E_field  = {E_BEST} V/cm")
print(f"  ne       = {NE_BEST:.2e} cm⁻³")
print()
print("Target Species:")
for sp_name, target in TARGETS.items():
    idx = ode.species.index(sp_name)
    density = y_final[idx]
    ratio = density / target
    status = "✓" if 0.6 <= ratio <= 1.4 else "✗"
    print(f"  {sp_name:4s} = {density:12.3e} cm⁻³  ({ratio:6.2f}×)  {status}")
print()
print("Charge Balance:")
print(f"  Ni/Ne = {ni_ne_ratio:.3f}  {'✓' if 2.0 <= ni_ne_ratio <= 6.0 else '✗'}")
print()
print("Objective:")
print(f"  f(x) = {objective:.2f}")
print()
print("="*80)
