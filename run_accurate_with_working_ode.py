#!/usr/bin/env python3
"""
Run high-accuracy simulation using the WORKING odefun.py class
"""

import json
import numpy as np
from scipy.integrate import solve_ivp
from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun import PlasmaODE  # Use the WORKING ODE class

print("="*80)
print("HIGH-ACCURACY SIMULATION WITH WORKING ODE CLASS")
print("="*80)
print()

# Load best checkpoint
with open('checkpoint_f3407.json', 'r') as f:
    best_result = json.load(f)

params_dict = best_result['params']
species = params_dict['species']

print(f"Using parameters from checkpoint_f3407.json")
print(f"Objective: f(x) = {best_result['objective']:.2f}")
print()
print("Plasma parameters:")
print(f"  Te: {params_dict['Te']:.3f} eV")
print(f"  Ne: {params_dict['ne']:.3e} cm^-3")
print(f"  E:  {params_dict['E_field']:.2f} V/cm")
print(f"  H drift gain: {params_dict['H_drift_gain']:.3e} cm^-3/s")
print()

# Build reactions
db = get_complete_rate_database()
k = define_rates_tunable(params_dict)

for name, val in params_dict.get('rate_values', {}).items():
    if name in k and name in db:
        k[name] = np.clip(val, db[name].min, db[name].max)

params = params_dict.copy()
params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Initial conditions
y0 = np.ones(len(species)) * 1e3

def set_density(name, value):
    try:
        y0[species.index(name)] = value
    except ValueError:
        pass

set_density('e', params['ne'])
set_density('Ar', 0.85 * 9.66e15)
set_density('CH4', 0.15 * 9.66e15)
set_density('ArPlus', 1e7)
set_density('CH4Plus', 1e5)
set_density('H2', 1e12)
set_density('H', 1e11)
set_density('C2', 5e7)
set_density('CH', 5e4)
set_density('CH3', 5e7)

print("Running HIGH-ACCURACY simulation with WORKING odefun.py class:")
print("  Method: BDF (implicit, stiff)")
print("  Tolerances: rtol=1e-9, atol=1e-12")
print("  Integration time: 0 to 1000 ms")
print("  Max step: 10 ms")
print()
print("This may take 2-5 minutes...")
print()

# HIGH-ACCURACY solver settings
sol = solve_ivp(
    PlasmaODE(params),  # Use WORKING class from odefun.py
    (0, 1000),
    y0,
    method='BDF',
    rtol=1e-9,
    atol=1e-12,
    max_step=10.0
)

if not sol.success:
    print(f"ERROR: Simulation failed - {sol.message}")
    exit(1)

print(f"SUCCESS! Simulation converged in {sol.nfev} function evaluations")
print()

y_final = sol.y[:, -1]

def get_density(name):
    try:
        return y_final[species.index(name)]
    except ValueError:
        return 0.0

# Calculate ion density and charge balance
ion_species = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
               'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus',
               'H3Plus', 'CHPlus', 'H2Plus']

print("="*80)
print("ALL SPECIES DENSITIES (cm^-3)")
print("="*80)
print()

# Group species by type
target_densities = {}
ion_densities = {}
neutral_species = []

for sp in species:
    dens = get_density(sp)

    if sp in ['H', 'CH', 'C2', 'C2H2']:
        target_densities[sp] = dens
    elif sp in ion_species:
        ion_densities[sp] = dens
    else:
        neutral_species.append((sp, dens))

# Print target species first
print("TARGET SPECIES:")
for sp in ['H', 'CH', 'C2', 'C2H2']:
    if sp in target_densities:
        print(f"  {sp:12s}: {target_densities[sp]:12.3e}")
print()

# Print ions
print("ION SPECIES:")
all_ions_positive = True
for sp in ion_species:
    if sp in ion_densities:
        dens = ion_densities[sp]
        status = ""
        if dens < 0:
            status = " ✗ NEGATIVE!"
            all_ions_positive = False
        print(f"  {sp:12s}: {dens:12.3e}{status}")
print()

if all_ions_positive:
    print("✓ All ions are POSITIVE!")
else:
    print("✗ WARNING: Some ions are negative - unphysical!")
print()

# Print neutrals
print("NEUTRAL SPECIES:")
for sp, dens in neutral_species:
    marker = ""
    if sp == 'e':
        marker = " ← Ne (electron density)"
    print(f"  {sp:12s}: {dens:12.3e}{marker}")

print()
print("="*80)
print("SUMMARY")
print("="*80)
print()

n_i = sum(ion_densities.values())
n_e = get_density('e')
charge_imbalance = abs(n_i - n_e) / n_e * 100 if n_e > 0 else 999

print(f"Total ion density (Ni):     {n_i:.3e} cm^-3")
print(f"Electron density (Ne):      {n_e:.3e} cm^-3")
print(f"Charge imbalance:           {charge_imbalance:.2f}%")
print()
print(f"Electric field (E):         {params_dict['E_field']:.2f} V/cm")
print(f"Electron temp (Te):         {params_dict['Te']:.3f} eV")
print(f"H drift gain:               {params_dict['H_drift_gain']:.3e} cm^-3/s")
print(f"Pressure (P):               {params_dict['P']} Torr")
print(f"Gas temperature (Tgas):     {params_dict['Tgas']} K")
print()

# Target comparison
TARGETS = {'H': 5.18e13, 'CH': 1.0e9, 'C2': 1.3e11}
print("TARGET SPECIES COMPARISON:")
print(f"  H:    {target_densities['H']:.3e} cm^-3  ({target_densities['H']/TARGETS['H']:.2f}× target, goal: 0.6-1.4×)")
print(f"  CH:   {target_densities['CH']:.3e} cm^-3  ({target_densities['CH']/TARGETS['CH']:.2f}× target, goal: 0.6-1.4×)")
print(f"  C2:   {target_densities['C2']:.3e} cm^-3  ({target_densities['C2']/TARGETS['C2']:.2f}× target, goal: 0.6-1.4×)")
print(f"  C2H2: {target_densities['C2H2']:.3e} cm^-3")
print()

# Status
h_ok = 0.6 <= target_densities['H']/TARGETS['H'] <= 1.4
ch_ok = 0.6 <= target_densities['CH']/TARGETS['CH'] <= 1.4
c2_ok = 0.6 <= target_densities['C2']/TARGETS['C2'] <= 1.4

print("STATUS:")
print(f"  H:  {'✓ IN RANGE' if h_ok else '✗ OUT OF RANGE'}")
print(f"  CH: {'✓ IN RANGE' if ch_ok else '✗ OUT OF RANGE'}")
print(f"  C2: {'✓ IN RANGE' if c2_ok else '✗ OUT OF RANGE'}")

if charge_imbalance < 5:
    print(f"  Charge balance: ✓ GOOD ({charge_imbalance:.2f}% imbalance)")
elif charge_imbalance < 20:
    print(f"  Charge balance: ~ ACCEPTABLE ({charge_imbalance:.2f}% imbalance)")
else:
    print(f"  Charge balance: ✗ POOR ({charge_imbalance:.2f}% imbalance)")

print()
print("="*80)

# Save to file
output = {
    'objective': best_result['objective'],
    'parameters': {
        'Te': params_dict['Te'],
        'Ne': params_dict['ne'],
        'E_field': params_dict['E_field'],
        'H_drift_gain': params_dict['H_drift_gain'],
        'P': params_dict['P'],
        'Tgas': params_dict['Tgas']
    },
    'all_densities': {sp: get_density(sp) for sp in species},
    'ion_densities': ion_densities,
    'target_densities': target_densities,
    'total_ion_density': n_i,
    'electron_density': n_e,
    'charge_imbalance_percent': charge_imbalance
}

with open('final_accurate_result.json', 'w') as f:
    json.dump(output, f, indent=2)

print("Results saved to: final_accurate_result.json")
print()
