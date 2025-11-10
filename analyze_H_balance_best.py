#!/usr/bin/env python3
"""
Compute the H atom balance for the best result to see where H is going.
"""

import numpy as np
import json
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from build_reactions import build_reactions

print("=" * 80)
print("H BALANCE ANALYSIS - BEST RESULT (f=41.3)")
print("=" * 80)
print()

# Load result
with open('optimization_results_charge_balanced/best_f41.3.json', 'r') as f:
    result = json.load(f)

# Setup
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

print("Building ODE with optimized rates...")
k = define_rates(params)

# Override with optimized rates
for name, val in result['rate_values'].items():
    if name in k:
        k[name] = val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)
ode = PlasmaODE_Optimized(params)

print(f"✓ ODE created")
print(f"  H_drift_gain = {ode.H_drift_gain:.2e} cm⁻³/s")
print()

# Build state vector
y = np.array([result['all_densities'][sp] for sp in species_list])
H_idx = species_list.index('H')

print(f"State:")
print(f"  H = {y[H_idx]:.2e} cm⁻³")
print(f"  Ne = {result['Ne']:.2e} cm⁻³")
print(f"  Te = {result['Te']:.3f} eV")
print()

# Evaluate dH/dt
dydt = ode(0, y)

print(f"dH/dt = {dydt[H_idx]:.2e} cm⁻³/s")

if abs(dydt[H_idx]) < 1e12:
    print("✓ At steady state")
else:
    print(f"✗ NOT at steady state (relative rate: {dydt[H_idx]/y[H_idx]:.2e} s⁻¹)")

print()
print("=" * 80)
print("H PRODUCTION AND CONSUMPTION BREAKDOWN")
print("=" * 80)
print()

# Compute reaction rates manually
from build_reactions import build_reactions

R = params['R']
tags = params['tags']

# Get densities dict
n = {sp: result['all_densities'][sp] for sp in species_list}

# H production and consumption
H_production = []
H_consumption = []

for i, (reactants, products, rate_const) in enumerate(R):
    tag = tags[i]

    # Compute rate
    rate_val = k.get(tag, 0.0)
    if rate_val == 0.0:
        continue

    # Compute reaction flux
    flux = rate_val
    for sp in reactants:
        flux *= n.get(sp, 0.0)

    # Check H stoichiometry
    H_in = reactants.count('H')
    H_out = products.count('H')
    H_net = H_out - H_in

    if H_net > 0:
        # H production
        H_production.append({
            'reaction': f"{'+'.join(reactants)} → {'+'.join(products)}",
            'tag': tag,
            'rate': rate_val,
            'flux': flux,
            'H_produced': flux * H_net
        })
    elif H_net < 0:
        # H consumption
        H_consumption.append({
            'reaction': f"{'+'.join(reactants)} → {'+'.join(products)}",
            'tag': tag,
            'rate': rate_val,
            'flux': flux,
            'H_consumed': flux * abs(H_net)
        })

# Sort by importance
H_production.sort(key=lambda x: x['H_produced'], reverse=True)
H_consumption.sort(key=lambda x: x['H_consumed'], reverse=True)

print("Top H Production Pathways:")
total_prod = sum(r['H_produced'] for r in H_production)
for i, r in enumerate(H_production[:10]):
    pct = r['H_produced'] / total_prod * 100 if total_prod > 0 else 0
    print(f"{i+1}. {r['reaction']}")
    print(f"   Rate: {r['H_produced']:.2e} cm⁻³/s ({pct:.1f}%)")
    print(f"   k = {r['rate']:.2e}, flux = {r['flux']:.2e}")

print()
print(f"Total chemical H production: {total_prod:.2e} cm⁻³/s")
print(f"H drift source: {ode.H_drift_gain:.2e} cm⁻³/s")
print(f"Total H sources: {total_prod + ode.H_drift_gain:.2e} cm⁻³/s")
print()

print("Top H Consumption Pathways:")
total_cons = sum(r['H_consumed'] for r in H_consumption)
for i, r in enumerate(H_consumption[:10]):
    pct = r['H_consumed'] / total_cons * 100 if total_cons > 0 else 0
    print(f"{i+1}. {r['reaction']}")
    print(f"   Rate: {r['H_consumed']:.2e} cm⁻³/s ({pct:.1f}%)")
    print(f"   k = {r['rate']:.2e}, flux = {r['flux']:.2e}")

print()
print(f"Total chemical H consumption: {total_cons:.2e} cm⁻³/s")

# Wall loss
if 'stick_H_9_1' in result['rate_values']:
    k_wall = result['rate_values']['stick_H_9_1']
    H_wall_loss = k_wall * y[H_idx]
    print(f"H wall loss: {H_wall_loss:.2e} cm⁻³/s")
    print()

    print(f"Total H sinks: {total_cons + H_wall_loss:.2e} cm⁻³/s")
    print()

    print("=" * 80)
    print("BALANCE CHECK")
    print("=" * 80)
    print()
    print(f"Sources:  {total_prod + ode.H_drift_gain:.2e} cm⁻³/s")
    print(f"Sinks:    {total_cons + H_wall_loss:.2e} cm⁻³/s")
    print(f"Net:      {(total_prod + ode.H_drift_gain) - (total_cons + H_wall_loss):.2e} cm⁻³/s")
    print()

    imbalance = abs((total_prod + ode.H_drift_gain) - (total_cons + H_wall_loss))
    rel_imbalance = imbalance / (total_prod + ode.H_drift_gain) * 100

    if rel_imbalance < 1.0:
        print(f"✓ Balanced to {rel_imbalance:.2f}%")
    else:
        print(f"✗ Imbalance: {rel_imbalance:.1f}%")
