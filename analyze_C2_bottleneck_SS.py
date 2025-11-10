#!/usr/bin/env python3
"""
Analyze C2 bottleneck at TRUE steady state, including C2H2 density
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from build_reactions import build_reactions

print("=" * 80)
print("C2 BOTTLENECK ANALYSIS AT TRUE STEADY STATE")
print("=" * 80)
print()

# Load best result
with open('optimization_results_charge_balanced/best_f31.3.json', 'r') as f:
    result = json.load(f)

# Setup ODE
species = ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
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
    'E_field': result['E_field'], 'L_discharge': 0.45,
    'ne': result['Ne'], 'Te': result['Te'], 'species': species,
    'T': 400.0, 'Tgas': 400.0, 'pressure': 500.0, 'mobilities': mobilities
}

k = define_rates(params)
for name, val in result['rate_values'].items():
    if name in k:
        k[name] = val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)
ode = PlasmaODE_Optimized(params)

# Integrate to TRUE steady state
y0 = np.array([result['all_densities'][sp] for sp in species])
H_idx = species.index('H')
y0[H_idx] = 1e11

sol = solve_ivp(ode, (0, 100), y0, method='BDF',
               rtol=1e-7, atol=1e-9, max_step=0.5)

y_ss = sol.y[:, -1]

# Extract key species
C2_idx = species.index('C2')
C2H2_idx = species.index('C2H2')
C2H_idx = species.index('C2H')
CH_idx = species.index('CH')
C_idx = species.index('C')
C2HPlus_idx = species.index('C2HPlus')
e_idx = species.index('e')

print("STEADY STATE DENSITIES:")
print(f"  C2   = {y_ss[C2_idx]:.2e} cm⁻³ (target: 5.60e11, {y_ss[C2_idx]/5.6e11*100:.3f}%)")
print(f"  C2H2 = {y_ss[C2H2_idx]:.2e} cm⁻³")
print(f"  C2H  = {y_ss[C2H_idx]:.2e} cm⁻³")
print(f"  CH   = {y_ss[CH_idx]:.2e} cm⁻³")
print(f"  C    = {y_ss[C_idx]:.2e} cm⁻³")
print(f"  C2H+ = {y_ss[C2HPlus_idx]:.2e} cm⁻³")
print(f"  H    = {y_ss[H_idx]:.2e} cm⁻³")
print(f"  e    = {y_ss[e_idx]:.2e} cm⁻³")
print()

c2h2_ratio = y_ss[C2H2_idx] / y_ss[C2_idx]
print(f"C2H2/C2 ratio: {c2h2_ratio:.1f}×")
print()

# Build densities dict
n = {sp: y_ss[species.index(sp)] for sp in species}

# Analyze C2 production and consumption
R = params['R']
tags = params['tags']

C2_production = []
C2_consumption = []

for i, rxn in enumerate(R):
    tag = tags[i]
    reactants = list(rxn.reactants)
    products = list(rxn.products)

    rate_val = k.get(tag, 0.0)
    if rate_val == 0.0:
        continue

    # Compute flux
    flux = rate_val
    for sp in reactants:
        flux *= n.get(sp, 0.0)

    # C2 stoichiometry
    C2_in = reactants.count('C2')
    C2_out = products.count('C2')
    C2_net = C2_out - C2_in

    if C2_net > 0:
        C2_production.append({
            'reaction': f"{'+'.join(reactants)} → {'+'.join(products)}",
            'tag': tag,
            'rate': rate_val,
            'flux': flux,
            'C2_produced': flux * C2_net
        })
    elif C2_net < 0:
        C2_consumption.append({
            'reaction': f"{'+'.join(reactants)} → {'+'.join(products)}",
            'tag': tag,
            'rate': rate_val,
            'flux': flux,
            'C2_consumed': flux * abs(C2_net)
        })

C2_production.sort(key=lambda x: x['C2_produced'], reverse=True)
C2_consumption.sort(key=lambda x: x['C2_consumed'], reverse=True)

print("=" * 80)
print("C2 PRODUCTION PATHWAYS (TRUE STEADY STATE)")
print("=" * 80)
total_prod = sum(r['C2_produced'] for r in C2_production)

if total_prod > 0:
    for i, r in enumerate(C2_production[:15]):
        pct = r['C2_produced'] / total_prod * 100
        print(f"{i+1}. {r['reaction']}")
        print(f"   Production: {r['C2_produced']:.2e} cm⁻³/s ({pct:.1f}%)")
        print(f"   k = {r['rate']:.2e}, flux = {r['flux']:.2e}")
else:
    print("NO C2 PRODUCTION!")

print()
print(f"TOTAL C2 PRODUCTION: {total_prod:.2e} cm⁻³/s")
print()

print("=" * 80)
print("C2 CONSUMPTION PATHWAYS (TRUE STEADY STATE)")
print("=" * 80)
total_cons = sum(r['C2_consumed'] for r in C2_consumption)

if total_cons > 0:
    for i, r in enumerate(C2_consumption[:15]):
        pct = r['C2_consumed'] / total_cons * 100
        print(f"{i+1}. {r['reaction']}")
        print(f"   Consumption: {r['C2_consumed']:.2e} cm⁻³/s ({pct:.1f}%)")
        if pct > 20:
            print(f"   >>> MAJOR SINK! <<<")
else:
    print("NO C2 CONSUMPTION!")

print()
print(f"TOTAL C2 CHEMICAL CONSUMPTION: {total_cons:.2e} cm⁻³/s")
print()

# Wall loss
if 'stick_C2_9_9' in result['rate_values']:
    k_wall = result['rate_values']['stick_C2_9_9']
    C2_wall = k_wall * y_ss[C2_idx]
    print(f"C2 wall loss: {C2_wall:.2e} cm⁻³/s")
    print()

    total_loss = total_cons + C2_wall
    print(f"TOTAL C2 LOSS: {total_loss:.2e} cm⁻³/s")
    print()

    if C2_wall > 0:
        print(f"Wall loss fraction: {C2_wall/total_loss*100:.1f}%")
        print(f"Chemical loss fraction: {total_cons/total_loss*100:.1f}%")
        print()

# Check key C2-producing reactions specifically
print("=" * 80)
print("KEY C2-PRODUCING REACTIONS - RATES AND FLUXES")
print("=" * 80)
print()

key_reactions = {
    'CH_CH_C2_H2_cm3_5_4': ('CH + CH → C2 + H2', n.get('CH', 0)**2),
    'C_CH_C2_H_cm3_7_4': ('C + CH → C2 + H', n.get('C', 0) * n.get('CH', 0)),
    'CH_C_C2_H_cm3_7_9': ('CH + C → C2 + H', n.get('CH', 0) * n.get('C', 0)),
    'e_C2H2_C2_H2_cm3_1_16': ('e + C2H2 → C2 + H2', n.get('e', 0) * n.get('C2H2', 0)),
    'C2HPlus_e_C2_H_cm3_6_18': ('C2H+ + e → C2 + H', n.get('C2HPlus', 0) * n.get('e', 0)),
    'C2H_H_C2_H2_cm3_7_47': ('C2H + H → C2 + H2', n.get('C2H', 0) * n.get('H', 0)),
}

for tag, (desc, density_product) in key_reactions.items():
    rate = k.get(tag, 0.0)
    flux = rate * density_product

    print(f"{desc}")
    print(f"  k = {rate:.2e} cm³/s")
    print(f"  flux = {flux:.2e} cm⁻³/s")

    if flux < 1e8:
        print(f"  >>> NEGLIGIBLE! <<<")
    elif flux > 1e12:
        print(f"  >>> SIGNIFICANT! <<<")
    print()

print("=" * 80)
print("DIAGNOSIS")
print("=" * 80)
print()

if y_ss[C2H2_idx] > 1e13:
    print(f"✓ C2H2 is HIGH ({y_ss[C2H2_idx]:.2e} cm⁻³)")
    print(f"  But e+C2H2→C2+H2 flux is only {k.get('e_C2H2_C2_H2_cm3_1_16', 0) * n.get('e', 0) * n.get('C2H2', 0):.2e}")
    print()
    print("  Possible issues:")
    print("  1. e+C2H2→C2+H2 rate constant too low?")
    print("  2. C2H2 being consumed by other pathways (not producing C2)?")
    print("  3. C2 sinks too strong?")
else:
    print(f"C2H2 is low ({y_ss[C2H2_idx]:.2e} cm⁻³) - not enough precursor")

print()
print(f"C2 balance:")
print(f"  Production:  {total_prod:.2e} cm⁻³/s")
print(f"  Loss:        {total_loss:.2e} cm⁻³/s")
print(f"  P/L ratio:   {total_prod/total_loss if total_loss > 0 else 0:.3e}")
print()

if total_prod > 0 and total_loss > 0:
    expected_C2 = total_prod / (total_loss / y_ss[C2_idx])
    print(f"Expected C2 at this P/L: {expected_C2:.2e} cm⁻³")
    print(f"Actual C2:               {y_ss[C2_idx]:.2e} cm⁻³")
    print(f"Match: {abs(expected_C2 - y_ss[C2_idx])/y_ss[C2_idx] < 0.1}")
