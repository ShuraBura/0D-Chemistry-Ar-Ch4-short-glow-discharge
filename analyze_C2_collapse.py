#!/usr/bin/env python3
"""
Analyze C2 balance at TRUE steady state for best_f49.5
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from build_reactions import build_reactions

print("=" * 80)
print("C2 BALANCE ANALYSIS AT TRUE STEADY STATE")
print("=" * 80)
print()

# Load best result
with open('optimization_results_charge_balanced/best_f49.5.json', 'r') as f:
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
C2_idx = species.index('C2')
CH_idx = species.index('CH')

y0[H_idx] = 1e11

sol = solve_ivp(ode, (0, 100), y0, method='BDF',
               rtol=1e-7, atol=1e-9, max_step=0.5)

y_ss = sol.y[:, -1]

print(f"TRUE Steady State:")
print(f"  H  = {y_ss[H_idx]:.2e} cm⁻³ (target: 2.52e14, {y_ss[H_idx]/2.52e14*100:.1f}%)")
print(f"  CH = {y_ss[CH_idx]:.2e} cm⁻³ (target: 1.00e9, {y_ss[CH_idx]/1e9:.1f}%)")
print(f"  C2 = {y_ss[C2_idx]:.2e} cm⁻³ (target: 5.60e11, {y_ss[C2_idx]/5.6e11*100:.3f}%)")
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
    reactants = rxn.reactants
    products = rxn.products

    rate_val = k.get(tag, 0.0)
    if rate_val == 0.0:
        continue

    # Compute flux
    flux = rate_val
    for sp in reactants:
        flux *= n.get(sp, 0.0)

    # C2 stoichiometry
    C2_in = list(reactants).count('C2')
    C2_out = list(products).count('C2')
    C2_net = C2_out - C2_in

    if C2_net > 0:
        C2_production.append({
            'reaction': f"{'+'.join(list(reactants))} → {'+'.join(list(products))}",
            'tag': tag,
            'rate': rate_val,
            'flux': flux,
            'C2_produced': flux * C2_net
        })
    elif C2_net < 0:
        C2_consumption.append({
            'reaction': f"{'+'.join(list(reactants))} → {'+'.join(list(products))}",
            'tag': tag,
            'rate': rate_val,
            'flux': flux,
            'C2_consumed': flux * abs(C2_net)
        })

C2_production.sort(key=lambda x: x['C2_produced'], reverse=True)
C2_consumption.sort(key=lambda x: x['C2_consumed'], reverse=True)

print("=" * 80)
print("C2 PRODUCTION")
print("=" * 80)
total_prod = sum(r['C2_produced'] for r in C2_production)

for i, r in enumerate(C2_production[:15]):
    pct = r['C2_produced'] / total_prod * 100 if total_prod > 0 else 0
    print(f"{i+1}. {r['reaction']}")
    print(f"   Rate: {r['C2_produced']:.2e} cm⁻³/s ({pct:.1f}%)")

print()
print(f"TOTAL C2 PRODUCTION: {total_prod:.2e} cm⁻³/s")
print()

print("=" * 80)
print("C2 CONSUMPTION")
print("=" * 80)
total_cons = sum(r['C2_consumed'] for r in C2_consumption)

for i, r in enumerate(C2_consumption[:15]):
    pct = r['C2_consumed'] / total_cons * 100 if total_cons > 0 else 0
    print(f"{i+1}. {r['reaction']}")
    print(f"   Rate: {r['C2_consumed']:.2e} cm⁻³/s ({pct:.1f}%)")
    if pct > 20:
        print(f"   >>> MAJOR SINK! <<<")

print()
print(f"TOTAL C2 CONSUMPTION: {total_cons:.2e} cm⁻³/s")
print()

# Wall loss
if 'stick_C2_9_9' in result['rate_values']:
    k_wall = result['rate_values']['stick_C2_9_9']
    C2_wall = k_wall * y_ss[C2_idx]
    print(f"C2 wall loss: {C2_wall:.2e} cm⁻³/s ({C2_wall/total_cons*100:.1f}% of total sinks)")
    print()

print("=" * 80)
print("BALANCE")
print("=" * 80)
print(f"Production:  {total_prod:.2e} cm⁻³/s")
print(f"Consumption: {total_cons:.2e} cm⁻³/s")
print(f"Net:         {total_prod - total_cons:.2e} cm⁻³/s")
print()

if total_cons > 0:
    print(f"P/C ratio: {total_prod/total_cons:.3f}")
    if total_prod / total_cons < 0.1:
        print(">>> C2 is SINK-LIMITED - production is much lower than consumption capacity!")
    elif total_prod / total_cons > 10:
        print(">>> C2 production is HUGE but consumption is weak")

print()
print("=" * 80)
print("DIAGNOSIS")
print("=" * 80)
print()

# Compare with target
print(f"C2 current: {y_ss[C2_idx]:.2e} cm⁻³")
print(f"C2 target:  5.60e11 cm⁻³")
print(f"Shortfall:  {5.6e11/y_ss[C2_idx]:.0f}× too low")
print()

# What production rate would we need?
C2_target = 5.6e11
needed_prod = C2_target * (total_cons / y_ss[C2_idx]) if y_ss[C2_idx] > 0 else 0

print(f"At current sink rate ({total_cons:.2e} cm⁻³/s),")
print(f"we need production rate: {needed_prod:.2e} cm⁻³/s")
print(f"Current production:      {total_prod:.2e} cm⁻³/s")
print(f"Need to increase production by: {needed_prod/total_prod:.0f}× ")
