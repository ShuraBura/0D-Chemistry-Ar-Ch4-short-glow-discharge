#!/usr/bin/env python3
"""
Analyze electron balance in the high-C2H2 state (Eval 2: C2H2=4.09e+12)
to understand why Ni/Ne = 215 (ions vastly outnumber electrons)
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

def pressure_to_density(pressure_mTorr, T_K=400):
    kB = 1.38064852e-23
    Torr_to_Pa = 133.322
    P_Pa = pressure_mTorr * 1e-3 * Torr_to_Pa
    n_m3 = P_Pa / (kB * T_K)
    return n_m3 * 1e-6

# Load high C2H2 result
with open('optimization_results_comprehensive_1e12/best_f13946912.1.json', 'r') as f:
    high = json.load(f)

print('=' * 80)
print('ELECTRON BALANCE ANALYSIS - HIGH C2H2 STATE')
print('=' * 80)
print()
print('This state has:')
print(f'  C2H2: {high["all_densities"]["C2H2"]:.2e} cm⁻³ (409% of 1e+12 target) ✓')
print(f'  C2:   {high["target_densities"]["C2"]:.2e} cm⁻³ (275% of target) ✓')
print(f'  Ne:   {high["Ne"]:.2e} cm⁻³')
print(f'  Ni:   {high["n_i_total"]:.2e} cm⁻³')
print(f'  Ni/Ne: {high["Ni_over_Ne"]:.1f} ✗ (should be 2-7)')
print()

# Setup params to run chemistry
params = {
    'P': 500.0,
    'Te': high['Te'],
    'ne': high['Ne'],
    'E_field': high['E_field'],
    'L_discharge': 0.45,
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
    },
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar',
                'C2H4', 'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C',
                'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H',
                'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus',
                'H2Plus', 'C2H2Star'],
}

k = define_rates(params)

# Apply tuned rates
for rate_name, rate_val in high.get('rate_values', {}).items():
    if rate_name in k:
        k[rate_name] = rate_val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)

species = params['species']
y0 = np.array([high['all_densities'].get(s, 1e3) for s in species])

# Evaluate chemistry at this state
ode_func = PlasmaODE_Optimized(params)
dydt = ode_func(0, y0)

idx_e = species.index('e')
de_dt = dydt[idx_e]

print(f'Net electron rate: de/dt = {de_dt:.2e} cm⁻³/s')
print()

# Now analyze all reactions affecting electrons
R = params['R']
tags = params['tags']

electron_production = []  # Reactions that create electrons
electron_loss = []        # Reactions that destroy electrons

for i, rxn in enumerate(R):
    tag = tags[i]

    # Get net electron change
    net_e = rxn.products[idx_e] - rxn.reactants[idx_e]

    if net_e != 0:
        # Calculate reaction rate
        rate_coeff = rxn.rate
        reactant_density = 1.0
        for j, stoich in enumerate(rxn.reactants):
            if stoich > 0:
                reactant_density *= y0[j]**stoich

        reaction_rate = rate_coeff * reactant_density
        electron_rate = reaction_rate * net_e

        if net_e > 0:
            electron_production.append({
                'rate': electron_rate,
                'net_e': net_e,
                'tag': tag,
                'rxn': rxn,
                'k': rate_coeff
            })
        else:
            electron_loss.append({
                'rate': -electron_rate,  # Make positive for loss
                'net_e': net_e,
                'tag': tag,
                'rxn': rxn,
                'k': rate_coeff
            })

# Sort by magnitude
electron_production.sort(key=lambda x: x['rate'], reverse=True)
electron_loss.sort(key=lambda x: x['rate'], reverse=True)

print('=' * 80)
print('ELECTRON PRODUCTION REACTIONS (Top 10)')
print('=' * 80)

total_production = sum(ep['rate'] for ep in electron_production)
for ep in electron_production[:10]:
    rxn = ep['rxn']
    reactants_list = [f"{species[j]}" for j, s in enumerate(rxn.reactants) if s > 0]
    products_list = [f"{species[j]}" for j, s in enumerate(rxn.products) if s > 0]

    rxn_str = ' + '.join(reactants_list) + ' → ' + ' + '.join(products_list)

    print(f"{rxn_str}")
    print(f"  Electron gain: {ep['net_e']:.0f} e per reaction")
    print(f"  Rate: {ep['rate']:.2e} cm⁻³/s ({ep['rate']/total_production*100:.1f}%)")
    print(f"  k = {ep['k']:.2e}, tag = {ep['tag']}")
    print()

print(f"TOTAL ELECTRON PRODUCTION: {total_production:.2e} cm⁻³/s")
print()

print('=' * 80)
print('ELECTRON LOSS REACTIONS (Top 10)')
print('=' * 80)

total_loss = sum(el['rate'] for el in electron_loss)
for el in electron_loss[:10]:
    rxn = el['rxn']
    reactants_list = [f"{species[j]}" for j, s in enumerate(rxn.reactants) if s > 0]
    products_list = [f"{species[j]}" for j, s in enumerate(rxn.products) if s > 0]

    rxn_str = ' + '.join(reactants_list) + ' → ' + ' + '.join(products_list)

    print(f"{rxn_str}")
    print(f"  Electron loss: {-el['net_e']:.0f} e per reaction")
    print(f"  Rate: {el['rate']:.2e} cm⁻³/s ({el['rate']/total_loss*100:.1f}%)")
    print(f"  k = {el['k']:.2e}, tag = {el['tag']}")
    print()

print(f"TOTAL ELECTRON LOSS: {total_loss:.2e} cm⁻³/s")
print()

print('=' * 80)
print('ELECTRON BALANCE SUMMARY')
print('=' * 80)
print(f'Production:  {total_production:.2e} cm⁻³/s')
print(f'Loss:        {total_loss:.2e} cm⁻³/s')
print(f'Net:         {total_production - total_loss:.2e} cm⁻³/s')
print(f'Loss/Prod:   {total_loss/total_production:.2f}×')
print()

if total_loss > total_production:
    print('PROBLEM: Electron loss exceeds production!')
    print(f'  Electrons being destroyed {total_loss/total_production:.2f}× faster than created')
    print()
    print('This explains why Ne is so low and Ni/Ne is so high!')
else:
    print('Electron production exceeds loss, but steady state is still low Ne.')
    print('This suggests electrons are being consumed or neutrals are being over-ionized.')
print()

print('KEY REACTIONS TO INVESTIGATE:')
print('  - Top electron consumers')
print('  - Whether any electron-producing reactions are missing from model')
print('  - Whether reaction rates are physically reasonable at high C2H2')
