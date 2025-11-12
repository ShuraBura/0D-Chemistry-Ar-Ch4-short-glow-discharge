#!/usr/bin/env python3
"""
Find alternative C2H2 production pathways that don't require high Ar*.

Current bottleneck: CH3 comes mainly from Ar* + CH4 → CH3
But if Ar* ~ 1e8 is realistic, we need other sources of:
1. CH3 (for 2CH3 → C2H2)
2. Direct C2H2 production
"""

import json
import numpy as np
from define_rates import define_rates
from build_reactions import build_reactions

# Load the stable result (Ar* ~ 2.58e8, realistic)
with open('optimization_results_fixed_ionization/best_f22.8.json', 'r') as f:
    data = json.load(f)

params = {
    'P': 500.0, 'Te': data['Te'], 'ne': data['Ne'], 'E_field': data['E_field'],
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
for rate_name, rate_val in data.get('rate_values', {}).items():
    if rate_name in k:
        k[rate_name] = rate_val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)

species = params['species']
y0 = np.array([data['all_densities'].get(s, 1e3) for s in species])

R = params['R']
tags = params['tags']

print('=' * 80)
print('ALTERNATIVE C2H2 PRODUCTION PATHWAYS')
print('=' * 80)
print()

# 1. Find all CH3 production reactions
idx_CH3 = species.index('CH3')
ch3_production = []

for i, rxn in enumerate(R):
    net_ch3 = rxn.products[idx_CH3] - rxn.reactants[idx_CH3]
    if net_ch3 > 0:
        rate_coeff = rxn.rate
        reactant_density = 1.0
        for j, stoich in enumerate(rxn.reactants):
            if stoich > 0:
                reactant_density *= y0[j]**stoich

        reaction_rate = rate_coeff * reactant_density
        ch3_rate = reaction_rate * net_ch3

        # Build reaction string
        reactants_list = [f"{stoich:.0f}{species[j]}" if stoich > 1 else species[j]
                         for j, stoich in enumerate(rxn.reactants) if stoich > 0]
        products_list = [f"{stoich:.0f}{species[j]}" if stoich > 1 else species[j]
                        for j, stoich in enumerate(rxn.products) if stoich > 0]
        rxn_str = ' + '.join(reactants_list) + ' → ' + ' + '.join(products_list)

        ch3_production.append({
            'rate': ch3_rate,
            'tag': tags[i],
            'k': rate_coeff,
            'rxn_str': rxn_str
        })

ch3_production.sort(key=lambda x: x['rate'], reverse=True)
total_ch3_prod = sum(cp['rate'] for cp in ch3_production)

print('CH3 PRODUCTION REACTIONS (Top 10)')
print('-' * 80)
for cp in ch3_production[:10]:
    print(f"{cp['rxn_str']}")
    print(f"  Rate: {cp['rate']:.2e} cm⁻³/s ({cp['rate']/total_ch3_prod*100:.1f}%)")
    print(f"  k = {cp['k']:.2e}, tag = {cp['tag']}")
    print()

print(f"TOTAL CH3 PRODUCTION: {total_ch3_prod:.2e} cm⁻³/s")
print()

# 2. Find all C2H2 production reactions
idx_C2H2 = species.index('C2H2')
c2h2_production = []

for i, rxn in enumerate(R):
    net_c2h2 = rxn.products[idx_C2H2] - rxn.reactants[idx_C2H2]
    if net_c2h2 > 0:
        rate_coeff = rxn.rate
        reactant_density = 1.0
        for j, stoich in enumerate(rxn.reactants):
            if stoich > 0:
                reactant_density *= y0[j]**stoich

        reaction_rate = rate_coeff * reactant_density
        c2h2_rate = reaction_rate * net_c2h2

        reactants_list = [f"{stoich:.0f}{species[j]}" if stoich > 1 else species[j]
                         for j, stoich in enumerate(rxn.reactants) if stoich > 0]
        products_list = [f"{stoich:.0f}{species[j]}" if stoich > 1 else species[j]
                        for j, stoich in enumerate(rxn.products) if stoich > 0]
        rxn_str = ' + '.join(reactants_list) + ' → ' + ' + '.join(products_list)

        c2h2_production.append({
            'rate': c2h2_rate,
            'tag': tags[i],
            'k': rate_coeff,
            'rxn_str': rxn_str
        })

c2h2_production.sort(key=lambda x: x['rate'], reverse=True)
total_c2h2_prod = sum(cp['rate'] for cp in c2h2_production)

print('=' * 80)
print('C2H2 PRODUCTION REACTIONS (Top 10)')
print('-' * 80)
for cp in c2h2_production[:10]:
    print(f"{cp['rxn_str']}")
    print(f"  Rate: {cp['rate']:.2e} cm⁻³/s ({cp['rate']/total_c2h2_prod*100:.1f}%)")
    print(f"  k = {cp['k']:.2e}, tag = {cp['tag']}")
    print()

print(f"TOTAL C2H2 PRODUCTION: {total_c2h2_prod:.2e} cm⁻³/s")
print()

# 3. Check densities of key precursors
print('=' * 80)
print('KEY PRECURSOR DENSITIES')
print('-' * 80)
precursors = ['CH3', 'CH2', 'CH', 'C', 'C2H4', 'C2H3', 'C2H5', 'H']
for sp in precursors:
    if sp in species:
        idx = species.index(sp)
        print(f"  {sp:6s}: {y0[idx]:.2e} cm⁻³")
print()

# 4. Identify potential missing reactions or underutilized pathways
print('=' * 80)
print('ANALYSIS: PATHWAYS TO BOOST C2H2')
print('=' * 80)
print()

# Check if electron-impact dissociation of CH4 is in the model
e_ch4_reactions = [tag for tag in tags if 'e_CH4' in tag and 'CH3' in tag]
print('Electron-impact CH4 dissociation reactions:')
if e_ch4_reactions:
    for tag in e_ch4_reactions:
        idx = tags.index(tag)
        print(f"  {tag}: k = {R[idx].rate:.2e}")
else:
    print('  NONE FOUND! This could be a missing pathway.')
print()

# Check C2H4 as C2H2 source
c2h4_to_c2h2 = [tag for tag in tags if 'C2H4' in tag and any(x in tag for x in ['C2H2', 'C2H'])]
print('C2H4 → C2H2 pathways:')
if c2h4_to_c2h2:
    for tag in c2h4_to_c2h2:
        idx = tags.index(tag)
        print(f"  {tag}: k = {R[idx].rate:.2e}")
else:
    print('  Limited C2H4 → C2H2 conversion')
print()

print('RECOMMENDATIONS:')
print('1. If e + CH4 → CH3 + H + e is missing, add it (important CH3 source)')
print('2. Check if C2H4 → C2H2 reactions can be boosted')
print('3. Look for alternative CH3 sources besides Ar* + CH4')
print('4. Consider ion-molecule reactions that produce C2H2 directly')
