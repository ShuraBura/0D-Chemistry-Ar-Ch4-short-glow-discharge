#!/usr/bin/env python3
"""
Re-analyze C2 production/loss with the correct perspective:
C2H2 is a MEANS to achieve C2, not the goal!

The pathway: CH3 → C2H2 → C2
- H + C2H2 → C2 + H2 is GOOD (produces C2)
- H + C2 → CH + C is BAD (destroys C2)

Question: Why can't we accumulate C2 when we have the production pathway?
"""

import json
import numpy as np
from define_rates import define_rates
from build_reactions import build_reactions

# Load stable result
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
print('C2 BOTTLENECK ANALYSIS')
print('=' * 80)
print()

print('Current State:')
print(f'  [C2]:   {y0[species.index("C2")]:.2e} cm⁻³ (target: 5.6e11, achieved: {y0[species.index("C2")]/5.6e11*100:.1f}%)')
print(f'  [C2H2]: {y0[species.index("C2H2")]:.2e} cm⁻³')
print(f'  [H]:    {y0[species.index("H")]:.2e} cm⁻³')
print(f'  [CH3]:  {y0[species.index("CH3")]:.2e} cm⁻³')
print()

# Analyze C2 balance
idx_C2 = species.index('C2')
c2_production = []
c2_loss = []

for i, rxn in enumerate(R):
    net_c2 = rxn.products[idx_C2] - rxn.reactants[idx_C2]
    if net_c2 != 0:
        rate_coeff = rxn.rate
        reactant_density = 1.0
        for j, stoich in enumerate(rxn.reactants):
            if stoich > 0:
                reactant_density *= y0[j]**stoich

        reaction_rate = rate_coeff * reactant_density
        c2_rate = reaction_rate * net_c2

        reactants_list = []
        for j, stoich in enumerate(rxn.reactants):
            if stoich > 0:
                reactants_list.append(f"{int(stoich)}{species[j]}" if stoich > 1 else species[j])

        products_list = []
        for j, stoich in enumerate(rxn.products):
            if stoich > 0:
                products_list.append(f"{int(stoich)}{species[j]}" if stoich > 1 else species[j])

        rxn_str = ' + '.join(reactants_list) + ' → ' + ' + '.join(products_list)

        if net_c2 > 0:
            c2_production.append({
                'rate': c2_rate,
                'tag': tags[i],
                'k': rate_coeff,
                'rxn_str': rxn_str
            })
        else:
            c2_loss.append({
                'rate': -c2_rate,
                'tag': tags[i],
                'k': rate_coeff,
                'rxn_str': rxn_str
            })

c2_production.sort(key=lambda x: x['rate'], reverse=True)
c2_loss.sort(key=lambda x: x['rate'], reverse=True)

total_prod = sum(cp['rate'] for cp in c2_production)
total_loss = sum(cl['rate'] for cl in c2_loss)

print('=' * 80)
print('C2 PRODUCTION REACTIONS (Top 5)')
print('=' * 80)
for i, cp in enumerate(c2_production[:5]):
    print(f"{i+1}. {cp['rxn_str']}")
    print(f"   Rate: {cp['rate']:.2e} cm⁻³/s ({cp['rate']/total_prod*100:.1f}%)")
    print(f"   k = {cp['k']:.2e}, tag = {cp['tag']}")
    print()

print(f"TOTAL C2 PRODUCTION: {total_prod:.2e} cm⁻³/s")
print()

print('=' * 80)
print('C2 LOSS REACTIONS (Top 5)')
print('=' * 80)
for i, cl in enumerate(c2_loss[:5]):
    print(f"{i+1}. {cl['rxn_str']}")
    print(f"   Rate: {cl['rate']:.2e} cm⁻³/s ({cl['rate']/total_loss*100:.1f}%)")
    print(f"   k = {cl['k']:.2e}, tag = {cl['tag']}")
    print()

print(f"TOTAL C2 LOSS: {total_loss:.2e} cm⁻³/s")
print()

print('=' * 80)
print('THE PATHWAY: C2H2 → C2')
print('=' * 80)
print()

# Find H + C2H2 → C2 reaction
h_c2h2_to_c2 = [cp for cp in c2_production if 'C2H2_H' in cp['tag'] or 'H_C2H2' in cp['tag']]
if h_c2h2_to_c2:
    print('H + C2H2 → C2 + H2 + H (the desired pathway!):')
    for rxn in h_c2h2_to_c2:
        print(f"  Rate: {rxn['rate']:.2e} cm⁻³/s ({rxn['rate']/total_prod*100:.1f}% of total C2 production)")
        print(f"  k = {rxn['k']:.2e}")
    print()

# Find H + C2 → CH destruction
h_c2_to_ch = [cl for cl in c2_loss if 'C2_H' in cl['tag'] and 'CH' in cl['tag']]
if h_c2_to_ch:
    print('H + C2 → CH + C (destroys the C2 we just made!):')
    for rxn in h_c2_to_ch:
        print(f"  Rate: {rxn['rate']:.2e} cm⁻³/s ({rxn['rate']/total_loss*100:.1f}% of total C2 loss)")
        print(f"  k = {rxn['k']:.2e}")
    print()

print('=' * 80)
print('THE VICIOUS CYCLE')
print('=' * 80)
print()

if h_c2h2_to_c2 and h_c2_to_ch:
    c2_from_c2h2 = sum(r['rate'] for r in h_c2h2_to_c2)
    c2_to_ch = sum(r['rate'] for r in h_c2_to_ch)

    print(f"Step 1: H + C2H2 → C2      (+{c2_from_c2h2:.2e} cm⁻³/s)")
    print(f"Step 2: H + C2 → CH + C    (-{c2_to_ch:.2e} cm⁻³/s)")
    print()
    print(f"Net C2 gain: {c2_from_c2h2 - c2_to_ch:.2e} cm⁻³/s")
    print(f"Destruction efficiency: {c2_to_ch/c2_from_c2h2*100:.1f}% of produced C2 is destroyed!")
    print()

print('=' * 80)
print('WHY C2 STAYS LOW')
print('=' * 80)
print()

print('The Problem:')
print('  1. H atoms convert C2H2 → C2 (GOOD!)')
print('  2. Same H atoms convert C2 → CH (BAD!)')
print('  3. With [H] = 2e14, BOTH reactions are very fast')
print('  4. Rate(C2 → CH) / Rate(C2H2 → C2) = {:.1f}%'.format(
    sum(r['rate'] for r in h_c2_to_ch) / sum(r['rate'] for r in h_c2h2_to_c2) * 100))
print('  5. Most of the C2 produced is immediately destroyed!')
print()

print('Why can\'t we just reduce H to stop C2 destruction?')
print('  - H is needed to convert C2H2 → C2 in the first place!')
print('  - AND H = 2.52e14 is a target (experimental measurement)')
print()

print('User\'s insight: "C2 can outnumber CH by 2 orders if C2H2 ≥ 5e12"')
print(f'  Current [C2H2]: {y0[species.index("C2H2")]:.2e}')
print(f'  Need [C2H2] ~ 5e12 (1400× higher!)')
print()
print('If we had C2H2 = 5e12:')
print('  - H + C2H2 → C2 rate would be 1400× higher')
print('  - This would saturate C2 production')
print('  - Even with high H + C2 → CH destruction, net C2 would be much higher')
print()

print('=' * 80)
print('CONCLUSION')
print('=' * 80)
print()
print('The bottleneck is NOT the C2H2 → C2 conversion (H does this well).')
print('The bottleneck is getting C2H2 high enough in the first place!')
print()
print('Current [C2H2] = 3.54e9 is 1400× too low to feed enough C2 production')
print('to overcome the C2 → CH destruction.')
print()
print('Back to the original question: Why can\'t we accumulate C2H2?')
print('  - We showed it\'s not CH3 production (boosted alternatives)')
print('  - We showed it\'s not Ar* (creates charge imbalance)')
print('  - H "destroys" C2H2 via H + C2H2 → C2, but this is DESIRED!')
print()
print('The question is: Can we get C2H2 ~ 5e12 with realistic physics?')
