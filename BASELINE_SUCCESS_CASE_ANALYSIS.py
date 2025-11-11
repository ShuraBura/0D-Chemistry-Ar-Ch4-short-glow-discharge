#!/usr/bin/env python3
"""
BASELINE SUCCESS CASE ANALYSIS

This result from optimization_results_comprehensive_1e12/best_f70.3.json
is our DEPARTING POINT going forward.

Result:
  H:    2.01e14 (79.9% of target)  ✓
  CH:   1.01e09 (100.7% of target) ✓✓✓ PERFECT!
  C2:   9.41e08 (16.8% of target)  - Need to improve
  Ni/Ne: 3.12 ✓ Good charge balance

Conditions:
  Te:   1.31 eV
  Ne:   1.22e08 cm⁻³
  E:    250 V/cm
  C2H2: 4.81e09 cm⁻³
  CH3:  6.85e11 cm⁻³
  Ar*:  4.23e08 cm⁻³

This is the ONLY result that hits H and CH simultaneously with good charge balance.
Question: What's limiting C2 to only 16.8%? Can we push it higher without breaking H/CH?
"""

import json
import numpy as np
from define_rates import define_rates
from build_reactions import build_reactions

# Load the baseline success case
with open('optimization_results_comprehensive_1e12/best_f70.3.json', 'r') as f:
    baseline = json.load(f)

print('=' * 80)
print('BASELINE SUCCESS CASE - DETAILED ANALYSIS')
print('=' * 80)
print()

print('TARGET DENSITIES:')
print(f"  H:  {baseline['target_densities']['H']:.2e} ({baseline['target_densities']['H']/2.52e14*100:.1f}%)")
print(f"  CH: {baseline['target_densities']['CH']:.2e} ({baseline['target_densities']['CH']/1.0e9*100:.1f}%)")
print(f"  C2: {baseline['target_densities']['C2']:.2e} ({baseline['target_densities']['C2']/5.6e11*100:.1f}%)")
print()

print('PLASMA CONDITIONS:')
print(f"  Te:        {baseline['Te']:.2f} eV")
print(f"  Ne:        {baseline['Ne']:.2e} cm⁻³")
print(f"  E-field:   {baseline['E_field']:.1f} V/cm")
print(f"  Ni/Ne:     {baseline['Ni_over_Ne']:.2f} ✓ (in range 2-7)")
print()

print('KEY SPECIES:')
species_of_interest = ['C2H2', 'CH3', 'ArStar', 'H', 'CH2', 'C']
for sp in species_of_interest:
    val = baseline['all_densities'].get(sp, 0)
    print(f"  {sp:8s}: {val:.2e} cm⁻³")
print()

print('CRITICAL RATE CONSTANTS:')
# Rebuild params to get rates
params = {
    'P': 500.0, 'Te': baseline['Te'], 'ne': baseline['Ne'], 'E_field': baseline['E_field'],
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
for rate_name, rate_val in baseline.get('rate_values', {}).items():
    if rate_name in k:
        k[rate_name] = rate_val

# Key reactions
key_reactions = {
    'C2H2 → C2 (H + C2H2 → C2)': baseline['rate_values'].get('C2H2_H_C2_H2_H_cm3_7_50', 'N/A'),
    'C2 → CH (H + C2 → CH)': baseline['rate_values'].get('C2_H_CH_C_cm3_7_6', 'N/A'),
    '2CH3 → C2H2': baseline['rate_values'].get('CH3_CH3_C2H2_H2_H2_cm3_7_49', 'N/A'),
    'e + CH4 → CH3': k.get('e_CH4_CH3_H_cm3_1_1', 'N/A'),
    'Ar* + CH4 → CH3': baseline['rate_values'].get('ArStar_CH4_CH3_H_cm3_3_1', 'N/A'),
}

for rxn, rate in key_reactions.items():
    if isinstance(rate, float):
        print(f"  {rxn:30s}: {rate:.2e}")
    else:
        print(f"  {rxn:30s}: {rate}")
print()

# Analyze C2 balance
params['k'] = k
params['R'], params['tags'] = build_reactions(params)
species = params['species']
y0 = np.array([baseline['all_densities'].get(s, 1e3) for s in species])

R = params['R']
tags = params['tags']
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

        if net_c2 > 0:
            c2_production.append({'rate': c2_rate, 'tag': tags[i]})
        else:
            c2_loss.append({'rate': -c2_rate, 'tag': tags[i]})

c2_production.sort(key=lambda x: x['rate'], reverse=True)
c2_loss.sort(key=lambda x: x['rate'], reverse=True)

total_prod = sum(cp['rate'] for cp in c2_production)
total_loss = sum(cl['rate'] for cl in c2_loss)

print('=' * 80)
print('C2 PRODUCTION/LOSS ANALYSIS')
print('=' * 80)
print()
print(f"Total C2 production: {total_prod:.2e} cm⁻³/s")
print(f"Total C2 loss:       {total_loss:.2e} cm⁻³/s")
print(f"Net rate:            {total_prod - total_loss:.2e} cm⁻³/s")
print(f"Steady state C2:     {baseline['target_densities']['C2']:.2e} cm⁻³")
print()

print('Top C2 production reactions:')
for i, cp in enumerate(c2_production[:5]):
    print(f"  {i+1}. {cp['tag']:40s} {cp['rate']:.2e} ({cp['rate']/total_prod*100:.1f}%)")
print()

print('Top C2 loss reactions:')
for i, cl in enumerate(c2_loss[:5]):
    print(f"  {i+1}. {cl['tag']:40s} {cl['rate']:.2e} ({cl['rate']/total_loss*100:.1f}%)")
print()

print('=' * 80)
print('WHY IS C2 ONLY AT 16.8%?')
print('=' * 80)
print()

# Check if it's production-limited or loss-limited
h_c2h2_rxn = [cp for cp in c2_production if 'C2H2_H' in cp['tag']]
if h_c2h2_rxn:
    print(f"H + C2H2 → C2 rate: {h_c2h2_rxn[0]['rate']:.2e} cm⁻³/s")
    print(f"  [C2H2] = {baseline['all_densities']['C2H2']:.2e}")
    print(f"  [H] = {baseline['all_densities']['H']:.2e}")
    print()

h_c2_rxn = [cl for cl in c2_loss if 'C2_H' in cl['tag'] and 'CH' in cl['tag']]
if h_c2_rxn:
    print(f"H + C2 → CH + C rate: {h_c2_rxn[0]['rate']:.2e} cm⁻³/s")
    destruction_fraction = h_c2_rxn[0]['rate'] / total_prod if total_prod > 0 else 0
    print(f"  Destroys {destruction_fraction*100:.1f}% of produced C2")
    print()

print('Possible bottlenecks:')
print(f"  1. C2H2 is only {baseline['all_densities']['C2H2']:.2e} (need ~5e12 per user)")
print(f"     → If C2H2 were 1000× higher, C2 production would be 1000× higher!")
print(f"  2. H + C2 → CH destruction: {destruction_fraction*100:.1f}% of production")
print(f"  3. Wall/volumetric losses")
print()

print('=' * 80)
print('STRATEGY TO IMPROVE C2')
print('=' * 80)
print()
print('Option 1: Increase C2H2 from 4.8e9 to ~5e12 (1000× higher)')
print('  - Need CH3 ~ 1000× higher? No, only ~30× since production ∝ [CH3]²')
print(f'  - Current [CH3] = {baseline["all_densities"]["CH3"]:.2e}')
print(f'  - Need [CH3] ~ {baseline["all_densities"]["CH3"]*30:.2e}')
print('  - BUT this might break the H/CH balance!')
print()
print('Option 2: Reduce C2 destruction (H + C2 → CH)')
print('  - Already tried this in earlier optimizers')
print('  - Limited by literature rates')
print()
print('Option 3: Reduce C2 losses (wall/volumetric)')
print(f'  - Current wall loss: {c2_loss[1]["rate"] if len(c2_loss) > 1 else 0:.2e}')
print(f'  - Might help but probably not enough')
print()
print('RECOMMENDATION:')
print('Start from these exact conditions and carefully tune to push C2H2 higher')
print('without breaking the H/CH balance that we achieved here.')
