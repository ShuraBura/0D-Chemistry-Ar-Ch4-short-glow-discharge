#!/usr/bin/env python3
"""
Detailed analysis of C2H2 loss mechanisms.

Question: Is C2H2 being produced but then rapidly destroyed?
Compare production vs loss rates and identify dominant loss pathways.
"""

import json
import numpy as np
from define_rates import define_rates
from build_reactions import build_reactions

def analyze_c2h2_balance(filename, label):
    with open(filename, 'r') as f:
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
    idx_C2H2 = species.index('C2H2')

    c2h2_production = []
    c2h2_loss = []

    for i, rxn in enumerate(R):
        net_c2h2 = rxn.products[idx_C2H2] - rxn.reactants[idx_C2H2]
        if net_c2h2 != 0:
            rate_coeff = rxn.rate
            reactant_density = 1.0
            for j, stoich in enumerate(rxn.reactants):
                if stoich > 0:
                    reactant_density *= y0[j]**stoich

            reaction_rate = rate_coeff * reactant_density
            c2h2_rate = reaction_rate * net_c2h2

            # Build reaction string
            reactants_list = []
            for j, stoich in enumerate(rxn.reactants):
                if stoich > 0:
                    if stoich > 1:
                        reactants_list.append(f"{int(stoich)}{species[j]}")
                    else:
                        reactants_list.append(species[j])

            products_list = []
            for j, stoich in enumerate(rxn.products):
                if stoich > 0:
                    if stoich > 1:
                        products_list.append(f"{int(stoich)}{species[j]}")
                    else:
                        products_list.append(species[j])

            rxn_str = ' + '.join(reactants_list) + ' → ' + ' + '.join(products_list)

            if net_c2h2 > 0:
                c2h2_production.append({
                    'rate': c2h2_rate,
                    'tag': tags[i],
                    'k': rate_coeff,
                    'rxn_str': rxn_str,
                    'reactant_density': reactant_density
                })
            else:
                c2h2_loss.append({
                    'rate': -c2h2_rate,
                    'tag': tags[i],
                    'k': rate_coeff,
                    'rxn_str': rxn_str,
                    'reactant_density': reactant_density
                })

    c2h2_production.sort(key=lambda x: x['rate'], reverse=True)
    c2h2_loss.sort(key=lambda x: x['rate'], reverse=True)

    total_prod = sum(cp['rate'] for cp in c2h2_production)
    total_loss = sum(cl['rate'] for cl in c2h2_loss)

    return {
        'label': label,
        'C2H2': data['all_densities']['C2H2'],
        'total_prod': total_prod,
        'total_loss': total_loss,
        'production': c2h2_production,
        'loss': c2h2_loss,
        'key_species': {
            'H': y0[species.index('H')],
            'CH3': y0[species.index('CH3')],
            'Ar*': y0[species.index('ArStar')],
            'Ni/Ne': data['Ni_over_Ne']
        }
    }

print('=' * 80)
print('C2H2 LOSS MECHANISM ANALYSIS')
print('=' * 80)
print()

# Analyze the stable, realistic result
stable = analyze_c2h2_balance('optimization_results_fixed_ionization/best_f22.8.json',
                               'STABLE (Fixed Ionization)')

print(f'{stable["label"]}:')
print('-' * 80)
print(f'  C2H2:     {stable["C2H2"]:.2e} cm⁻³')
print(f'  Ar*:      {stable["key_species"]["Ar*"]:.2e} cm⁻³')
print(f'  CH3:      {stable["key_species"]["CH3"]:.2e} cm⁻³')
print(f'  H:        {stable["key_species"]["H"]:.2e} cm⁻³')
print(f'  Ni/Ne:    {stable["key_species"]["Ni/Ne"]:.2f}')
print()
print(f'  C2H2 production: {stable["total_prod"]:.2e} cm⁻³/s')
print(f'  C2H2 loss:       {stable["total_loss"]:.2e} cm⁻³/s')
print(f'  Net:             {stable["total_prod"] - stable["total_loss"]:.2e} cm⁻³/s')
print(f'  Production/Loss: {stable["total_prod"]/stable["total_loss"]:.3f}')
print()

print('=' * 80)
print('TOP C2H2 PRODUCTION REACTIONS')
print('=' * 80)
for i, cp in enumerate(stable['production'][:10]):
    print(f"{i+1}. {cp['rxn_str']}")
    print(f"   Rate: {cp['rate']:.2e} cm⁻³/s ({cp['rate']/stable['total_prod']*100:.1f}%)")
    print(f"   k = {cp['k']:.2e}, tag = {cp['tag']}")
    print()

print('=' * 80)
print('TOP C2H2 LOSS REACTIONS')
print('=' * 80)
for i, cl in enumerate(stable['loss'][:10]):
    print(f"{i+1}. {cl['rxn_str']}")
    print(f"   Rate: {cl['rate']:.2e} cm⁻³/s ({cl['rate']/stable['total_loss']*100:.1f}%)")
    print(f"   k = {cl['k']:.2e}, tag = {cl['tag']}")

    # Identify loss type
    if 'stick' in cl['tag']:
        loss_type = 'WALL LOSS'
    elif 'loss' in cl['tag']:
        loss_type = 'VOLUMETRIC LOSS'
    elif 'e_C2H2' in cl['tag']:
        loss_type = 'ELECTRON IMPACT'
    elif 'H_C2H2' in cl['tag'] or 'C2H2_H' in cl['tag']:
        loss_type = 'H ATOM REACTIONS'
    else:
        loss_type = 'OTHER'

    print(f"   Type: {loss_type}")
    print()

print('=' * 80)
print('C2H2 LOSS BREAKDOWN BY TYPE')
print('=' * 80)

wall_loss = sum(cl['rate'] for cl in stable['loss'] if 'stick' in cl['tag'])
vol_loss = sum(cl['rate'] for cl in stable['loss'] if 'loss' in cl['tag'])
electron_loss = sum(cl['rate'] for cl in stable['loss'] if cl['tag'].startswith('e_C2H2'))
h_reactions = sum(cl['rate'] for cl in stable['loss']
                  if 'H_C2H2' in cl['tag'] or 'C2H2_H' in cl['tag'])
other_loss = stable['total_loss'] - (wall_loss + vol_loss + electron_loss + h_reactions)

print(f"Wall sticking:      {wall_loss:.2e} cm⁻³/s ({wall_loss/stable['total_loss']*100:.1f}%)")
print(f"Volumetric loss:    {vol_loss:.2e} cm⁻³/s ({vol_loss/stable['total_loss']*100:.1f}%)")
print(f"Electron impact:    {electron_loss:.2e} cm⁻³/s ({electron_loss/stable['total_loss']*100:.1f}%)")
print(f"H atom reactions:   {h_reactions:.2e} cm⁻³/s ({h_reactions/stable['total_loss']*100:.1f}%)")
print(f"Other:              {other_loss:.2e} cm⁻³/s ({other_loss/stable['total_loss']*100:.1f}%)")
print()

print('=' * 80)
print('ANALYSIS: WHAT LIMITS C2H2?')
print('=' * 80)
print()

# Calculate C2H2 lifetime
if stable['total_loss'] > 0:
    lifetime = stable['C2H2'] / stable['total_loss']
    print(f"C2H2 lifetime: {lifetime*1000:.2f} ms")
    print()

# Check if production/loss are balanced
if abs(stable['total_prod'] - stable['total_loss']) / stable['total_prod'] < 0.01:
    print("✓ Production and loss are balanced (steady state)")
else:
    print(f"✗ NOT at steady state: prod/loss = {stable['total_prod']/stable['total_loss']:.2f}")
print()

# Identify bottleneck
print("BOTTLENECK ANALYSIS:")
print()

# Check if we could reduce losses
print("Could we reduce C2H2 losses?")
print(f"  - Wall sticking coefficient: {stable['loss'][0]['k']:.2e} s⁻¹")
print(f"    (Is this realistic? Could it be lower?)")
print()
print(f"  - Volumetric loss rate: {stable['loss'][1]['k']:.2e} s⁻¹ if second loss is volumetric")
print(f"    (What physical process is this?)")
print()

# Check H + C2H2 reaction
h_c2h2_rxns = [cl for cl in stable['loss'] if 'H_C2H2' in cl['tag'] or 'C2H2_H' in cl['tag']]
if h_c2h2_rxns:
    print(f"  - H + C2H2 reactions: {sum(r['rate'] for r in h_c2h2_rxns):.2e} cm⁻³/s")
    print(f"    With [H] = {stable['key_species']['H']:.2e}, this is significant!")
    print(f"    Could these rates be overestimated?")
    print()

# Check if production could be increased
print("Could we increase C2H2 production?")
print(f"  - Main producer: 2CH3 → C2H2 (98% of production)")
print(f"  - Current [CH3]: {stable['key_species']['CH3']:.2e} cm⁻³")
print(f"  - Production ∝ [CH3]²")
print(f"  - To get 10× more C2H2 production, need [CH3] ~ {stable['key_species']['CH3']*np.sqrt(10):.2e}")
print(f"    (3.2× higher CH3)")
print()

print("RECOMMENDATION:")
print("Examine the loss rate constants in define_rates.py:")
print("  1. Are wall sticking coefficients realistic?")
print("  2. What is the 'volumetric loss' term physically?")
print("  3. Are H + C2H2 reaction rates from literature?")
print("  4. Could any of these be overestimated?")
