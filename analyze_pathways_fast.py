#!/usr/bin/env python3
"""
FAST pathway analysis - uses already-computed densities from optimization
No ODE solving needed - just calculate reaction rates at final state
"""

import json
import numpy as np
from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions

print("="*80)
print("FAST REACTION PATHWAY ANALYSIS FOR C2 AND CH")
print("="*80)
print()

# Load best checkpoint with densities already computed
with open('checkpoint_f3407.json', 'r') as f:
    best_result = json.load(f)

params_dict = best_result['params']
species = params_dict['species']

# Get densities from the optimization (already computed!)
dens_from_opt = best_result['densities']

print(f"Using pre-computed densities from checkpoint_f3407.json (f(x) = {best_result['objective']:.2f})")
print(f"  H:  {dens_from_opt['H']:.3e} cm^-3")
print(f"  CH: {dens_from_opt['CH']:.3e} cm^-3")
print(f"  C2: {dens_from_opt['C2']:.3e} cm^-3")
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

# We need all densities - reconstruct from params
# The optimization only saved H, CH, C2, C2H2
# We need to use reasonable estimates for others or load from all_densities if available

# Try to load from original starting point to get all species
with open('optimization_results_C2H2_boost/best_f259.0_Te1.09.json', 'r') as f:
    start_data = json.load(f)

y_final = np.ones(len(species)) * 1e3

# Try to get densities from all_densities
if 'all_densities' in start_data:
    for sp in species:
        if sp in start_data['all_densities']:
            y_final[species.index(sp)] = start_data['all_densities'][sp]

# Override with optimized target species
y_final[species.index('H')] = dens_from_opt['H']
y_final[species.index('CH')] = dens_from_opt['CH']
y_final[species.index('C2')] = dens_from_opt['C2']
y_final[species.index('C2H2')] = dens_from_opt['C2H2']
y_final[species.index('e')] = params_dict['ne']
y_final[species.index('Ar')] = 0.85 * 9.66e15
y_final[species.index('CH4')] = 0.15 * 9.66e15

y_final = np.maximum(y_final, 1e-6)

print("Calculating reaction rates at these densities...")
print()

C2_idx = species.index('C2')
CH_idx = species.index('CH')
H_idx = species.index('H')
e_idx = species.index('e')

# Calculate all reaction rates
reaction_rates = []

for rxn_idx, reaction in enumerate(params['R']):
    rate_constant = params['k'][params['tags'][rxn_idx]]
    tag = params['tags'][rxn_idx]

    # Calculate rate
    if 'drift' in tag:
        rate = rate_constant * y_final[e_idx]
    elif 'loss' in tag or 'stick' in tag:
        react_species = np.where(reaction.reactants > 0)[0]
        if len(react_species) > 0:
            rate = rate_constant * y_final[react_species[0]]
        else:
            rate = 0.0
    else:
        rate = rate_constant
        react_species = np.where(reaction.reactants > 0)[0]
        for sp_idx in react_species:
            rate *= y_final[sp_idx] ** reaction.reactants[sp_idx]

    # Build reaction string
    reactants_str = []
    for i, sp in enumerate(species):
        if reaction.reactants[i] > 0:
            if reaction.reactants[i] == 1:
                reactants_str.append(sp)
            else:
                reactants_str.append(f"{int(reaction.reactants[i])}{sp}")

    products_str = []
    for i, sp in enumerate(species):
        if reaction.products[i] > 0:
            if reaction.products[i] == 1:
                products_str.append(sp)
            else:
                products_str.append(f"{int(reaction.products[i])}{sp}")

    reaction_str = " + ".join(reactants_str) + " → " + " + ".join(products_str)

    reaction_rates.append({
        'idx': rxn_idx,
        'tag': tag,
        'reaction': reaction_str,
        'rate': rate,
        'rate_constant': rate_constant,
        'C2_produced': reaction.products[C2_idx],
        'C2_consumed': reaction.reactants[C2_idx],
        'CH_produced': reaction.products[CH_idx],
        'CH_consumed': reaction.reactants[CH_idx],
        'H_produced': reaction.products[H_idx],
        'H_consumed': reaction.reactants[H_idx],
    })

# Analyze C2 pathways
print("="*80)
print("C2 PRODUCTION REACTIONS (sorted by net C2 production rate)")
print("="*80)
print()

c2_production = []
for rxn in reaction_rates:
    net_c2 = rxn['C2_produced'] - rxn['C2_consumed']
    if net_c2 > 0:
        net_rate = net_c2 * rxn['rate']
        c2_production.append({
            'reaction': rxn['reaction'],
            'tag': rxn['tag'],
            'net_c2': net_c2,
            'rate': rxn['rate'],
            'net_rate': net_rate,
            'rate_constant': rxn['rate_constant']
        })

c2_production.sort(key=lambda x: x['net_rate'], reverse=True)

total_c2_production = sum(r['net_rate'] for r in c2_production)

print(f"Total C2 production: {total_c2_production:.3e} cm^-3/s")
print()
print("Top 15 C2 production reactions:")
print()

for i, rxn in enumerate(c2_production[:15]):
    if total_c2_production > 0:
        fraction = rxn['net_rate'] / total_c2_production * 100
    else:
        fraction = 0
    print(f"{i+1:2d}. [{fraction:5.1f}%] {rxn['net_rate']:10.3e} cm^-3/s")
    print(f"    {rxn['reaction']}")
    print(f"    tag: {rxn['tag']}")
    print(f"    k = {rxn['rate_constant']:.3e}")
    print()

print()
print("="*80)
print("C2 CONSUMPTION REACTIONS (sorted by net C2 consumption rate)")
print("="*80)
print()

c2_consumption = []
for rxn in reaction_rates:
    net_c2 = rxn['C2_consumed'] - rxn['C2_produced']
    if net_c2 > 0:
        net_rate = net_c2 * rxn['rate']
        c2_consumption.append({
            'reaction': rxn['reaction'],
            'tag': rxn['tag'],
            'net_c2': net_c2,
            'rate': rxn['rate'],
            'net_rate': net_rate,
            'rate_constant': rxn['rate_constant']
        })

c2_consumption.sort(key=lambda x: x['net_rate'], reverse=True)

total_c2_consumption = sum(r['net_rate'] for r in c2_consumption)

print(f"Total C2 consumption: {total_c2_consumption:.3e} cm^-3/s")
print()
print("Top 15 C2 consumption reactions:")
print()

for i, rxn in enumerate(c2_consumption[:15]):
    if total_c2_consumption > 0:
        fraction = rxn['net_rate'] / total_c2_consumption * 100
    else:
        fraction = 0
    print(f"{i+1:2d}. [{fraction:5.1f}%] {rxn['net_rate']:10.3e} cm^-3/s")
    print(f"    {rxn['reaction']}")
    print(f"    tag: {rxn['tag']}")
    print(f"    k = {rxn['rate_constant']:.3e}")
    print()

print()
print(f"NET C2 RATE: {total_c2_production - total_c2_consumption:.3e} cm^-3/s")
print()

# Analyze CH pathways
print("="*80)
print("CH PRODUCTION REACTIONS (sorted by net CH production rate)")
print("="*80)
print()

ch_production = []
for rxn in reaction_rates:
    net_ch = rxn['CH_produced'] - rxn['CH_consumed']
    if net_ch > 0:
        net_rate = net_ch * rxn['rate']
        ch_production.append({
            'reaction': rxn['reaction'],
            'tag': rxn['tag'],
            'net_ch': net_ch,
            'rate': rxn['rate'],
            'net_rate': net_rate,
            'rate_constant': rxn['rate_constant']
        })

ch_production.sort(key=lambda x: x['net_rate'], reverse=True)

total_ch_production = sum(r['net_rate'] for r in ch_production)

print(f"Total CH production: {total_ch_production:.3e} cm^-3/s")
print()
print("Top 15 CH production reactions:")
print()

for i, rxn in enumerate(ch_production[:15]):
    if total_ch_production > 0:
        fraction = rxn['net_rate'] / total_ch_production * 100
    else:
        fraction = 0
    print(f"{i+1:2d}. [{fraction:5.1f}%] {rxn['net_rate']:10.3e} cm^-3/s")
    print(f"    {rxn['reaction']}")
    print(f"    tag: {rxn['tag']}")
    print(f"    k = {rxn['rate_constant']:.3e}")
    print()

print()
print("="*80)
print("CH CONSUMPTION REACTIONS (sorted by net CH consumption rate)")
print("="*80)
print()

ch_consumption = []
for rxn in reaction_rates:
    net_ch = rxn['CH_consumed'] - rxn['CH_produced']
    if net_ch > 0:
        net_rate = net_ch * rxn['rate']
        ch_consumption.append({
            'reaction': rxn['reaction'],
            'tag': rxn['tag'],
            'net_ch': net_ch,
            'rate': rxn['rate'],
            'net_rate': net_rate,
            'rate_constant': rxn['rate_constant']
        })

ch_consumption.sort(key=lambda x: x['net_rate'], reverse=True)

total_ch_consumption = sum(r['net_rate'] for r in ch_consumption)

print(f"Total CH consumption: {total_ch_consumption:.3e} cm^-3/s")
print()
print("Top 15 CH consumption reactions:")
print()

for i, rxn in enumerate(ch_consumption[:15]):
    if total_ch_consumption > 0:
        fraction = rxn['net_rate'] / total_ch_consumption * 100
    else:
        fraction = 0
    print(f"{i+1:2d}. [{fraction:5.1f}%] {rxn['net_rate']:10.3e} cm^-3/s")
    print(f"    {rxn['reaction']}")
    print(f"    tag: {rxn['tag']}")
    print(f"    k = {rxn['rate_constant']:.3e}")
    print()

print()
print(f"NET CH RATE: {total_ch_production - total_ch_consumption:.3e} cm^-3/s")
print()

# Save results
output = {
    'C2_production': c2_production[:20],
    'C2_consumption': c2_consumption[:20],
    'CH_production': ch_production[:20],
    'CH_consumption': ch_consumption[:20],
    'totals': {
        'C2_production': total_c2_production,
        'C2_consumption': total_c2_consumption,
        'CH_production': total_ch_production,
        'CH_consumption': total_ch_consumption
    }
}

with open('pathway_analysis_C2_CH.json', 'w') as f:
    json.dump(output, f, indent=2)

print("="*80)
print("Results saved to: pathway_analysis_C2_CH.json")
print("="*80)
