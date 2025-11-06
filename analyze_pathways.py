#!/usr/bin/env python3
"""
Analyze reaction pathways for C2 and CH production/consumption
to identify which reactions dominate and what to tune
"""

import json
import numpy as np
from scipy.integrate import solve_ivp
from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun import PlasmaODE

print("="*80)
print("REACTION PATHWAY ANALYSIS FOR C2 AND CH")
print("="*80)
print()

# Load best checkpoint
with open('checkpoint_f3407.json', 'r') as f:
    best_result = json.load(f)

params_dict = best_result['params']
species = params_dict['species']

print(f"Using parameters from checkpoint_f3407.json (f(x) = {best_result['objective']:.2f})")
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

print("Running SHORT simulation for approximate pathway analysis...")
print("(Not full steady-state, but good enough to identify dominant reactions)")
sol = solve_ivp(PlasmaODE(params), (0, 10), y0, method='BDF', rtol=1e-4, atol=1e-6, max_step=5.0)

if not sol.success:
    print(f"ERROR: Simulation failed - {sol.message}")
    exit(1)

print("SUCCESS!")
print()

y_final = sol.y[:, -1]
y_final = np.maximum(y_final, 1e-6)  # Prevent negative

# Now calculate reaction rates at steady state
print("Calculating reaction rates at steady state...")
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

    # Store reaction info
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
    fraction = rxn['net_rate'] / total_c2_production * 100
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
    fraction = rxn['net_rate'] / total_c2_consumption * 100
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
    fraction = rxn['net_rate'] / total_ch_production * 100
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
    fraction = rxn['net_rate'] / total_ch_consumption * 100
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
