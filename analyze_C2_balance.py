#!/usr/bin/env python3
"""
Analyze C2 production and loss - the H+C2H2→C2 reaction produces 4.09e+14 cm⁻³/s
but C2 is only 1.04e+10! Where is it going?
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

# Load C2H2 optimizer result
with open('optimization_results_C2H2_tunable_loss/best_f248.8.json', 'r') as f:
    best = json.load(f)

print("=" * 80)
print("C2 PRODUCTION/LOSS ANALYSIS")
print("=" * 80)
print()

# Setup params
params = {
    'P': 500.0,
    'Te': best['Te'],
    'ne': best['Ne'],
    'E_field': best['E_field'],
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

n_total = pressure_to_density(500.0)
k = define_rates(params)

# Apply tuned rates
for rate_name, rate_val in best.get('rate_values', {}).items():
    if rate_name in k:
        k[rate_name] = rate_val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Set initial conditions from best result
species = params['species']
y0 = np.array([best['all_densities'].get(s, 1e3) for s in species])

# Run to steady state
print("Running to t=2000s steady state...")
ode_func = PlasmaODE_Optimized(params)

sol = solve_ivp(
    ode_func,
    (0, 2000),
    y0,
    method='BDF',
    rtol=1e-7,
    atol=1e-9,
    max_step=1.0
)

y_final = sol.y[:, -1]
idx_C2 = species.index('C2')
C2_final = y_final[idx_C2]

print(f"C2 at steady state: {C2_final:.2e} cm⁻³")
print(f"C2 target: 5.6e+11 cm⁻³ ({C2_final/5.6e11*100:.1f}%)")
print()

# Analyze C2 production and loss
print("Searching for all C2-producing and C2-consuming reactions...")
print()

R = params['R']
tags = params['tags']

c2_production = []
c2_loss = []
c2_wall_loss = []

for i, rxn in enumerate(R):
    tag = tags[i]

    # Get net C2 change
    idx_C2 = species.index('C2')
    net_c2 = rxn.products[idx_C2] - rxn.reactants[idx_C2]

    if net_c2 > 0:
        c2_production.append({
            'index': i,
            'tag': tag,
            'net_c2': net_c2,
            'rxn': rxn
        })
    elif net_c2 < 0:
        if tag.startswith('stick_C2') or tag.startswith('loss_C2'):
            c2_wall_loss.append({
                'index': i,
                'tag': tag,
                'net_c2': net_c2,
                'rxn': rxn
            })
        else:
            c2_loss.append({
                'index': i,
                'tag': tag,
                'net_c2': net_c2,
                'rxn': rxn
            })

print(f"Found {len(c2_production)} C2-producing reactions")
print(f"Found {len(c2_loss)} C2-consuming reactions (chemistry)")
print(f"Found {len(c2_wall_loss)} C2 wall/volumetric loss terms")
print()

# Calculate production rates
print("=" * 80)
print("C2 PRODUCTION REACTIONS (Top 10)")
print("=" * 80)

production_details = []
for rxn_info in c2_production:
    rxn = rxn_info['rxn']
    rate_coeff = rxn.rate

    # Calculate concentrations
    reactant_density = 1.0
    for j, stoich in enumerate(rxn.reactants):
        if stoich > 0:
            reactant_density *= y_final[j]**stoich

    reaction_rate = rate_coeff * reactant_density
    production_rate = reaction_rate * rxn_info['net_c2']

    production_details.append({
        'rate': production_rate,
        'tag': rxn_info['tag'],
        'rxn': rxn,
        'k': rate_coeff
    })

# Sort by rate
production_details.sort(key=lambda x: x['rate'], reverse=True)
total_production = sum(pd['rate'] for pd in production_details)

for pd in production_details[:10]:  # Top 10
    rxn = pd['rxn']
    reactants_list = [f"{int(s) if s > 1 else ''}{species[j]}".strip()
                      for j, s in enumerate(rxn.reactants) if s > 0]
    products_list = [f"{int(s) if s > 1 else ''}{species[j]}".strip()
                     for j, s in enumerate(rxn.products) if s > 0]

    reactants_str = ' + '.join(reactants_list)
    products_str = ' + '.join(products_list)

    print(f"{reactants_str} → {products_str}")
    print(f"  Rate: {pd['rate']:.2e} cm⁻³/s ({pd['rate']/total_production*100:.1f}%)")
    print(f"  k = {pd['k']:.2e}, tag = {pd['tag']}")
    print()

print(f"TOTAL C2 PRODUCTION: {total_production:.2e} cm⁻³/s")
print()

# Calculate loss rates
print("=" * 80)
print("C2 LOSS REACTIONS (Chemistry) - Top 10")
print("=" * 80)

loss_details = []
for rxn_info in c2_loss:
    rxn = rxn_info['rxn']
    rate_coeff = rxn.rate

    # Calculate concentrations
    reactant_density = 1.0
    for j, stoich in enumerate(rxn.reactants):
        if stoich > 0:
            reactant_density *= y_final[j]**stoich

    reaction_rate = rate_coeff * reactant_density
    loss_rate = -reaction_rate * rxn_info['net_c2']  # Positive for loss

    loss_details.append({
        'rate': loss_rate,
        'tag': rxn_info['tag'],
        'rxn': rxn,
        'k': rate_coeff
    })

# Sort by rate
loss_details.sort(key=lambda x: x['rate'], reverse=True)
total_loss_chem = sum(ld['rate'] for ld in loss_details)

for ld in loss_details[:10]:  # Top 10
    rxn = ld['rxn']
    reactants_list = [f"{int(s) if s > 1 else ''}{species[j]}".strip()
                      for j, s in enumerate(rxn.reactants) if s > 0]
    products_list = [f"{int(s) if s > 1 else ''}{species[j]}".strip()
                     for j, s in enumerate(rxn.products) if s > 0]

    reactants_str = ' + '.join(reactants_list)
    products_str = ' + '.join(products_list)

    print(f"{reactants_str} → {products_str}")
    print(f"  Rate: {ld['rate']:.2e} cm⁻³/s ({ld['rate']/total_loss_chem*100:.1f}%)")
    print(f"  k = {ld['k']:.2e}, tag = {ld['tag']}")
    print()

print(f"TOTAL C2 CHEMICAL LOSS: {total_loss_chem:.2e} cm⁻³/s")
print()

# Calculate wall loss
print("=" * 80)
print("C2 WALL/VOLUMETRIC LOSS")
print("=" * 80)

total_loss_wall = 0.0
for rxn_info in c2_wall_loss:
    rxn = rxn_info['rxn']
    rate_coeff = rxn.rate
    loss_rate = rate_coeff * C2_final  # s⁻¹ * cm⁻³ = cm⁻³/s
    total_loss_wall += loss_rate

    print(f"{rxn_info['tag']}")
    print(f"  Rate: {loss_rate:.2e} cm⁻³/s")
    print(f"  k = {rate_coeff:.2e} s⁻¹")
    print()

print(f"TOTAL C2 WALL LOSS: {total_loss_wall:.2e} cm⁻³/s")
print()

# Summary
print("=" * 80)
print("C2 BALANCE SUMMARY")
print("=" * 80)
print(f"Production:       {total_production:.2e} cm⁻³/s")
print(f"Chemical loss:    {total_loss_chem:.2e} cm⁻³/s ({total_loss_chem/total_production:.1f}× production)")
print(f"Wall loss:        {total_loss_wall:.2e} cm⁻³/s ({total_loss_wall/total_production:.1f}× production)")
print(f"Total loss:       {total_loss_chem + total_loss_wall:.2e} cm⁻³/s")
print()
print(f"Net rate:         {total_production - total_loss_chem - total_loss_wall:.2e} cm⁻³/s")
print(f"Loss/Production:  {(total_loss_chem + total_loss_wall)/total_production:.2f}×")
print()
print(f"KEY INSIGHT:")
print(f"  C2 is at {C2_final:.2e} (target: 5.6e+11)")
print(f"  Need {5.6e11/C2_final:.1f}× more C2 to reach target")
print(f"  To achieve this, need to reduce C2 LOSS or increase C2 production")
print()
