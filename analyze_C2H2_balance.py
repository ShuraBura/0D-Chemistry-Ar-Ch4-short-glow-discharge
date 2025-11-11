#!/usr/bin/env python3
"""
Analyze C2H2 production and loss pathways to understand why it can't accumulate
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
print("C2H2 PRODUCTION/LOSS ANALYSIS")
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
idx_C2H2 = species.index('C2H2')
C2H2_final = y_final[idx_C2H2]

print(f"C2H2 at steady state: {C2H2_final:.2e} cm⁻³")
print()

# Now analyze C2H2 production and loss
print("Searching for all C2H2-producing and C2H2-consuming reactions...")
print()

R = params['R']
tags = params['tags']

c2h2_production = []
c2h2_loss = []
c2h2_wall_loss = []

for i, rxn in enumerate(R):
    tag = tags[i]

    # Get net C2H2 change from stoichiometry vectors
    idx_C2H2 = species.index('C2H2')
    net_c2h2 = rxn.products[idx_C2H2] - rxn.reactants[idx_C2H2]

    if net_c2h2 > 0:
        # Production reaction
        c2h2_production.append({
            'index': i,
            'tag': tag,
            'net_c2h2': net_c2h2,
            'rxn': rxn
        })
    elif net_c2h2 < 0:
        # Loss reaction
        if tag.startswith('stick_C2H2') or tag.startswith('loss_C2H2'):
            c2h2_wall_loss.append({
                'index': i,
                'tag': tag,
                'net_c2h2': net_c2h2,
                'rxn': rxn
            })
        else:
            c2h2_loss.append({
                'index': i,
                'tag': tag,
                'net_c2h2': net_c2h2,
                'rxn': rxn
            })

print(f"Found {len(c2h2_production)} C2H2-producing reactions")
print(f"Found {len(c2h2_loss)} C2H2-consuming reactions (chemistry)")
print(f"Found {len(c2h2_wall_loss)} C2H2 wall/volumetric loss terms")
print()

# Calculate rates
print("=" * 80)
print("C2H2 PRODUCTION REACTIONS")
print("=" * 80)

total_production = 0.0
for rxn_info in c2h2_production:
    rxn = rxn_info['rxn']

    # Calculate reaction rate using rxn.rate which is already the rate coefficient
    rate_coeff = rxn.rate

    # Calculate concentrations - multiply densities of reactants
    reactant_density = 1.0
    for j, stoich in enumerate(rxn.reactants):
        if stoich > 0:
            reactant_density *= y_final[j]**stoich

    reaction_rate = rate_coeff * reactant_density  # cm⁻³/s
    production_rate = reaction_rate * rxn_info['net_c2h2']
    total_production += production_rate

    # Format reactants and products
    reactants_list = [f"{int(s) if s > 1 else ''}{species[j]}".strip()
                      for j, s in enumerate(rxn.reactants) if s > 0]
    products_list = [f"{int(s) if s > 1 else ''}{species[j]}".strip()
                     for j, s in enumerate(rxn.products) if s > 0]

    reactants_str = ' + '.join(reactants_list)
    products_str = ' + '.join(products_list)

    print(f"{reactants_str} → {products_str}")
    print(f"  Rate: {production_rate:.2e} cm⁻³/s")
    print(f"  k = {rate_coeff:.2e}, tag = {rxn_info['tag']}")
    print()

print(f"TOTAL C2H2 PRODUCTION: {total_production:.2e} cm⁻³/s")
print()

print("=" * 80)
print("C2H2 LOSS REACTIONS (Chemistry)")
print("=" * 80)

total_loss_chem = 0.0
loss_details = []

for rxn_info in c2h2_loss:
    rxn = rxn_info['rxn']

    # Calculate reaction rate
    rate_coeff = rxn.rate

    # Calculate concentrations
    reactant_density = 1.0
    for j, stoich in enumerate(rxn.reactants):
        if stoich > 0:
            reactant_density *= y_final[j]**stoich

    reaction_rate = rate_coeff * reactant_density
    loss_rate = -reaction_rate * rxn_info['net_c2h2']  # Positive for loss
    total_loss_chem += loss_rate

    loss_details.append({
        'rate': loss_rate,
        'tag': rxn_info['tag'],
        'rxn': rxn,
        'k': rate_coeff
    })

# Sort by rate
loss_details.sort(key=lambda x: x['rate'], reverse=True)

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

print(f"TOTAL C2H2 CHEMICAL LOSS: {total_loss_chem:.2e} cm⁻³/s")
print()

print("=" * 80)
print("C2H2 WALL/VOLUMETRIC LOSS")
print("=" * 80)

total_loss_wall = 0.0
for rxn_info in c2h2_wall_loss:
    rxn = rxn_info['rxn']
    rate_coeff = rxn.rate
    loss_rate = rate_coeff * C2H2_final  # s⁻¹ * cm⁻³ = cm⁻³/s
    total_loss_wall += loss_rate

    print(f"{rxn_info['tag']}")
    print(f"  Rate: {loss_rate:.2e} cm⁻³/s")
    print(f"  k = {rate_coeff:.2e} s⁻¹ (TUNABLE)")
    print()

print(f"TOTAL C2H2 WALL LOSS: {total_loss_wall:.2e} cm⁻³/s")
print()

print("=" * 80)
print("C2H2 BALANCE SUMMARY")
print("=" * 80)
print(f"Production:       {total_production:.2e} cm⁻³/s")
print(f"Chemical loss:    {total_loss_chem:.2e} cm⁻³/s ({total_loss_chem/total_production:.1f}× production)")
print(f"Wall loss:        {total_loss_wall:.2e} cm⁻³/s ({total_loss_wall/total_production:.1f}× production)")
print(f"Total loss:       {total_loss_chem + total_loss_wall:.2e} cm⁻³/s")
print()
print(f"Net rate:         {total_production - total_loss_chem - total_loss_wall:.2e} cm⁻³/s")
print(f"Loss/Production:  {(total_loss_chem + total_loss_wall)/total_production:.2f}×")
print()
print(f"To reach C2H2 = 5e+12 from {C2H2_final:.2e}:")
print(f"  Need {5e12/C2H2_final:.1f}× more C2H2")
print()
