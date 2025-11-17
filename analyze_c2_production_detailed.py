#!/usr/bin/env python3
"""
Detailed analysis of C2 production pathways to find the 1263× gap

Current: C2 = 4.43×10⁸ cm⁻³
Target:  C2 = 5.6×10¹¹ cm⁻³
Gap:     1263× too low

Goal: Identify which production pathway(s) can account for this gap
"""

import numpy as np
from scipy.integrate import solve_ivp
from build_reactions import build_reactions
from define_rates import define_rates
from odefun_optimized import PlasmaODE_Optimized

# Set conditions
PRESSURE_MTORR = 400
TGAS_K = 570
NE = 2.3e9
TE_EV = 1.3

def pressure_to_density(pressure_mTorr, T_K):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

n_total = pressure_to_density(PRESSURE_MTORR, TGAS_K)

params = {
    'P': PRESSURE_MTORR,
    'Tgas': TGAS_K,
    'Te': TE_EV,
    'ne': NE,
    'E_field': 150,
    'L_discharge': 0.45,
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
    },
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4',
                'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3',
                'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3',
                'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'C3H5', 'HMinus', 'C2H5Plus',
                'CH2Plus', 'C4H', 'H2Plus', 'C2H6Plus', 'C2H2Star', 'C3H6'],
    'stick_coeffs': {
        'ArPlus': 1.0, 'CH4Plus': 1.0, 'CH3Plus': 1.0, 'CH5Plus': 1.0,
        'ArHPlus': 1.0, 'CH3Minus': 1.0, 'C2': 0.001, 'CH': 0.001, 'H': 0.0,
        'C': 0.01, 'CH2': 0.001, 'CH3': 0.001, 'C2H': 0.001, 'C2H2': 0.001,
        'C2H3': 0.001, 'C2H4': 0.001, 'C2H5': 0.001, 'C2H6': 0.001, 'C3H2': 0.001,
        'CHPlus': 1.0, 'C3H': 0.001, 'C4H2': 0.001, 'C3H3': 0.001, 'C3H4': 0.001,
        'C3': 0.001, 'C2H4Plus': 1.0, 'C2H3Plus': 1.0, 'C2HPlus': 1.0,
        'C3H5': 0.001, 'HMinus': 1.0, 'C2H5Plus': 1.0, 'CH2Plus': 1.0,
        'C4H': 0.001, 'H3Plus': 1.0, 'H2Plus': 1.0, 'C2H6Plus': 1.0
    },
    'n_Ar': n_total * 0.97,
    'n_CH4': n_total * 0.03,
}

# Build reactions
k = define_rates(params)
params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Initialize
species = params['species']
y0 = np.ones(len(species)) * 1e3
y0[species.index('Ar')] = n_total * 0.97
y0[species.index('CH4')] = n_total * 0.03
y0[species.index('e')] = NE

# Run simulation
ode_func = PlasmaODE_Optimized(params)
ode_func.H_drift_gain = 5.7e16

print("="*80)
print("DETAILED C2 PRODUCTION ANALYSIS")
print("="*80)
print()

sol = solve_ivp(ode_func, (0, 500), y0, method='BDF',
               rtol=1e-7, atol=1e-9, max_step=1.0)

if not sol.success:
    print(f"ERROR: Solver failed: {sol.message}")
    exit(1)

y_final = sol.y[:, -1]

# Get species indices
C2_idx = species.index('C2')
n_C2 = y_final[C2_idx]

print(f"Final C2 density: {n_C2:.2e} cm⁻³")
print(f"Target C2 density: 5.6e11 cm⁻³")
print(f"Gap: {5.6e11/n_C2:.1f}× too low")
print()
print("="*80)
print()

# Analyze all reactions involving C2
production_pathways = []
destruction_pathways = []

for rxn_idx, reaction in enumerate(params['R']):
    c2_change = reaction.products[C2_idx] - reaction.reactants[C2_idx]

    if c2_change == 0:
        continue

    # Get reaction string
    reactants = []
    products = []
    for sp_idx, stoich in enumerate(reaction.reactants):
        if stoich > 0:
            if stoich == 1:
                reactants.append(species[sp_idx])
            else:
                reactants.append(f"{int(stoich)}{species[sp_idx]}")
    for sp_idx, stoich in enumerate(reaction.products):
        if stoich > 0:
            if stoich == 1:
                products.append(species[sp_idx])
            else:
                products.append(f"{int(stoich)}{species[sp_idx]}")

    rxn_string = " + ".join(reactants) + " → " + " + ".join(products)

    # Calculate rate
    rate_constant = k[params['tags'][rxn_idx]]
    rate = rate_constant
    for sp_idx, stoich in enumerate(reaction.reactants):
        if stoich > 0:
            rate *= y_final[sp_idx] ** stoich

    # Store pathway
    pathway = {
        'name': params['tags'][rxn_idx],
        'reaction': rxn_string,
        'k': rate_constant,
        'rate': rate * abs(c2_change),
        'net_rate': rate * c2_change,
        'reactant_densities': {}
    }

    # Store reactant densities
    for sp_idx, stoich in enumerate(reaction.reactants):
        if stoich > 0:
            pathway['reactant_densities'][species[sp_idx]] = y_final[sp_idx]

    if c2_change > 0:
        production_pathways.append(pathway)
    else:
        destruction_pathways.append(pathway)

# Sort by rate
production_pathways.sort(key=lambda x: x['rate'], reverse=True)
destruction_pathways.sort(key=lambda x: x['rate'], reverse=True)

# Total rates
total_production = sum(p['rate'] for p in production_pathways)
total_destruction = sum(p['rate'] for p in destruction_pathways)

print("C2 PRODUCTION PATHWAYS")
print("="*80)
print()

for i, pathway in enumerate(production_pathways[:15]):  # Top 15
    percent = pathway['rate'] / total_production * 100
    print(f"{i+1}. {pathway['reaction']}")
    print(f"   Tag: {pathway['name']}")
    print(f"   Rate constant: k = {pathway['k']:.2e} cm³/s")
    print(f"   Production rate: {pathway['rate']:.2e} cm⁻³/s ({percent:.1f}%)")

    # Show reactant densities
    for species_name, density in pathway['reactant_densities'].items():
        print(f"   [{species_name}] = {density:.2e} cm⁻³")

    print()

print(f"Total C2 production: {total_production:.2e} cm⁻³/s")
print()
print("="*80)
print()

print("C2 DESTRUCTION PATHWAYS")
print("="*80)
print()

for i, pathway in enumerate(destruction_pathways[:10]):  # Top 10
    percent = pathway['rate'] / total_destruction * 100
    print(f"{i+1}. {pathway['reaction']}")
    print(f"   Tag: {pathway['name']}")
    print(f"   Rate constant: k = {pathway['k']:.2e} cm³/s")
    print(f"   Destruction rate: {pathway['rate']:.2e} cm⁻³/s ({percent:.1f}%)")
    print()

print(f"Total C2 destruction: {total_destruction:.2e} cm⁻³/s")
print()
print("="*80)
print()

# Calculate what production rate we need
target_c2 = 5.6e11
current_c2 = n_C2
needed_factor = target_c2 / current_c2

# If we increase production by needed_factor, what happens?
# New steady state: production = destruction
# destruction scales with C2 density, so we need even more production

print("PRODUCTION/DESTRUCTION BALANCE")
print("="*80)
print()
print(f"Current C2: {current_c2:.2e} cm⁻³")
print(f"Current production: {total_production:.2e} cm⁻³/s")
print(f"Current destruction: {total_destruction:.2e} cm⁻³/s")
print(f"Net rate: {total_production - total_destruction:.2e} cm⁻³/s")
print()
print(f"To reach C2 = {target_c2:.2e} cm⁻³:")
print(f"  Need {needed_factor:.0f}× more C2")
print()

# Check precursor densities
print("="*80)
print("KEY PRECURSOR DENSITIES")
print("="*80)
print()

precursors = ['C2H2', 'C', 'CH', 'CH2', 'H', 'C2H', 'C2H3', 'e']
for sp in precursors:
    if sp in species:
        idx = species.index(sp)
        print(f"[{sp}] = {y_final[idx]:.2e} cm⁻³")

print()
print("="*80)
print()

# Identify rate-limiting steps
print("RATE-LIMITING ANALYSIS")
print("="*80)
print()

print("For each major production pathway, calculate how much it would increase")
print("if the rate constant were increased by 1263×:")
print()

for i, pathway in enumerate(production_pathways[:5]):
    current_rate = pathway['rate']
    current_k = pathway['k']

    # If we increase k by 1263×, rate increases by 1263×
    new_rate = current_rate * needed_factor

    percent_of_total = current_rate / total_production * 100
    new_total = total_production - current_rate + new_rate
    new_percent = new_rate / new_total * 100

    print(f"{i+1}. {pathway['name']}")
    print(f"   Current k: {current_k:.2e}")
    print(f"   Current rate: {current_rate:.2e} cm⁻³/s ({percent_of_total:.1f}% of total)")
    print(f"   If k × 1263: {current_k * needed_factor:.2e}")
    print(f"   New rate: {new_rate:.2e} cm⁻³/s ({new_percent:.1f}% of new total)")
    print(f"   → Would increase total production by {new_total/total_production:.0f}×")
    print()

print("="*80)
