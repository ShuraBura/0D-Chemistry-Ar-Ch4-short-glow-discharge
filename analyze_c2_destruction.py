#!/usr/bin/env python3
"""
Analyze C2 DESTRUCTION pathways with corrected chemistry
We removed C2 + H → CH + C (93% of destruction) - what's left?
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
        'ArHPlus': 1.0, 'CH3Minus': 1.0, 'C2': 0.01, 'CH': 0.001, 'H': 0.0,
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

print("="*80)
print("C2 DESTRUCTION PATHWAY ANALYSIS")
print("="*80)
print()

# Build reactions
k = define_rates(params)
params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Check if C2 + H was disabled
c2_h_rate = k.get('C2_H_CH_C_cm3_7_6', None)
print(f"C2 + H → CH + C rate: {c2_h_rate} (should be 0.0)")
print()

# Initialize
species = params['species']
y0 = np.ones(len(species)) * 1e3
y0[species.index('Ar')] = n_total * 0.97
y0[species.index('CH4')] = n_total * 0.03
y0[species.index('e')] = NE

# Run simulation
ode_func = PlasmaODE_Optimized(params)
ode_func.H_drift_gain = 5.7e16

print("Solving for steady state...")
sol = solve_ivp(ode_func, (0, 500), y0, method='BDF',
               rtol=1e-7, atol=1e-9, max_step=1.0)

if not sol.success:
    print(f"❌ Solver failed: {sol.message}")
    exit(1)

# Get final densities
y_final = sol.y[:, -1]

# Get C2 index
C2_idx = species.index('C2')
n_C2 = y_final[C2_idx]

print(f"✓ Converged!")
print(f"  C2 = {n_C2:.2e} cm⁻³")
print()

# Analyze C2 production and destruction
print("="*80)
print("C2 PRODUCTION & DESTRUCTION BREAKDOWN")
print("="*80)
print()

production_pathways = []
destruction_pathways = []

for rxn_idx, reaction in enumerate(params['R']):
    # Check if C2 is involved
    c2_change = reaction.products[C2_idx] - reaction.reactants[C2_idx]

    if c2_change == 0:
        continue  # C2 not involved

    # Get rate constant
    rate_constant = params['k'][params['tags'][rxn_idx]]

    # Calculate reaction rate: k * [A] * [B] * ...
    rate = rate_constant
    reactant_names = []
    for sp_idx, stoich in enumerate(reaction.reactants):
        if stoich > 0:
            rate *= y_final[sp_idx] ** stoich
            reactant_names.append(f"{species[sp_idx]}^{int(stoich)}" if stoich > 1 else species[sp_idx])

    # Get product names
    product_names = []
    for sp_idx, stoich in enumerate(reaction.products):
        if stoich > 0:
            product_names.append(f"{species[sp_idx]}^{int(stoich)}" if stoich > 1 else species[sp_idx])

    rxn_string = " + ".join(reactant_names) + " → " + " + ".join(product_names)

    if c2_change > 0:
        # C2 production
        production_pathways.append({
            'name': params['tags'][rxn_idx],
            'reaction': rxn_string,
            'rate': rate * c2_change,
            'k': rate_constant
        })
    else:
        # C2 destruction
        destruction_pathways.append({
            'name': params['tags'][rxn_idx],
            'reaction': rxn_string,
            'rate': rate * abs(c2_change),
            'k': rate_constant
        })

# Sort by rate
production_pathways.sort(key=lambda x: x['rate'], reverse=True)
destruction_pathways.sort(key=lambda x: x['rate'], reverse=True)

# Calculate totals
total_prod = sum(p['rate'] for p in production_pathways)
total_dest = sum(p['rate'] for p in destruction_pathways)

print(f"Total C2 production:  {total_prod:.2e} cm⁻³/s")
print(f"Total C2 destruction: {total_dest:.2e} cm⁻³/s")
print(f"Net C2 rate:          {total_prod - total_dest:.2e} cm⁻³/s")
print(f"Lifetime (τ = [C2] / destruction): {n_C2/total_dest:.2e} s")
print()

print("="*80)
print("TOP 10 C2 DESTRUCTION PATHWAYS")
print("="*80)
print()
print(f"{'#':<3} {'Reaction':<60} {'Rate (cm⁻³/s)':<15} {'%':<6}")
print("-"*80)

for i, pathway in enumerate(destruction_pathways[:10], 1):
    pct = 100 * pathway['rate'] / total_dest if total_dest > 0 else 0
    print(f"{i:<3} {pathway['reaction']:<60} {pathway['rate']:<15.2e} {pct:<6.2f}")

print()
print("="*80)
print("TOP 10 C2 PRODUCTION PATHWAYS")
print("="*80)
print()
print(f"{'#':<3} {'Reaction':<60} {'Rate (cm⁻³/s)':<15} {'%':<6}")
print("-"*80)

for i, pathway in enumerate(production_pathways[:10], 1):
    pct = 100 * pathway['rate'] / total_prod if total_prod > 0 else 0
    print(f"{i:<3} {pathway['reaction']:<60} {pathway['rate']:<15.2e} {pct:<6.2f}")

print()
print("="*80)
print("KEY FINDINGS")
print("="*80)
print()

# Find the C2 + H reaction
c2_h_found = False
for pathway in destruction_pathways:
    if 'C2_H_CH_C' in pathway['name']:
        c2_h_found = True
        print(f"⚠️ WARNING: C2 + H → CH + C is STILL ACTIVE!")
        print(f"   Rate: {pathway['rate']:.2e} cm⁻³/s ({100*pathway['rate']/total_dest:.1f}%)")
        break

if not c2_h_found:
    print("✓ C2 + H → CH + C successfully disabled")

print()

# Identify wall loss
wall_dest = 0
for pathway in destruction_pathways:
    if 'stick_C2' in pathway['name'] or 'loss_C2' in pathway['name']:
        wall_dest += pathway['rate']

print(f"Wall & diffusion losses: {wall_dest:.2e} cm⁻³/s ({100*wall_dest/total_dest:.1f}% of destruction)")
print()

# Chemical destruction
chem_dest = total_dest - wall_dest
print(f"Chemical destruction:    {chem_dest:.2e} cm⁻³/s ({100*chem_dest/total_dest:.1f}% of destruction)")
print()

print("="*80)
print("COMPARISON WITH OLD MODEL")
print("="*80)
print()
print("OLD MODEL (with C2 + H → CH + C):")
print("  C2 = 2.75×10⁸ cm⁻³")
print("  C2 + H destruction was 93% of total")
print()
print("NEW MODEL (without C2 + H → CH + C):")
print(f"  C2 = {n_C2:.2e} cm⁻³")
if len(destruction_pathways) > 0:
    print(f"  Top destruction: {destruction_pathways[0]['reaction'][:50]}")
    print(f"                   {100*destruction_pathways[0]['rate']/total_dest:.1f}% of total")
print()

if n_C2 < 2.75e8:
    print("⚠️ C2 DECREASED after removing fake destruction!")
    print("   This means production also decreased (lost fake pathways)")
    print()
    print("Lost production pathways:")
    print("  - CH + CH2 → C2 (was 21%)")
    print("  - CH2 + CH2 → C2 (was 16%, now makes C2H2)")
    print("  - CH + CH → C2 (now makes C2H2)")
    print("  - H + C2H2 → C2 (rate ~99.9% slower)")
    print()
    print("Net effect: Lost more production than destruction!")

print("="*80)
