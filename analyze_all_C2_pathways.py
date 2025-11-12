"""
Comprehensive analysis of ALL C2 production and destruction pathways

Goal: Identify every pathway that affects C2 to find opportunities for improvement
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

def pressure_to_density(pressure_mTorr, T_K=300):
    """Convert pressure to total density (cm⁻³)"""
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23  # J/K
    n_m3 = pressure_Pa / (k_B * T_K)  # m⁻³
    return n_m3 * 1e-6  # Convert to cm⁻³

# Load baseline
with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

targets = {
    'H': 2.52e14,
    'CH': 1.0e9,
    'C2': 5.6e11,
}

# Run best configuration to get steady state
P = 500.0
n_total = pressure_to_density(P)
ne_frac = baseline['Ne'] / pressure_to_density(500.0)
ne = ne_frac * n_total

params = {
    'P': P,
    'Te': baseline['Te'],
    'ne': ne,
    'E_field': baseline['E_field'],
    'L_discharge': 0.45,
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
    },
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH',
                'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C',
                'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H',
                'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus',
                'H2Plus', 'C2H2Star'],
}

# Get rates with best configuration
k = define_rates(params)
for rate_name, rate_value in baseline['rate_values'].items():
    if rate_name in k:
        k[rate_name] = rate_value

# Apply best multipliers (WITHOUT C2+H elimination to see full picture)
rate_multipliers = {
    'e_CH4_CH3_HMinus_cm3_8_1': 200.0,
    'ArStar_CH4_CH3_H_cm3_3_1': 200.0,
    'e_CH4Plus_CH3_H_cm3_6_4': 200.0,
    'stick_CH3_9_2': 0.01,
    'stick_C2H2_9_11': 0.01,
    'loss_C2H2_11_19': 0.01,
}

for rate_name, mult in rate_multipliers.items():
    if rate_name in k:
        k[rate_name] *= mult

params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Run simulation
species = params['species']
y0 = np.ones(len(species)) * 1e3
y0[species.index('Ar')] = n_total * 0.85
y0[species.index('CH4')] = n_total * 0.15
y0[species.index('e')] = ne

ode_func = PlasmaODE_Optimized(params)
sol = solve_ivp(
    ode_func, (0, 500), y0,
    method='BDF', rtol=1e-7, atol=1e-9, max_step=1.0
)

y_final = sol.y[:, -1]

# Get final densities
densities = {sp: y_final[i] for i, sp in enumerate(species)}

print("="*80)
print("COMPREHENSIVE C2 PATHWAY ANALYSIS")
print("="*80)
print(f"\nCurrent state (200× CH3, 99% loss reduction):")
print(f"  C2: {densities['C2']:.2e} cm⁻³ ({densities['C2']/targets['C2']*100:.2f}%)")
print(f"  H:  {densities['H']:.2e} cm⁻³ ({densities['H']/targets['H']*100:.1f}%)")
print(f"  CH: {densities['CH']:.2e} cm⁻³ ({densities['CH']/targets['CH']*100:.1f}%)")

# Search for all C2-related reactions in rate constants
print("\n" + "="*80)
print("ALL C2 PRODUCTION REACTIONS")
print("="*80)

c2_production = []

# Scan through all rate constants for C2 production
c2_production_reactions = {
    # Direct C2 formation
    'CH_CH_C2_H2_cm3_5_4': ('CH + CH → C2 + H2', ['CH', 'CH']),
    'CH_CH_C2_H2_cm3_7_44': ('CH + CH → C2 + H2 (alt)', ['CH', 'CH']),
    'C_C_M_C2_M_cm6_7_64': ('C + C + M → C2 + M', ['C', 'C']),
    'H_C2H2_C2_H2_cm3_7_22': ('H + C2H2 → C2 + H2', ['H', 'C2H2']),
    'C_C2H_C2_CH_cm3_7_53': ('C + C2H → C2 + CH', ['C', 'C2H']),
    'C_C2H2_C2_CH2_cm3_7_16': ('C + C2H2 → C2 + CH2', ['C', 'C2H2']),
    'C_C2H3_C2_CH3_cm3_7_17': ('C + C2H3 → C2 + CH3', ['C', 'C2H3']),
}

print("\nReaction                                Rate Const      Reactant Dens      Rate (cm⁻³/s)")
print("-"*80)

for rate_key, (desc, reactants) in c2_production_reactions.items():
    if rate_key in k and k[rate_key] > 0:
        rate_const = k[rate_key]
        if len(reactants) == 2 and reactants[0] == reactants[1]:
            # Bimolecular with same species: rate = k * n²
            n = densities.get(reactants[0], 0)
            rate = rate_const * n * n
            print(f"{desc:<40} {rate_const:>10.2e}  [{reactants[0]}]²={n:.2e}  {rate:>10.2e}")
        elif len(reactants) == 2:
            # Bimolecular: rate = k * n1 * n2
            n1 = densities.get(reactants[0], 0)
            n2 = densities.get(reactants[1], 0)
            rate = rate_const * n1 * n2
            print(f"{desc:<40} {rate_const:>10.2e}  [{reactants[0]}]×[{reactants[1]}]  {rate:>10.2e}")
        else:
            rate = 0
            print(f"{desc:<40} {rate_const:>10.2e}  (complex)        {rate:>10.2e}")

        c2_production.append((desc, rate, rate_key))

# Sort by rate
c2_production.sort(key=lambda x: x[1], reverse=True)

total_production = sum(r[1] for r in c2_production)
print(f"\n{'TOTAL C2 PRODUCTION:':<40} {total_production:>40.2e} cm⁻³/s")

# Now destruction
print("\n" + "="*80)
print("ALL C2 DESTRUCTION REACTIONS")
print("="*80)

c2_destruction = []

c2_destruction_reactions = {
    'C2_H_CH_C_cm3_7_6': ('C2 + H → CH + C', ['C2', 'H']),
    'e_C2_C_C_e_cm3_1_24': ('e + C2 → C + C + e', ['e', 'C2']),
    'C2_CH_C_C_CH_cm3_7_8': ('C2 + CH → C + C + CH', ['C2', 'CH']),
    'C2_CH2_C_C_CH2_cm3_7_9': ('C2 + CH2 → C + C + CH2', ['C2', 'CH2']),
    'C2_CH3_C_C_CH3_cm3_7_11': ('C2 + CH3 → C + C + CH3', ['C2', 'CH3']),
    'C2_C2H_C_C_C2H_cm3_7_54': ('C2 + C2H → C + C + C2H', ['C2', 'C2H']),
    'C2_H2_C2H_H_cm3_7_21': ('C2 + H2 → C2H + H', ['C2', 'H2']),
    'ArStar_C2_C_C_cm3_3_8': ('Ar* + C2 → C + C', ['ArStar', 'C2']),
    'stick_C2_9_15': ('C2 wall loss', ['C2']),
}

print("\nReaction                                Rate Const      Reactant Dens      Rate (cm⁻³/s)")
print("-"*80)

for rate_key, (desc, reactants) in c2_destruction_reactions.items():
    if rate_key in k and k[rate_key] > 0:
        rate_const = k[rate_key]
        if len(reactants) == 1:
            # Loss term
            n = densities.get(reactants[0], 0)
            rate = rate_const * n
            print(f"{desc:<40} {rate_const:>10.2e}  [{reactants[0]}]={n:.2e}  {rate:>10.2e}")
        elif len(reactants) == 2:
            n1 = densities.get(reactants[0], 0)
            n2 = densities.get(reactants[1], 0)
            rate = rate_const * n1 * n2
            print(f"{desc:<40} {rate_const:>10.2e}  [{reactants[0]}]×[{reactants[1]}]  {rate:>10.2e}")
        else:
            rate = 0

        c2_destruction.append((desc, rate, rate_key))

c2_destruction.sort(key=lambda x: x[1], reverse=True)

total_destruction = sum(r[1] for r in c2_destruction)
print(f"\n{'TOTAL C2 DESTRUCTION:':<40} {total_destruction:>40.2e} cm⁻³/s")

# Net balance
print("\n" + "="*80)
print("C2 BALANCE")
print("="*80)
print(f"\nTotal C2 production:  {total_production:>15.2e} cm⁻³/s")
print(f"Total C2 destruction: {total_destruction:>15.2e} cm⁻³/s")
print(f"NET C2 rate:          {total_production - total_destruction:>+15.2e} cm⁻³/s")

# Rank opportunities
print("\n" + "="*80)
print("OPTIMIZATION OPPORTUNITIES")
print("="*80)

print("\nTOP C2 PRODUCTION REACTIONS TO BOOST:")
print("-"*80)
for i, (desc, rate, key) in enumerate(c2_production[:5], 1):
    print(f"{i}. {desc:<45} {rate:>12.2e} cm⁻³/s")
    print(f"   Rate key: {key}")
    print(f"   Current k: {k[key]:.2e}")

print("\nTOP C2 DESTRUCTION REACTIONS TO SUPPRESS:")
print("-"*80)
for i, (desc, rate, key) in enumerate(c2_destruction[:5], 1):
    print(f"{i}. {desc:<45} {rate:>12.2e} cm⁻³/s")
    print(f"   Rate key: {key}")
    print(f"   Current k: {k[key]:.2e}")

# Calculate potential improvements
print("\n" + "="*80)
print("POTENTIAL IMPROVEMENTS")
print("="*80)

print("\nIf we boost H + C2H2 → C2 by 10×:")
for desc, rate, key in c2_production:
    if 'C2H2' in desc and 'H +' in desc:
        new_rate = rate * 10
        improvement = new_rate / total_production
        print(f"  Current:  {rate:.2e} cm⁻³/s")
        print(f"  New:      {new_rate:.2e} cm⁻³/s")
        print(f"  Production increases by {improvement:.1f}×")
        break

print("\nIf we suppress C2 + H → CH + C by 100×:")
for desc, rate, key in c2_destruction:
    if desc == 'C2 + H → CH + C':
        new_rate = rate * 0.01
        saved = rate - new_rate
        print(f"  Current:  {rate:.2e} cm⁻³/s")
        print(f"  New:      {new_rate:.2e} cm⁻³/s")
        print(f"  Saves:    {saved:.2e} cm⁻³/s")
        break

# Key recommendations
print("\n" + "="*80)
print("KEY RECOMMENDATIONS FOR C2 IMPROVEMENT")
print("="*80)

print("\n1. BOOST PRODUCTION:")
for desc, rate, key in c2_production[:3]:
    if rate > 1e9:
        print(f"   ✓ {desc}")
        print(f"     Rate key: {key}")

print("\n2. SUPPRESS DESTRUCTION:")
for desc, rate, key in c2_destruction[:3]:
    if rate > 1e12:
        print(f"   ✓ {desc}")
        print(f"     Rate key: {key}")

print("\n3. BOOST C2 PRECURSORS:")
print("   ✓ C2H2 production (currently limited by collision frequency)")
print("   ✓ C2H species (C2H, C2H3) that could convert to C2")

print("\n" + "="*80)
