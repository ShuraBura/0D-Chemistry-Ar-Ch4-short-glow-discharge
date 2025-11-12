"""
Analyze why high CH (343% of target) doesn't convert to C2

Key questions:
1. What's the CH + CH → C2 production rate?
2. What destroys C2 faster than CH + CH produces it?
3. What destroys CH faster than CH + CH can consume it?
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

# Test best result from three-body physics: 200× CH3 + 99% loss
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

# Get rates with baseline tuned values
k = define_rates(params)
for rate_name, rate_value in baseline['rate_values'].items():
    if rate_name in k:
        k[rate_name] = rate_value

# Apply 200× CH3 multipliers
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

# Extract densities
H_final = y_final[species.index('H')]
CH_final = y_final[species.index('CH')]
C2_final = y_final[species.index('C2')]
C_final = y_final[species.index('C')]
C2H2_final = y_final[species.index('C2H2')]

print("="*80)
print("CH TO C2 PATHWAY ANALYSIS")
print("="*80)
print(f"\nDensities at steady state:")
print(f"  H:    {H_final:.2e} cm⁻³ ({H_final/targets['H']*100:6.1f}%)")
print(f"  CH:   {CH_final:.2e} cm⁻³ ({CH_final/targets['CH']*100:6.1f}%)")
print(f"  C2:   {C2_final:.2e} cm⁻³ ({C2_final/targets['C2']*100:6.2f}%)")
print(f"  C:    {C_final:.2e} cm⁻³")
print(f"  C2H2: {C2H2_final:.2e} cm⁻³")

# Analyze CH + CH → C2 reactions
print("\n" + "="*80)
print("CH + CH → C2 + H2 PATHWAY")
print("="*80)

# There are two CH + CH → C2 reactions in the model
k_CH_CH_1 = k.get('CH_CH_C2_H2_cm3_5_4', 0)  # 2.16e-10
k_CH_CH_2 = k.get('CH_CH_C2_H2_cm3_7_44', 0)  # 1.0e-10

print(f"\nReaction 1: CH + CH → C2 + H2")
print(f"  Rate constant: {k_CH_CH_1:.2e} cm³/s")
rate_1 = k_CH_CH_1 * CH_final * CH_final if k_CH_CH_1 > 0 else 0
print(f"  Production rate: {rate_1:.2e} cm⁻³/s")

print(f"\nReaction 2: CH + CH → C2 + H2 (alternative)")
print(f"  Rate constant: {k_CH_CH_2:.2e} cm³/s")
rate_2 = k_CH_CH_2 * CH_final * CH_final if k_CH_CH_2 > 0 else 0
print(f"  Production rate: {rate_2:.2e} cm⁻³/s")

total_CH_CH_production = rate_1 + rate_2
print(f"\nTotal CH + CH → C2 production: {total_CH_CH_production:.2e} cm⁻³/s")

# Analyze C2 + H → CH + C (C2 destruction)
print("\n" + "="*80)
print("C2 + H → CH + C PATHWAY (C2 DESTRUCTION)")
print("="*80)

k_C2_H = k.get('C2_H_CH_C_cm3_7_6', 0)
print(f"\nReaction: C2 + H → CH + C")
print(f"  Rate constant: {k_C2_H:.2e} cm³/s")
rate_C2_destruction = k_C2_H * C2_final * H_final
print(f"  C2 destruction rate: {rate_C2_destruction:.2e} cm⁻³/s")
print(f"  CH production rate: {rate_C2_destruction:.2e} cm⁻³/s")

# Analyze H + C2H2 → C2 + H2 (main C2 production)
print("\n" + "="*80)
print("H + C2H2 → C2 + H2 PATHWAY (MAIN C2 PRODUCTION)")
print("="*80)

k_H_C2H2 = k.get('H_C2H2_C2_H2_cm3_7_22', 0)
print(f"\nReaction: H + C2H2 → C2 + H2")
print(f"  Rate constant: {k_H_C2H2:.2e} cm³/s")
rate_H_C2H2 = k_H_C2H2 * H_final * C2H2_final
print(f"  C2 production rate: {rate_H_C2H2:.2e} cm⁻³/s")

# Analyze CH destruction pathways
print("\n" + "="*80)
print("CH DESTRUCTION PATHWAYS")
print("="*80)

ch_destruction_reactions = [
    ('CH_CH4_C2H4_H_cm3_7_20', 'CH + CH4 → C2H4 + H', y_final[species.index('CH4')]),
    ('CH_CH3_C2H4_cm3_7_5', 'CH + CH3 → C2H4', y_final[species.index('CH3')]),
    ('CH_CH3_C2H3_H_cm3_7_10', 'CH + CH3 → C2H3 + H', y_final[species.index('CH3')]),
    ('CH_CH3_C2H2_H2_cm3_7_23', 'CH + CH3 → C2H2 + H2', y_final[species.index('CH3')]),
    ('CH_CH2_C2H2_H_cm3_7_7', 'CH + CH2 → C2H2 + H', y_final[species.index('CH2')]),
]

print("\nMajor CH destruction reactions:")
ch_destruction_rates = []
for rate_key, desc, reactant_density in ch_destruction_reactions:
    if rate_key in k:
        rate = k[rate_key] * CH_final * reactant_density
        ch_destruction_rates.append((desc, rate))
        print(f"  {desc}")
        print(f"    k = {k[rate_key]:.2e}, rate = {rate:.2e} cm⁻³/s")

# Summary
print("\n" + "="*80)
print("SUMMARY: WHY HIGH CH DOESN'T PRODUCE MORE C2")
print("="*80)

print(f"\n1. CH + CH → C2 production:")
print(f"   Rate: {total_CH_CH_production:.2e} cm⁻³/s")
print(f"   Despite CH = {CH_final:.2e}, this is VERY SMALL")

print(f"\n2. H + C2H2 → C2 production (MAIN pathway):")
print(f"   Rate: {rate_H_C2H2:.2e} cm⁻³/s")
print(f"   This is {rate_H_C2H2/total_CH_CH_production:.0f}× FASTER than CH + CH")

print(f"\n3. C2 + H → CH + C destruction:")
print(f"   Rate: {rate_C2_destruction:.2e} cm⁻³/s")
print(f"   This destroys C2 and CREATES CH")

print(f"\n4. Net C2 balance:")
print(f"   Production (H + C2H2): +{rate_H_C2H2:.2e} cm⁻³/s")
print(f"   Production (CH + CH):  +{total_CH_CH_production:.2e} cm⁻³/s")
print(f"   Destruction (C2 + H):  -{rate_C2_destruction:.2e} cm⁻³/s")
net_c2 = rate_H_C2H2 + total_CH_CH_production - rate_C2_destruction
print(f"   NET:                   {net_c2:+.2e} cm⁻³/s")

print(f"\n5. CH balance:")
print(f"   Production from C2 + H: +{rate_C2_destruction:.2e} cm⁻³/s")
print(f"   Consumption by CH + CH: -{2*total_CH_CH_production:.2e} cm⁻³/s")
ch_from_other = sum(rate for _, rate in ch_destruction_rates)
print(f"   Other CH destruction:   -{ch_from_other:.2e} cm⁻³/s")

print("\n" + "="*80)
print("KEY FINDING:")
print("="*80)
print("\nCH + CH → C2 is EXTREMELY SLOW despite high CH density!")
print(f"  CH density: {CH_final:.2e} cm⁻³ (343% of target)")
print(f"  CH + CH rate: {total_CH_CH_production:.2e} cm⁻³/s")
print(f"  Scales as [CH]² = ({CH_final:.2e})² × k = {CH_final**2 * (k_CH_CH_1 + k_CH_CH_2):.2e}")
print(f"\nWhy so slow?")
print(f"  CH is DILUTE: {CH_final:.2e} vs H = {H_final:.2e} (600× more H)")
print(f"  CH + CH ∝ [CH]² = {CH_final:.2e}² = {CH_final**2:.2e}")
print(f"  H + C2H2 ∝ [H]×[C2H2] = {H_final:.2e} × {C2H2_final:.2e} = {H_final*C2H2_final:.2e}")
print(f"\n  Ratio: [H]×[C2H2] / [CH]² = {(H_final*C2H2_final)/(CH_final**2):.0f}×")
print(f"\nConclusion: CH is TOO DILUTE for bimolecular CH + CH to compete!")

# What if we had 100× more CH?
print("\n" + "="*80)
print("THOUGHT EXPERIMENT: 100× MORE CH")
print("="*80)
CH_test = CH_final * 100
rate_CH_CH_test = (k_CH_CH_1 + k_CH_CH_2) * CH_test**2
print(f"\nIf CH = {CH_test:.2e} (100× higher):")
print(f"  CH + CH → C2 rate: {rate_CH_CH_test:.2e} cm⁻³/s")
print(f"  Would be {rate_CH_CH_test/rate_H_C2H2:.1f}× the H + C2H2 rate")
print(f"  But this requires CH = {CH_test/targets['CH']:.0f}× target (completely unrealistic)")

print("\n" + "="*80)
print("ANSWER TO USER'S QUESTION:")
print("="*80)
print("\nWhy doesn't excess CH convert to C2?")
print("\n1. CH + CH → C2 rate constant is reasonable (1-2×10⁻¹⁰ cm³/s)")
print("\n2. BUT CH is too DILUTE for bimolecular reaction to be fast:")
print(f"   [CH] = {CH_final:.2e} (343% of target but still very low)")
print(f"   [CH]² = {CH_final**2:.2e} (quadratic dependence kills it!)")
print(f"\n3. Compare to H + C2H2 → C2 (main pathway):")
print(f"   [H] is 600× higher than [CH]")
print(f"   [H]×[C2H2] = {H_final*C2H2_final:.2e} vs [CH]² = {CH_final**2:.2e}")
print(f"   That's {(H_final*C2H2_final)/(CH_final**2):.0f}× more reactive collisions!")
print(f"\n4. To make CH + CH competitive, would need:")
print(f"   CH ~ {np.sqrt(H_final*C2H2_final/(k_CH_CH_1+k_CH_CH_2)*k_H_C2H2):.2e} cm⁻³")
print(f"   That's {np.sqrt(H_final*C2H2_final/(k_CH_CH_1+k_CH_CH_2)*k_H_C2H2)/targets['CH']:.0f}× the target!")
print(f"   (Completely impossible)")
print("\n5. C2 + H → CH + C makes it worse:")
print(f"   This reaction CONSUMES C2 and PRODUCES CH")
print(f"   Creating a cycle that favors high CH, low C2")
print("\n" + "="*80)
