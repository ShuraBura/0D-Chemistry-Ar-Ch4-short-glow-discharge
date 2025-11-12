"""
Check current ne density and list all CH3 production pathways in the model
"""

import json
from define_rates import define_rates

def pressure_to_density(pressure_mTorr, T_K=300):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

# Load baseline
with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

P = 500.0
n_total = pressure_to_density(P)
ne_frac = baseline['Ne'] / pressure_to_density(500.0)
ne = ne_frac * n_total

print("="*80)
print("CURRENT ELECTRON DENSITY")
print("="*80)
print(f"\nPressure: {P} mTorr")
print(f"Total density: {n_total:.2e} cm⁻³")
print(f"Electron density (ne): {ne:.2e} cm⁻³")
print(f"Electron fraction: {ne/n_total:.2e} ({ne/n_total*100:.4f}%)")

print("\n" + "="*80)
print("ALL CH3 PRODUCTION PATHWAYS IN MODEL")
print("="*80)

# Get rate constants
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
    'species': [],
}

k = define_rates(params)

# Find all reactions producing CH3
ch3_production = {
    'e_CH4_CH3_HMinus_cm3_8_1': ('e + CH4 → CH3 + H⁻', 'electron impact'),
    'ArStar_CH4_CH3_H_cm3_3_1': ('Ar* + CH4 → CH3 + H', 'metastable'),
    'e_CH4Plus_CH3_H_cm3_6_4': ('e + CH4⁺ → CH3 + H', 'dissoc. recomb.'),
    'e_C2H5_CH3_CH2_cm3_1_5': ('e + C2H5 → CH3 + CH2', 'electron impact'),
    'e_C2H6_CH3_CH3_cm3_1_7': ('e + C2H6 → CH3 + CH3', 'electron impact'),
    'ArStar_C2H6_CH3_CH3_cm3_3_6': ('Ar* + C2H6 → CH3 + CH3', 'metastable'),
    'C2H5Plus_e_CH3_CH2_cm3_6_10': ('C2H5⁺ + e → CH3 + CH2', 'dissoc. recomb.'),
    'CH_C2H6_C2H2_CH3_H_cm3_7_57': ('CH + C2H6 → C2H2 + CH3 + H', 'neutral-neutral'),
}

print("\nReaction pathways currently in model:")
print(f"\n{'Reaction':<50} {'Type':<20} {'Rate (cm³/s)':<15}")
print("-"*90)

for key, (reaction, rxn_type) in ch3_production.items():
    if key in k:
        rate = k[key]
        print(f"{reaction:<50} {rxn_type:<20} {rate:.2e}")
    else:
        print(f"{reaction:<50} {rxn_type:<20} NOT FOUND")

print("\n" + "="*80)
print("POTENTIALLY MISSING CH3 PRODUCTION PATHWAYS")
print("="*80)

print("\n1. Additional Electron Impact Reactions:")
print("   - e + C2H4 → CH3 + CH")
print("   - e + C2H5 → CH3 + CH2 (already included)")
print("   - e + C3H8 → CH3 + C2H5")
print("   - e + C3H6 → CH3 + C2H3")

print("\n2. Additional Metastable Reactions:")
print("   - Ar* + C2H4 → CH3 + CH")
print("   - Ar* + C2H5 → CH3 + CH2")

print("\n3. Ion-Neutral Reactions:")
print("   - H⁺ + CH4 → CH3⁺ + H2 (then CH3⁺ + e → CH3)")
print("   - Ar⁺ + CH4 → CH3⁺ + Ar + H")

print("\n4. Neutral-Neutral Reactions:")
print("   - H + C2H5 → CH3 + CH2")
print("   - CH + CH4 → CH3 + CH (H-abstraction)")
print("   - CH2 + CH2 → CH3 + CH")

print("\n5. Surface/Wall Reactions:")
print("   - CH4 + surface → CH3 (ads) + H (ads)")
print("   - CH3 (ads) → CH3 (gas) (desorption)")
print("   - Surface-assisted dissociation")

print("\n6. Three-Body Reactions:")
print("   - CH2 + H + M → CH3 + M")
print("   - CH + H2 + M → CH3 + M")

print("\n" + "="*80)
print("RECOMMENDATION: ADD HIGH-PRIORITY MISSING PATHWAYS")
print("="*80)

print("\nBased on literature, add these reactions:")
print("\n1. e + C2H4 → CH3 + CH")
print("   Typical k: 1-5 × 10⁻¹¹ cm³/s at Te ~ 1.5 eV")
print("   Source: C2H4 is abundant intermediate")

print("\n2. H + C2H5 → CH3 + CH2")
print("   Typical k: 5-10 × 10⁻¹¹ cm³/s")
print("   Source: H is very abundant (2×10¹⁴ cm⁻³)")

print("\n3. CH2 + H + M → CH3 + M")
print("   Typical k: 1-5 × 10⁻³⁰ cm⁶/s")
print("   Source: Three-body stabilization")

print("\n4. Ar⁺ + CH4 → CH3⁺ + Ar + H (Penning ionization)")
print("   Typical k: 1-5 × 10⁻⁹ cm³/s")
print("   Then: CH3⁺ + e → CH3 (fast dissoc. recomb.)")

print("\n" + "="*80)
