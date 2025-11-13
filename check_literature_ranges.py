"""
Check if 200× CH3 boost is within literature ranges

User question: Is the 200× CH3 boost achieved within literature range of rates?

We've been multiplying:
- e + CH4 → CH3 + H⁻ by 200×
- Ar* + CH4 → CH3 + H by 200×
- e + CH4⁺ → CH3 + H by 200×

Are the resulting rate constants physically plausible?
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

print("="*80)
print("CHECK IF 200× CH3 BOOST IS WITHIN LITERATURE RANGES")
print("="*80)

# Get default rates at baseline conditions
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
    'species': [],
}

k = define_rates(params)

# The three CH3 production reactions we boost
ch3_reactions = {
    'e_CH4_CH3_HMinus_cm3_8_1': 'e + CH4 → CH3 + H⁻',
    'ArStar_CH4_CH3_H_cm3_3_1': 'Ar* + CH4 → CH3 + H',
    'e_CH4Plus_CH3_H_cm3_6_4': 'e + CH4⁺ → CH3 + H',
}

print("\n1. CH3 PRODUCTION REACTIONS")
print("="*80)
print(f"\n{'Reaction':<40} {'Default k':<15} {'After 200×':<15} {'Status':<15}")
print("-"*80)

for rate_key, description in ch3_reactions.items():
    default_k = k.get(rate_key, 0)

    # Apply baseline tuning if it exists
    if rate_key in baseline.get('rate_values', {}):
        tuned_k = baseline['rate_values'][rate_key]
    else:
        tuned_k = default_k

    boosted_k = tuned_k * 200.0

    print(f"{description:<40}")
    print(f"  Default:     {default_k:.2e} cm³/s")
    print(f"  Tuned:       {tuned_k:.2e} cm³/s")
    print(f"  After 200×:  {boosted_k:.2e} cm³/s")
    print()

print("\n2. LITERATURE RANGES FOR COMPARISON")
print("="*80)

print("\nElectron impact reactions (e + molecule → products):")
print("  Typical range: 1e-13 to 1e-8 cm³/s at Te ~ 1-3 eV")
print("  Our e + CH4 → CH3 + H⁻:")
print(f"    Default:  {k['e_CH4_CH3_HMinus_cm3_8_1']:.2e} cm³/s")
print(f"    After 200×: {k['e_CH4_CH3_HMinus_cm3_8_1'] * 200:.2e} cm³/s")

if k['e_CH4_CH3_HMinus_cm3_8_1'] * 200 > 1e-8:
    print("    ⚠️  WARNING: Exceeds typical upper limit (1e-8 cm³/s)")
elif k['e_CH4_CH3_HMinus_cm3_8_1'] * 200 < 1e-13:
    print("    ⚠️  WARNING: Below typical lower limit (1e-13 cm³/s)")
else:
    print("    ✓ Within typical literature range")

print("\nMetastable reactions (Ar* + molecule → products):")
print("  Typical range: 1e-11 to 1e-9 cm³/s")
print("  Our Ar* + CH4 → CH3 + H:")
print(f"    Default:  {k['ArStar_CH4_CH3_H_cm3_3_1']:.2e} cm³/s")
print(f"    After 200×: {k['ArStar_CH4_CH3_H_cm3_3_1'] * 200:.2e} cm³/s")

if k['ArStar_CH4_CH3_H_cm3_3_1'] * 200 > 1e-9:
    print("    ⚠️  WARNING: Exceeds typical upper limit (1e-9 cm³/s)")
elif k['ArStar_CH4_CH3_H_cm3_3_1'] * 200 < 1e-11:
    print("    ⚠️  WARNING: Below typical lower limit (1e-11 cm³/s)")
else:
    print("    ✓ Within typical literature range")

print("\nDissociative recombination (e + ion → neutrals):")
print("  Typical range: 1e-8 to 1e-6 cm³/s at Te ~ 1 eV")
print("  Our e + CH4⁺ → CH3 + H:")
print(f"    Default:  {k['e_CH4Plus_CH3_H_cm3_6_4']:.2e} cm³/s")
print(f"    After 200×: {k['e_CH4Plus_CH3_H_cm3_6_4'] * 200:.2e} cm³/s")

if k['e_CH4Plus_CH3_H_cm3_6_4'] * 200 > 1e-6:
    print("    ⚠️  WARNING: Exceeds typical upper limit (1e-6 cm³/s)")
elif k['e_CH4Plus_CH3_H_cm3_6_4'] * 200 < 1e-8:
    print("    ⚠️  WARNING: Below typical lower limit (1e-8 cm³/s)")
else:
    print("    ✓ Within typical literature range")

print("\n" + "="*80)
print("3. PHYSICAL INTERPRETATION")
print("="*80)

print("\nWhat does 200× boost ACTUALLY mean?")
print("\nOption A: Rate constant multiplier")
print("  - We're saying k_effective = 200 × k_literature")
print("  - This would be NON-PHYSICAL (violates quantum mechanics)")
print("  - Rate constants are determined by cross-sections and collision theory")
print("")
print("Option B: Enhanced production via other mechanisms")
print("  - 200× could represent:")
print("    * Higher electron density")
print("    * More Ar* metastables")
print("    * Additional production pathways not in model")
print("    * Surface reactions producing CH3")
print("  - This is more realistic interpretation")

print("\n" + "="*80)
print("4. RECOMMENDATION")
print("="*80)

print("\nThe 200× 'boost' is likely NOT realistic as a rate constant multiplier.")
print("\nInstead, consider it represents:")
print("  1. Missing CH3 production pathways in the model")
print("  2. Surface-assisted CH3 generation")
print("  3. Enhanced metastable production not captured")
print("  4. Uncertainties in the baseline rate constants")
print("\nTo be physically realistic, we should:")
print("  - Keep rate constants within literature ranges")
print("  - Add missing production pathways explicitly")
print("  - Model surface chemistry if relevant")
print("  - Increase electron/metastable densities if justified")

print("\n" + "="*80)
print("ANSWER TO USER'S QUESTION:")
print("="*80)
print("\nIs 200× CH3 boost within literature range?")
print("\nNO - if interpreted as rate constant multiplier:")

for rate_key, description in ch3_reactions.items():
    default_k = k.get(rate_key, 0)
    if rate_key in baseline.get('rate_values', {}):
        tuned_k = baseline['rate_values'][rate_key]
    else:
        tuned_k = default_k
    boosted_k = tuned_k * 200.0

    print(f"\n  {description}")
    print(f"    Boosted k: {boosted_k:.2e} cm³/s")

    # Check against typical ranges
    if rate_key == 'e_CH4_CH3_HMinus_cm3_8_1':
        if boosted_k > 1e-8:
            print(f"    ⚠️  Exceeds typical e-impact range (1e-13 to 1e-8)")
        else:
            print(f"    ✓ Within typical range")
    elif rate_key == 'ArStar_CH4_CH3_H_cm3_3_1':
        if boosted_k > 1e-9:
            print(f"    ⚠️  Exceeds typical metastable range (1e-11 to 1e-9)")
        else:
            print(f"    ✓ Within typical range")
    elif rate_key == 'e_CH4Plus_CH3_H_cm3_6_4':
        if boosted_k > 1e-6:
            print(f"    ⚠️  Exceeds typical dissoc. recomb. range (1e-8 to 1e-6)")
        else:
            print(f"    ✓ Within typical range")

print("\n" + "="*80)
