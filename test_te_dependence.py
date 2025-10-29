"""
Test script to verify Te-dependent rates implementation
"""
import numpy as np
from define_rates_tunable import define_rates_tunable

# Test parameters
test_params = {
    'L_discharge': 0.4,  # cm
    'mobilities': {
        'ArPlus': 1.5, 'CH4Plus': 2.3, 'CH3Plus': 2.3, 'CH5Plus': 2.3,
        'ArHPlus': 2.0, 'CH2Plus': 2.5, 'C2H5Plus': 2.0, 'C2H4Plus': 2.0,
        'C2H3Plus': 2.0, 'H3Plus': 2.5, 'CHPlus': 2.5, 'CH3Minus': 2.0,
        'C2HPlus': 2.0
    },
    'Tgas': 400,
    'L_diff': 0.1,
    'scale_e_impact': 1.0,
}

# Test different electron temperatures
Te_values = [0.5, 1.0, 2.0, 5.0, 7.0]

print("Testing Te-dependent rate implementation")
print("=" * 80)

# Test a few representative reactions
test_reactions = [
    ('e_CH4_CH3_H_cm3_1_1', 'Electron-impact dissociation'),
    ('e_Ar_ArPlus_cm3_2_3', 'Electron-impact ionization'),
    ('ArPlus_e_Ar_cm3_6_1', 'Dissociative recombination'),
    ('CH_CH3_C2H4_cm3_7_5', 'Neutral-neutral (should be constant)')
]

results = {}
for reaction_name, description in test_reactions:
    results[reaction_name] = []
    print(f"\n{description}: {reaction_name}")
    print(f"{'Te (eV)':<10} {'Rate (cm³/s)':<15} {'Relative to Te=1'}")
    print("-" * 45)

    ref_rate = None
    for Te in Te_values:
        params = test_params.copy()
        params['Te'] = Te
        k = define_rates_tunable(params)
        rate = k.get(reaction_name, 0)
        results[reaction_name].append(rate)

        # Calculate relative to Te=1 eV
        if Te == 1.0:
            ref_rate = rate
        relative = rate / ref_rate if ref_rate and ref_rate != 0 else 1.0

        print(f"{Te:<10.1f} {rate:<15.3e} {relative:.3f}")

print("\n" + "=" * 80)
print("VALIDATION CHECKS:")
print("-" * 80)

# Check 1: Electron-impact rates should increase with Te
e_impact_rate_05 = results['e_CH4_CH3_H_cm3_1_1'][0]
e_impact_rate_70 = results['e_CH4_CH3_H_cm3_1_1'][-1]
if e_impact_rate_70 > e_impact_rate_05:
    print("✓ Electron-impact dissociation rates increase with Te")
else:
    print("✗ ERROR: Electron-impact dissociation rates should increase with Te")

# Check 2: Ionization rates should increase with Te
ion_rate_05 = results['e_Ar_ArPlus_cm3_2_3'][0]
ion_rate_70 = results['e_Ar_ArPlus_cm3_2_3'][-1]
if ion_rate_70 > ion_rate_05:
    print("✓ Ionization rates increase with Te")
else:
    print("✗ ERROR: Ionization rates should increase with Te")

# Check 3: Recombination rates should decrease with Te
rec_rate_05 = results['ArPlus_e_Ar_cm3_6_1'][0]
rec_rate_70 = results['ArPlus_e_Ar_cm3_6_1'][-1]
if rec_rate_70 < rec_rate_05:
    print("✓ Recombination rates decrease with Te")
else:
    print("✗ ERROR: Recombination rates should decrease with Te")

# Check 4: Neutral-neutral rates should be constant
neutral_rates = results['CH_CH3_C2H4_cm3_7_5']
if all(abs(r - neutral_rates[0]) < 1e-20 for r in neutral_rates):
    print("✓ Neutral-neutral rates are constant (independent of Te)")
else:
    print("✗ ERROR: Neutral-neutral rates should be independent of Te")

print("\n" + "=" * 80)
print("Test completed successfully!")
print("=" * 80)
