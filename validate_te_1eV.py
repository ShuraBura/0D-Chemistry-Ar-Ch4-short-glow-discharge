"""
Validation script: Compare rates from define_rates.py with define_rates_tunable.py at Te=1 eV
They should be EXACTLY the same (with scale_e_impact=1.0)
"""
from define_rates import define_rates
from define_rates_tunable import define_rates_tunable

# Test parameters
base_params = {
    'E_field': 800,
    'L_discharge': 0.4,
    'mobilities': {
        'ArPlus': 1.5, 'CH4Plus': 2.3, 'CH3Plus': 2.3, 'CH5Plus': 2.3,
        'ArHPlus': 2.0, 'CH2Plus': 2.5, 'C2H5Plus': 2.0, 'C2H4Plus': 2.0,
        'C2H3Plus': 2.0, 'H3Plus': 2.5, 'CHPlus': 2.5, 'CH3Minus': 2.0,
        'C2HPlus': 2.0
    },
}

# Get reference rates from original define_rates
k_original = define_rates(base_params)

# Get rates from Te-dependent version with Te=1.0 eV
tunable_params = base_params.copy()
tunable_params.update({
    'Te': 1.0,
    'Tgas': 400,
    'L_diff': 0.1,
    'scale_e_impact': 1.0,
    'P': 0.4,
})
k_tunable = define_rates_tunable(tunable_params)

# Groups to check
groups_to_check = {
    'Group 1: Electron-Impact Excitation/Dissociation': [
        'e_CH4_CH3_H_cm3_1_1', 'e_CH4_CH2_H2_cm3_1_2', 'e_CH4_CH_H2_H_vib_cm3_1_3',
        'e_H2_H_H_cm3_1_4', 'e_CH3_CH2_H_cm3_1_5', 'e_C2H4_C2H2_H2_cm3_1_6',
        'e_Ar_ArStar_cm3_1_7', 'e_C2H6_C2H4_H2_cm3_1_8', 'e_C2H6_C2H4_H2_e_cm3_1_9',
        'e_CH4_CH3Minus_H_cm3_1_10', 'e_CH4_CH_H_H2_cm3_1_11', 'e_CH_CH_C_H_e_cm3_1_12',
        'e_H2_HMinus_H_cm3_1_13', 'e_CH3_CH3Minus_cm3_1_14', 'e_C2H4_C2H2_H2_cm3_1_15',
        'e_C2H2_C2_H2_cm3_1_16', 'e_C2H4_C2H2_H_H_cm3_1_17', 'e_C2H6_C2H2_2H2_cm3_1_18',
        'e_C2H2_C2H_H_cm3_1_19', 'e_C2H4_C2H3_H_cm3_1_20', 'e_C2H6_C2H5_H_cm3_1_21',
    ],
    'Group 2: Electron-Impact Ionization': [
        'e_CH4_CH3Plus_H_cm3_2_1', 'e_CH4_CH4Plus_cm3_2_2', 'e_Ar_ArPlus_cm3_2_3',
        'e_ArStar_ArPlus_cm3_2_4', 'e_C2H6_C2H5Plus_H_2e_cm3_2_5', 'e_C2H4_C2H4Plus_2e_cm3_2_6',
        'e_C2H4_C2H3Plus_H_2e_cm3_2_7', 'e_C2H2_C2HPlus_2e_cm3_2_8',
    ],
    'Group 6: Dissociative Recombination': [
        'ArPlus_e_Ar_cm3_6_1', 'CH3Plus_e_CH3_cm3_6_2', 'CH5Plus_e_CH4_H_cm3_6_3',
        'e_CH4Plus_CH3_H_cm3_6_4', 'CH3Minus_ArPlus_CH3_Ar_cm3_6_5', 'CH3Minus_CH4Plus_CH4_CH3_cm3_6_6',
        'CH3Minus_CH3Plus_CH4_CH2_cm3_6_7', 'CH5Plus_e_CH3_H2_cm3_6_8', 'e_CH4Plus_CH2_H2_cm3_6_9',
        'CH5Plus_e_CH2_H2_H_cm3_6_10', 'e_CH4Plus_CH_H2_H_cm3_6_11', 'CH5Plus_e_CH3_2H_cm3_6_12',
        'e_CH4Plus_C_2H2_cm3_6_13', 'C2H5Plus_e_C2H4_H_cm3_6_14', 'C2H4Plus_e_C2H2_H2_cm3_6_15',
        'C2H3Plus_e_C2H2_H_cm3_6_16', 'HMinus_ArPlus_H_Ar_cm3_6_17', 'C2HPlus_e_C2_H_cm3_6_18',
        'HMinus_CH5Plus_CH4_H2_H_cm3_6_19', 'CH4Plus_HMinus_CH4_H_cm3_6_20', 'CH3Plus_HMinus_CH4_H2_cm3_6_21',
        'C2H5Plus_HMinus_C2H6_H_cm3_6_22', 'ArHPlus_HMinus_Ar_H2_H_cm3_6_23', 'CH5Plus_CH3Minus_CH4_CH4_H_cm3_6_24',
        'CH4Plus_CH3Minus_CH4_CH3_H_cm3_6_25', 'CH3Plus_CH3Minus_CH4_CH2_H_cm3_6_26', 'C2H5Plus_CH3Minus_C2H6_H_cm3_6_27',
        'C2H5Plus_e_C2H4_H_cm3_6_28',
    ],
}

print("=" * 100)
print("VALIDATION: Comparing define_rates.py vs define_rates_tunable.py at Te=1.0 eV")
print("=" * 100)

total_checked = 0
mismatches = []

for group_name, reactions in groups_to_check.items():
    print(f"\n{group_name}")
    print("-" * 100)
    print(f"{'Reaction':<40} {'Original':<15} {'Te=1 eV':<15} {'Ratio':<15} {'Status'}")
    print("-" * 100)

    for reaction in reactions:
        if reaction in k_original and reaction in k_tunable:
            orig = k_original[reaction]
            tune = k_tunable[reaction]
            ratio = tune / orig if orig != 0 else float('inf')

            # Check if they match within floating point tolerance
            relative_error = abs(ratio - 1.0)
            if relative_error < 1e-10:
                status = "✓ MATCH"
            else:
                status = f"✗ MISMATCH ({relative_error:.2e})"
                mismatches.append((reaction, orig, tune, ratio))

            print(f"{reaction:<40} {orig:<15.3e} {tune:<15.3e} {ratio:<15.6f} {status}")
            total_checked += 1
        else:
            print(f"{reaction:<40} {'MISSING IN ONE FILE':<60}")

print("\n" + "=" * 100)
print("SUMMARY")
print("=" * 100)
print(f"Total reactions checked: {total_checked}")
print(f"Matching rates: {total_checked - len(mismatches)}")
print(f"Mismatched rates: {len(mismatches)}")

if mismatches:
    print("\n⚠ FAILURES DETECTED:")
    print("-" * 100)
    for reaction, orig, tune, ratio in mismatches:
        print(f"  {reaction}: {orig:.3e} → {tune:.3e} (ratio: {ratio:.6f})")
    print("\nThe rates at Te=1 eV should EXACTLY match the original rates!")
    print("Check the scaling functions and ensure they return k_ref when Te=Te_ref=1.0")
else:
    print("\n✓ SUCCESS! All rates at Te=1 eV exactly match the original reference rates.")
    print("The Te-dependent implementation is correctly calibrated.")
