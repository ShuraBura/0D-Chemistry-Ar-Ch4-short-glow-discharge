"""
Quick test to verify three-body electron-ion recombination reactions were added correctly
"""

import numpy as np
from define_rates import define_rates

# Minimal params for rate definition
params = {
    'P': 500,
    'Te': 1.3,
    'ne': 1e8,
    'E_field': 250,
    'L_discharge': 0.45,
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
    },
}

k = define_rates(params)

print("Three-body electron-ion recombination rates added:")
print("="*70)

three_body_e_ion = [
    'e_ArPlus_M_Ar_M_cm6_8_4',
    'e_CH4Plus_M_CH4_M_cm6_8_5',
    'e_CH3Plus_M_CH3_M_cm6_8_6',
    'e_CH5Plus_M_CH5_M_cm6_8_7',
    'e_ArHPlus_M_ArH_M_cm6_8_8',
    'e_C2H5Plus_M_C2H5_M_cm6_8_9',
]

for rate_name in three_body_e_ion:
    if rate_name in k:
        print(f"✓ {rate_name:50s} = {k[rate_name]:.2e} cm⁶/s")
    else:
        print(f"✗ {rate_name} MISSING!")

print("\nCompare to existing three-body reactions:")
print("-"*70)
existing_three_body = [
    'H_H_M_H2_M_cm6_8_1',
    'CH3_CH3_M_C2H6_M_cm6_8_2',
    'CH3_H_M_CH4_M_cm6_8_3',
]

for rate_name in existing_three_body:
    if rate_name in k:
        print(f"  {rate_name:50s} = {k[rate_name]:.2e} cm⁶/s")

print("\n" + "="*70)
print("Verification complete!")

# Calculate impact at different densities
print("\nExpected impact at different electron densities:")
print("-"*70)

n_total = 1.2e16  # cm⁻³ at 500 mTorr
k_3body = 1e-25  # cm⁶/s
k_2body = 1.5e-7  # cm³/s (dissociative recombination)

for ne, ni in [(1e8, 1e8), (1e9, 1e9), (1e10, 1e10)]:
    rate_3body = k_3body * ne * ni * n_total
    rate_2body = k_2body * ne * ni
    ratio = rate_3body / rate_2body

    print(f"ne = ni = {ne:.1e}:")
    print(f"  Two-body rate:   {rate_2body:.2e} cm⁻³/s")
    print(f"  Three-body rate: {rate_3body:.2e} cm⁻³/s")
    print(f"  Ratio (3-body/2-body): {ratio:.2e}")
    print()
