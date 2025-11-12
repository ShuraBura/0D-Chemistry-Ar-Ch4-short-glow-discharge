#!/usr/bin/env python3
"""
Analyze the boosted metastables result to understand why C2H2 remains low.

Compare Ar* levels between:
1. Fixed ionization (stable, C2H2=3.54e+09)
2. Boosted metastables (stable, C2H2=5.53e+08)
3. Unstable high-ionization (C2H2=4.09e+12, Ni/Ne=215)
"""

import json
import numpy as np
from define_rates import define_rates

# Load results
with open('optimization_results_fixed_ionization/best_f22.8.json', 'r') as f:
    fixed_ion = json.load(f)

with open('optimization_results_boosted_metastables/best_f20.9.json', 'r') as f:
    boosted_meta = json.load(f)

with open('optimization_results_comprehensive_1e12/best_f13946912.1.json', 'r') as f:
    unstable = json.load(f)

print('=' * 80)
print('COMPARISON: Fixed Ionization vs Boosted Metastables vs Unstable')
print('=' * 80)
print()

states = [
    ('FIXED IONIZATION (baseline)', fixed_ion),
    ('BOOSTED METASTABLES (attempted)', boosted_meta),
    ('UNSTABLE HIGH-IONIZATION', unstable),
]

for label, data in states:
    print(f'{label}:')
    print('-' * 80)
    print(f'  C2H2:     {data["all_densities"]["C2H2"]:.2e} cm⁻³')
    print(f'  Ar*:      {data["all_densities"]["ArStar"]:.2e} cm⁻³ (metastable)')
    print(f'  CH3:      {data["all_densities"]["CH3"]:.2e} cm⁻³')
    print(f'  Ni/Ne:    {data["Ni_over_Ne"]:.2f}', '✓' if 2 <= data["Ni_over_Ne"] <= 7 else '✗')
    print(f'  Te:       {data["Te"]:.2f} eV')
    print(f'  E-field:  {data["E_field"]:.1f} V/cm')

    # Check Ar* excitation rate
    params = {
        'P': 500.0, 'Te': data['Te'], 'ne': data['Ne'], 'E_field': data['E_field'],
        'L_discharge': 0.45,
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
        },
    }
    k = define_rates(params)

    # Check what Ar* rates are available
    ar_star_keys = [key for key in k.keys() if 'ArStar' in key or 'Ar_Ar' in key]
    print(f'  \n  Ar* excitation rate (e + Ar → Ar*):')

    found_excitation = False
    for key in ar_star_keys:
        if 'e_Ar' in key and 'ArStar' in key:
            rate = data.get('rate_values', {}).get(key, k.get(key, 'N/A'))
            print(f'    {key}: {rate:.2e}' if isinstance(rate, float) else f'    {key}: {rate}')
            found_excitation = True

    if not found_excitation:
        print('    NOT FOUND IN MODEL!')

    # Check Ar* + CH4 rates
    ar_star_ch4_key = 'ArStar_CH4_CH3_H_cm3_3_1'
    if ar_star_ch4_key in data.get('rate_values', {}):
        print(f'  \n  Ar* + CH4 → CH3 + H (Penning dissociation):')
        print(f'    k = {data["rate_values"][ar_star_ch4_key]:.2e} cm³/s')

    # Check Ar* wall loss
    ar_star_wall_key = 'stick_ArStar_9_5'
    if ar_star_wall_key in data.get('rate_values', {}):
        print(f'  \n  Ar* wall loss:')
        print(f'    k = {data["rate_values"][ar_star_wall_key]:.2e} s⁻¹')

    print()

print('=' * 80)
print('KEY FINDINGS')
print('=' * 80)
print()

# Compare Ar* levels
print('Ar* Density Comparison:')
print(f'  Fixed ionization:     {fixed_ion["all_densities"]["ArStar"]:.2e} cm⁻³')
print(f'  Boosted metastables:  {boosted_meta["all_densities"]["ArStar"]:.2e} cm⁻³')
print(f'  Ratio:                {boosted_meta["all_densities"]["ArStar"]/fixed_ion["all_densities"]["ArStar"]:.2f}×')
print()

print('CH3 Density Comparison:')
print(f'  Fixed ionization:     {fixed_ion["all_densities"]["CH3"]:.2e} cm⁻³')
print(f'  Boosted metastables:  {boosted_meta["all_densities"]["CH3"]:.2e} cm⁻³')
print(f'  Ratio:                {boosted_meta["all_densities"]["CH3"]/fixed_ion["all_densities"]["CH3"]:.2f}×')
print()

print('C2H2 Density Comparison:')
print(f'  Fixed ionization:     {fixed_ion["all_densities"]["C2H2"]:.2e} cm⁻³')
print(f'  Boosted metastables:  {boosted_meta["all_densities"]["C2H2"]:.2e} cm⁻³')
print(f'  Ratio:                {boosted_meta["all_densities"]["C2H2"]/fixed_ion["all_densities"]["C2H2"]:.2f}×')
print()

print('Comparison with Unstable State:')
print(f'  Unstable Ar*: {unstable["all_densities"]["ArStar"]:.2e} cm⁻³')
print(f'  Unstable is {unstable["all_densities"]["ArStar"]/boosted_meta["all_densities"]["ArStar"]:.0f}× higher than boosted!')
print()

print('=' * 80)
print('CONCLUSION')
print('=' * 80)
print()
print('The "boosted metastables" optimizer did NOT significantly increase Ar*!')
print()
print('Possible reasons:')
print('1. e + Ar → Ar* rate constant may not be in the model (naming issue?)')
print('2. Optimizer found that boosting Ar* ruins charge balance via Penning ionization')
print('3. Ar* + CH4 has TWO channels: dissociation (→CH3) and ionization (→CH4+)')
print('   Boosting Ar* increases BOTH, so high Ar* → high ionization → bad Ni/Ne')
print()
print('This explains why the unstable state had Ar* 1000× higher:')
print('  It was using unphysical ionization rates AND high Ar*')
print('  Both contributed to extreme charge imbalance (Ni/Ne=215)')
