#!/usr/bin/env python3
"""
Verify all rate constants against literature ranges from rate_database_complete.py
"""

import numpy as np
from define_rates import define_rates
from rate_database_complete import get_complete_rate_database

# Test parameters at different temperatures
test_conditions = [
    {'Te': 0.5, 'label': 'Te=0.5 eV (low)'},
    {'Te': 1.0, 'label': 'Te=1.0 eV (reference)'},
    {'Te': 2.0, 'label': 'Te=2.0 eV (moderate)'},
    {'Te': 3.0, 'label': 'Te=3.0 eV (high)'},
]

params_base = {
    'E_field': 50,
    'L_discharge': 0.45,
    'Tgas': 400,
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
    }
}

# Get literature database
db = get_complete_rate_database()

print("="*80)
print(" RATE CONSTANT VERIFICATION AGAINST LITERATURE")
print("="*80)
print()

# Track statistics
total_rates_checked = 0
out_of_range_count = {}
for cond in test_conditions:
    out_of_range_count[cond['label']] = 0

# Check each temperature condition
for condition in test_conditions:
    print(f"\n{'='*80}")
    print(f" Checking {condition['label']}")
    print(f"{'='*80}")

    params = params_base.copy()
    params['Te'] = condition['Te']
    k = define_rates(params)

    out_of_range = []
    warnings = []

    for rate_name in sorted(k.keys()):
        if rate_name in db:
            total_rates_checked += 1
            rate_val = k[rate_name]
            rate_info = db[rate_name]

            # Check if within range
            if not rate_info.is_within_range(rate_val):
                out_of_range.append({
                    'name': rate_name,
                    'value': rate_val,
                    'min': rate_info.min,
                    'max': rate_info.max,
                    'source': rate_info.source
                })
                out_of_range_count[condition['label']] += 1

    # Report findings for this condition
    if len(out_of_range) == 0:
        print(f"\n✓ ALL rates within literature ranges!")
    else:
        print(f"\n✗ Found {len(out_of_range)} rates outside literature range:")
        print(f"\n{'Rate Name':<40} {'Value':>12} {'Min':>12} {'Max':>12} {'Source'}")
        print("-"*80)
        for r in out_of_range:
            print(f"{r['name']:<40} {r['value']:>12.2e} {r['min']:>12.2e} {r['max']:>12.2e} {r['source']}")

# Summary
print(f"\n{'='*80}")
print(" SUMMARY")
print(f"{'='*80}")
print(f"\nTotal rates in define_rates: {len(k)}")
print(f"Total rates checked against database: {total_rates_checked}")
print(f"Rates with literature data: {len(db)}")
print()

for cond in test_conditions:
    status = "✓" if out_of_range_count[cond['label']] == 0 else "✗"
    print(f"{status} {cond['label']:30} {out_of_range_count[cond['label']]:3} out of range")

# Check rates NOT in database
print(f"\n{'='*80}")
print(" RATES WITHOUT LITERATURE DATA")
print(f"{'='*80}")
print()

rates_without_lit = []
for rate_name in sorted(k.keys()):
    if rate_name not in db:
        rates_without_lit.append(rate_name)

if len(rates_without_lit) > 0:
    print(f"Found {len(rates_without_lit)} rates without literature ranges:")
    for i, rate_name in enumerate(rates_without_lit, 1):
        print(f"  {i:3}. {rate_name}")
else:
    print("✓ All rates have literature data!")

# Temperature dependence verification
print(f"\n{'='*80}")
print(" TEMPERATURE DEPENDENCE VERIFICATION")
print(f"{'='*80}")
print()

# Check that electron-impact rates increase with Te
params_low = params_base.copy()
params_low['Te'] = 1.0
k_low = define_rates(params_low)

params_high = params_base.copy()
params_high['Te'] = 2.0
k_high = define_rates(params_high)

print("Electron-impact dissociation (should INCREASE with Te):")
e_impact_tests = [
    'e_CH4_CH3_H_cm3_1_1',
    'e_H2_H_H_cm3_1_4',
    'e_C2H4_C2H2_H2_cm3_1_6'
]
for rate_name in e_impact_tests:
    ratio = k_high[rate_name] / k_low[rate_name]
    status = "✓" if ratio > 1.0 else "✗"
    print(f"  {status} {rate_name:35} {ratio:6.2f}× increase")

print("\nIonization (should INCREASE with Te):")
ion_tests = [
    'e_CH4_CH3Plus_H_cm3_2_1',
    'e_Ar_ArPlus_cm3_2_3',
    'e_H2_H2Plus_2e_cm3_2_9'
]
for rate_name in ion_tests:
    ratio = k_high[rate_name] / k_low[rate_name]
    status = "✓" if ratio > 1.0 else "✗"
    print(f"  {status} {rate_name:35} {ratio:6.2f}× increase")

print("\nRecombination (should DECREASE with Te):")
rec_tests = [
    'ArPlus_e_Ar_cm3_6_1',
    'CH3Plus_e_CH3_cm3_6_2',
    'H3Plus_e_H2_H_cm3_6_30'
]
for rate_name in rec_tests:
    ratio = k_high[rate_name] / k_low[rate_name]
    status = "✓" if ratio < 1.0 else "✗"
    print(f"  {status} {rate_name:35} {ratio:6.2f}× decrease")

print("\nNeutral-neutral (should be CONSTANT):")
neutral_tests = [
    'CH_CH3_C2H4_cm3_7_5',
    'CH_H_C_H2_cm3_7_3',
    'CH3_CH3_C2H6_cm3_7_40'
]
for rate_name in neutral_tests:
    diff = abs(k_high[rate_name] - k_low[rate_name]) / k_low[rate_name]
    status = "✓" if diff < 1e-10 else "✗"
    print(f"  {status} {rate_name:35} {diff:6.2e} relative diff")

print(f"\n{'='*80}")
print(" VERIFICATION COMPLETE")
print(f"{'='*80}")
