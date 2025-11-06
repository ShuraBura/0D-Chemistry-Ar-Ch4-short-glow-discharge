#!/usr/bin/env python3
"""
Test that stick and loss coefficients are properly overridden from rate_values
"""

import json
from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
import numpy as np

# Load checkpoint
with open('checkpoint_f3407.json', 'r') as f:
    data = json.load(f)

params = data['params']

print("="*80)
print("TESTING RATE OVERRIDE MECHANISM")
print("="*80)
print()

# Get baseline rates (without overrides)
k_baseline = define_rates_tunable(params)

print("BASELINE VALUES (from define_rates_tunable hardcoded values):")
print(f"  stick_C2_9_9:        {k_baseline['stick_C2_9_9']:.2e}")
print(f"  stick_CH_9_3:        {k_baseline['stick_CH_9_3']:.2e}")
print(f"  loss_C2_11_3:        {k_baseline['loss_C2_11_3']:.2e}")
print(f"  stick_C2H2_9_11:     {k_baseline['stick_C2H2_9_11']:.2e}")
print(f"  C2_H_CH_C_cm3_7_6:   {k_baseline['C2_H_CH_C_cm3_7_6']:.2e}")
print()

# Apply rate_values overrides (as optimizer does)
k_override = define_rates_tunable(params)
db = get_complete_rate_database()

applied_count = 0
for name, val in params.get('rate_values', {}).items():
    if name in k_override and name in db:
        k_override[name] = np.clip(val, db[name].min, db[name].max)
        applied_count += 1

print(f"Applied {applied_count} overrides from rate_values")
print()

print("VALUES AFTER RATE_VALUES OVERRIDE:")
print(f"  stick_C2_9_9:        {k_override['stick_C2_9_9']:.2e} (from rate_values: {params['rate_values']['stick_C2_9_9']:.2e})")
print(f"  stick_CH_9_3:        {k_override['stick_CH_9_3']:.2e} (from rate_values: {params['rate_values']['stick_CH_9_3']:.2e})")
print(f"  loss_C2_11_3:        {k_override['loss_C2_11_3']:.2e} (from rate_values: {params['rate_values']['loss_C2_11_3']:.2e})")
print(f"  stick_C2H2_9_11:     {k_override['stick_C2H2_9_11']:.2e} (from rate_values: {params['rate_values']['stick_C2H2_9_11']:.2e})")
print(f"  C2_H_CH_C_cm3_7_6:   {k_override['C2_H_CH_C_cm3_7_6']:.2e} (NOT in rate_values - using baseline)")
print()

# Check if they changed
print("VERIFICATION:")
changed = []
unchanged = []

test_params = [
    ('stick_C2_9_9', 'stick_C2_9_9'),
    ('stick_CH_9_3', 'stick_CH_9_3'),
    ('loss_C2_11_3', 'loss_C2_11_3'),
    ('stick_C2H2_9_11', 'stick_C2H2_9_11'),
]

for display_name, key in test_params:
    baseline_val = k_baseline[key]
    override_val = k_override[key]

    if abs(baseline_val - override_val) > 1e-6:
        changed.append(display_name)
        print(f"  ✓ {display_name}: CHANGED from {baseline_val:.2e} to {override_val:.2e}")
    else:
        unchanged.append(display_name)
        print(f"  ✗ {display_name}: UNCHANGED at {baseline_val:.2e}")

print()
print("="*80)
print(f"RESULT: {len(changed)}/{len(test_params)} parameters properly overridden")
print("="*80)

if len(changed) == len(test_params):
    print("✓ ALL TESTED PARAMETERS ARE BEING APPLIED CORRECTLY")
else:
    print(f"✗ {len(unchanged)} PARAMETERS NOT BEING APPLIED: {unchanged}")
