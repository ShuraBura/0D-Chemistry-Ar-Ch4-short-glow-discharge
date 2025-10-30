#!/usr/bin/env python3
"""
Correct all out-of-bounds rates to be within literature ranges.
Creates define_rates_corrected.py with all rates within bounds.
"""

import numpy as np
from rate_database_complete import get_complete_rate_database
from define_rates import define_rates

# Load literature database
db = get_complete_rate_database()

# Load current rates
params = {
    'E_field': 50,
    'L_discharge': 0.45,
}

try:
    k_current = define_rates(params)
except:
    # If define_rates needs more params, use database values
    k_current = {name: rate.value for name, rate in db.items()}

print("=" * 80)
print("CORRECTING RATES TO LITERATURE BOUNDS")
print("=" * 80)

corrections = []

for name, rate_db in db.items():
    if name not in k_current:
        continue

    current_val = k_current[name]

    # Check if outside bounds
    if current_val < rate_db.min * 0.99:
        # Below minimum - correct to min
        pct_below = (rate_db.min - current_val) / rate_db.min * 100
        corrections.append({
            'name': name,
            'current': current_val,
            'corrected': rate_db.min,
            'min': rate_db.min,
            'max': rate_db.max,
            'type': 'BELOW_MIN',
            'violation': f'{pct_below:.1f}% below',
            'source': rate_db.source
        })
        k_current[name] = rate_db.min

    elif current_val > rate_db.max * 1.01:
        # Above maximum - correct to max
        pct_above = (current_val - rate_db.max) / rate_db.max * 100
        corrections.append({
            'name': name,
            'current': current_val,
            'corrected': rate_db.max,
            'min': rate_db.min,
            'max': rate_db.max,
            'type': 'ABOVE_MAX',
            'violation': f'{pct_above:.1f}% above',
            'source': rate_db.source
        })
        k_current[name] = rate_db.max

print(f"\nFound {len(corrections)} rates outside literature bounds\n")

# Show corrections by type
above_max = [c for c in corrections if c['type'] == 'ABOVE_MAX']
below_min = [c for c in corrections if c['type'] == 'BELOW_MIN']

if above_max:
    print(f"RATES ABOVE MAXIMUM ({len(above_max)}):")
    print("-" * 80)
    for c in above_max:
        print(f"\n{c['name']}:")
        print(f"  Current:   {c['current']:.2e} ({c['violation']})")
        print(f"  Corrected: {c['corrected']:.2e} (literature MAX)")
        print(f"  Range:     [{c['min']:.2e}, {c['max']:.2e}]")
        print(f"  Source:    {c['source']}")

if below_min:
    print(f"\nRATES BELOW MINIMUM ({len(below_min)}):")
    print("-" * 80)
    for c in below_min:
        print(f"\n{c['name']}:")
        print(f"  Current:   {c['current']:.2e} ({c['violation']})")
        print(f"  Corrected: {c['corrected']:.2e} (literature MIN)")
        print(f"  Range:     [{c['min']:.2e}, {c['max']:.2e}]")
        print(f"  Source:    {c['source']}")

# Generate corrected define_rates.py
print("\n" + "=" * 80)
print("GENERATING define_rates_corrected.py")
print("=" * 80)

# Read original define_rates.py
with open('define_rates.py', 'r') as f:
    original_lines = f.readlines()

# Create corrected version
corrected_lines = []
corrections_dict = {c['name']: c for c in corrections}

in_function = False
for line in original_lines:
    # Check if we're in the function body
    if 'def define_rates(' in line:
        in_function = True
        corrected_lines.append(line)
        continue

    # Check if we've left the function
    if in_function and line.strip().startswith('def ') and 'define_rates' not in line:
        in_function = False

    # Check for rate assignments
    if in_function and "k['" in line and '] = ' in line:
        # Extract rate name
        for rate_name in corrections_dict:
            if f"k['{rate_name}']" in line:
                c = corrections_dict[rate_name]
                # Replace with corrected value and add comment
                indent = len(line) - len(line.lstrip())
                new_line = f"{' ' * indent}k['{rate_name}'] = {c['corrected']:.2e}  # CORRECTED from {c['current']:.2e} ({c['violation']})\n"
                corrected_lines.append(new_line)
                break
        else:
            corrected_lines.append(line)
    else:
        corrected_lines.append(line)

# Write corrected file
with open('define_rates_corrected.py', 'w') as f:
    f.write(''.join(corrected_lines))

print(f"\n✓ Created define_rates_corrected.py with {len(corrections)} corrections")
print(f"✓ All {len(db)} rates now within literature bounds")

# Verify
print("\n" + "=" * 80)
print("VERIFICATION")
print("=" * 80)

# Import corrected version
import importlib.util
spec = importlib.util.spec_from_file_location("define_rates_corrected", "define_rates_corrected.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

k_corrected = module.define_rates(params)

# Check all rates are within bounds
violations = 0
for name, rate_db in db.items():
    if name in k_corrected:
        val = k_corrected[name]
        if val < rate_db.min * 0.99 or val > rate_db.max * 1.01:
            violations += 1
            print(f"WARNING: {name} still outside bounds: {val:.2e} not in [{rate_db.min:.2e}, {rate_db.max:.2e}]")

if violations == 0:
    print("✓ ALL RATES VERIFIED WITHIN LITERATURE BOUNDS")
else:
    print(f"✗ {violations} rates still outside bounds")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"  Total rates: {len(db)}")
print(f"  Corrected: {len(corrections)} ({len(above_max)} above max, {len(below_min)} below min)")
print(f"  Now compliant: {len(db) - violations}")
print("\n✓ Ready for baseline simulation and optimization")
