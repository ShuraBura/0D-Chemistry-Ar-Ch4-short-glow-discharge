#!/usr/bin/env python3
"""
Check which rates in define_rates.py are outside literature bounds
"""

from rate_database_complete import get_complete_rate_database
from define_rates import define_rates

# Get database with literature bounds
db = get_complete_rate_database()

# Get current rates
params = {
    'pressure': 400e-3,  # Torr
    'T': 400,            # K
    'L': 0.05,           # m
    'R': 0.05,           # m
    'E_field': 100,      # V/cm (approximate)
}

try:
    k = define_rates(params)
except Exception as e:
    print(f"Error loading rates: {e}")
    print("Using database values directly...")
    k = {name: rate.value for name, rate in db.items()}

print("=" * 80)
print("RATES OUTSIDE LITERATURE BOUNDS")
print("=" * 80)

outside_bounds = []

for name, rate_db in db.items():
    if name in k:
        current_val = k[name]

        # Check if outside bounds (with 1% tolerance for floating point)
        if current_val < rate_db.min * 0.99:
            pct_below = (rate_db.min - current_val) / rate_db.min * 100
            outside_bounds.append({
                'name': name,
                'current': current_val,
                'min': rate_db.min,
                'max': rate_db.max,
                'violation': f'{pct_below:.1f}% BELOW minimum',
                'source': rate_db.source
            })
        elif current_val > rate_db.max * 1.01:
            pct_above = (current_val - rate_db.max) / rate_db.max * 100
            outside_bounds.append({
                'name': name,
                'current': current_val,
                'min': rate_db.min,
                'max': rate_db.max,
                'violation': f'{pct_above:.1f}% ABOVE maximum',
                'source': rate_db.source
            })

if outside_bounds:
    print(f"\nFound {len(outside_bounds)} rates outside literature bounds:\n")
    for item in outside_bounds:
        print(f"{item['name']}:")
        print(f"  Current: {item['current']:.2e}")
        print(f"  Literature range: [{item['min']:.2e}, {item['max']:.2e}]")
        print(f"  VIOLATION: {item['violation']}")
        print(f"  Source: {item['source']}")
        print()
else:
    print("\n✓ All rates are within literature bounds!")

print("=" * 80)
print("RATES AT BOUNDARIES (within 5%)")
print("=" * 80)

at_boundaries = []

for name, rate_db in db.items():
    if name in k:
        current_val = k[name]

        # Check if at min or max (within 5%)
        if abs(current_val - rate_db.min) / rate_db.min < 0.05:
            at_boundaries.append((name, 'MIN', current_val, rate_db.min, rate_db.max, rate_db.source))
        elif abs(current_val - rate_db.max) / rate_db.max < 0.05:
            at_boundaries.append((name, 'MAX', current_val, rate_db.min, rate_db.max, rate_db.source))

# Filter for CH-related rates
ch_at_boundaries = [x for x in at_boundaries if 'CH' in x[0]]

print(f"\nFound {len(ch_at_boundaries)} CH-related rates at boundaries:\n")
for name, boundary, current, min_val, max_val, source in ch_at_boundaries:
    range_factor = max_val / min_val if min_val > 0 else 1
    print(f"{name}: at {boundary}")
    print(f"  Current: {current:.2e}, Range: [{min_val:.2e}, {max_val:.2e}] ({range_factor:.1f}x)")
    print(f"  Source: {source}")
    print()
