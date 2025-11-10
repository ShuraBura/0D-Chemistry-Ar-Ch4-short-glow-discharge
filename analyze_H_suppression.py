#!/usr/bin/env python3
"""
Analyze why the optimizer suppresses H even with equal weights.
Compare optimized rates to literature values and identify what's being tuned.
"""

import numpy as np
import json
from define_rates import define_rates

print("=" * 80)
print("H SUPPRESSION ANALYSIS")
print("=" * 80)
print()

# Load best result with equal weights
with open('optimization_results_charge_balanced/best_f41.3.json', 'r') as f:
    result = json.load(f)

print("Result summary:")
print(f"  H = {result['all_densities']['H']:.2e} cm⁻³ (target: 2.52e14, ratio: {result['all_densities']['H']/2.52e14:.3f})")
print(f"  CH = {result['all_densities']['CH']:.2e} cm⁻³ (target: 1.0e9, ratio: {result['all_densities']['CH']/1e9:.2f})")
print(f"  C2 = {result['all_densities']['C2']:.2e} cm⁻³ (target: 5.6e11, ratio: {result['all_densities']['C2']/5.6e11:.3f})")
print(f"  Te = {result['Te']:.3f} eV")
print(f"  Ne = {result['Ne']:.2e} cm⁻³")
print()

# Get default rates for comparison
mobilities = {
    'ArPlus': 1.54e3, 'CH4Plus': 1.54e3, 'CH3Plus': 1.54e3,
    'CH5Plus': 1.54e3, 'ArHPlus': 1.54e3, 'H3Plus': 1.54e3,
    'CH2Plus': 1.54e3, 'C2H5Plus': 1.54e3, 'C2H4Plus': 1.54e3,
    'C2H3Plus': 1.54e3, 'C2HPlus': 1.54e3, 'H2Plus': 1.54e3,
    'CHPlus': 1.54e3, 'CH3Minus': 1.54e3, 'HMinus': 1.54e3
}

params = {
    'Te': result['Te'],
    'ne': result['Ne'],
    'E_field': result['E_field'],
    'pressure': 500.0,
    'T': 400.0,
    'Tgas': 400.0,
    'L_discharge': 0.45,
    'mobilities': mobilities
}

default_rates = define_rates(params)

print("=" * 80)
print("OPTIMIZED VS DEFAULT RATES")
print("=" * 80)
print()

# Compare each optimized rate to its default
rate_comparisons = []

for name, optimized_val in result['rate_values'].items():
    if name.startswith('stick_'):
        # Sticking coefficient - should be <= 1.0
        species = name.split('_')[1]
        print(f"{name}:")
        print(f"  Optimized: {optimized_val:.1f} s⁻¹")

        # These are loss rates, not coefficients
        # For 500 mTorr, 400K, L=0.45cm discharge
        # Typical rate would be v̄/L where v̄ ~ 1000 m/s
        # k_loss ~ 100,000 / L(cm) ~ 220,000 s⁻¹ for unit sticking
        typical_max = 220000.0

        if optimized_val > typical_max:
            print(f"  → UNPHYSICAL! ({optimized_val/typical_max:.1f}× max possible)")
        elif optimized_val < 10:
            print(f"  → Nearly zero sticking")
        else:
            coeff = optimized_val / typical_max
            print(f"  → Sticking coefficient ~ {coeff:.4f}")

        if species == 'H':
            print(f"  *** H STICKING - KEY PARAMETER ***")

    elif name.startswith('loss_'):
        # Other loss mechanisms
        species = name.split('_')[1]
        print(f"{name}:")
        print(f"  Optimized: {optimized_val:.1f} s⁻¹")

    else:
        # Chemical reaction rate
        if name in default_rates:
            default_val = default_rates[name]
            ratio = optimized_val / default_val

            rate_comparisons.append({
                'name': name,
                'default': default_val,
                'optimized': optimized_val,
                'ratio': ratio
            })

            if abs(ratio - 1.0) > 0.5:  # Changed by more than 50%
                print(f"{name}:")
                print(f"  Default:   {default_val:.2e} cm³/s")
                print(f"  Optimized: {optimized_val:.2e} cm³/s")
                print(f"  Ratio: {ratio:.2f}×")

                if 'CH4' in name and 'dissociation' not in name:
                    print(f"  → Affects CH4 chemistry")
                if '_CH_' in name or name.startswith('e_CH_'):
                    print(f"  *** AFFECTS CH ***")
                if 'H_' in name or '_H_' in name:
                    print(f"  *** AFFECTS H ***")
        else:
            print(f"{name}: NOT IN DEFAULT RATES (optimized: {optimized_val:.2e})")

print()
print("=" * 80)
print("KEY FINDINGS")
print("=" * 80)
print()

# Find the most changed rates
rate_comparisons.sort(key=lambda x: abs(np.log10(x['ratio'])), reverse=True)

print("Most tuned rates:")
for rc in rate_comparisons[:10]:
    print(f"  {rc['name']}: {rc['ratio']:.2f}× (def: {rc['default']:.2e}, opt: {rc['optimized']:.2e})")

print()

# Check H sticking specifically
if 'stick_H_9_1' in result['rate_values']:
    stick_H = result['rate_values']['stick_H_9_1']
    typical_max = 220000.0
    coeff = stick_H / typical_max

    print(f"H sticking analysis:")
    print(f"  Optimized k_wall(H) = {stick_H:.1f} s⁻¹")
    print(f"  Effective sticking coefficient γ ≈ {coeff:.4f}")
    print(f"  H lifetime = {1/stick_H*1000:.2f} ms")
    print()

    if coeff < 0.01:
        print("  → VERY LOW sticking (γ < 0.01)")
        print("  → Optimizer is minimizing H wall loss")
        print("  → But H is still low! There must be strong chemical sinks.")
    elif coeff > 0.5:
        print("  → HIGH sticking (γ > 0.5)")
        print("  → Optimizer is maximizing H wall loss to suppress H")

print()
print("=" * 80)
print("HYPOTHESIS")
print("=" * 80)
print()
print("The optimizer faces a trade-off:")
print("  • Higher H improves H target")
print("  • But higher H drives CH2+H→CH+H2, overproducing CH")
print("  • Optimizer sacrifices H to keep CH under control")
print()
print("Even with equal weights (H=20, CH=20), the constraint is:")
print("  • CH is very sensitive to H (via CH2+H reaction)")
print("  • Small increase in H causes large increase in CH")
print("  • Optimizer cannot increase H without breaking CH")
