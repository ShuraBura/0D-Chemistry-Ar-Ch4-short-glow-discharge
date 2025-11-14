#!/usr/bin/env python3
"""
Comprehensive review of all rate constants in the model.
Focus on finding issues like the H + C2H2 problem:
1. Reactions with missing temperature dependencies
2. Reactions using wrong temperature ranges
3. Key reactions for C2, CH, H chemistry
"""

import numpy as np
from pathlib import Path
import sys
import re

sys.path.insert(0, str(Path(__file__).parent))

from define_rates import define_rates

# Test at two temperatures to see which rates have T-dependence
T_low = 300  # K
T_high = 1500  # K

print("="*80)
print("COMPREHENSIVE RATE CONSTANT REVIEW")
print("="*80)
print()
print("Strategy: Compare rates at T=300K vs T=1500K to identify T-dependence")
print()

# Get rates at both temperatures
params_low = {
    'Tgas': T_low, 'Te': 1.5, 'ne': 2.3e9, 'P': 400, 'E_field': 150, 'L_discharge': 0.45,
    'mobilities': {'ArPlus': 3057.28}  # Minimal mobilities
}
params_high = {
    'Tgas': T_high, 'Te': 1.5, 'ne': 2.3e9, 'P': 400, 'E_field': 150, 'L_discharge': 0.45,
    'mobilities': {'ArPlus': 3057.28}
}

k_low = define_rates(params_low)
k_high = define_rates(params_high)

# Categorize reactions by species involvement
key_species = {
    'H': [],
    'CH': [],
    'C2': [],
    'C2H2': [],
    'CH2': [],
    'C': [],
}

constant_rates = []
temp_dependent_rates = []

for key in k_low.keys():
    ratio = k_high[key] / k_low[key] if k_low[key] > 0 else 1.0

    # Classify by temperature dependence
    if abs(ratio - 1.0) < 0.01:  # Less than 1% change
        constant_rates.append((key, k_low[key]))
    else:
        temp_dependent_rates.append((key, k_low[key], k_high[key], ratio))

    # Categorize by species
    for species in key_species.keys():
        if species in key.upper():
            key_species[species].append(key)

print("="*80)
print("STATISTICS")
print("="*80)
print(f"\nTotal reactions: {len(k_low)}")
print(f"Constant rates (T-independent): {len(constant_rates)}")
print(f"Temperature-dependent rates: {len(temp_dependent_rates)}")
print()

# Focus on key reactions for C2, CH, H
print("="*80)
print("KEY REACTIONS: C2 CHEMISTRY")
print("="*80)
print()

c2_reactions = [k for k in k_low.keys() if 'C2' in k and 'C2H' not in k]
print(f"Found {len(c2_reactions)} reactions involving C2")
print()

# Highlight potentially problematic ones
print("Checking for missing T-dependence in C2 reactions:")
print("-" * 80)

problematic = []

for key in c2_reactions:
    ratio = k_high[key] / k_low[key] if k_low[key] > 0 else 1.0

    # Flag reactions that should probably have T-dependence but don't
    if abs(ratio - 1.0) < 0.01:  # Constant
        # Check if it involves barrier-crossing reactions
        if any(substring in key for substring in ['_H_', '_C_', '_CH_', '_O_']):
            print(f"⚠  {key}")
            print(f"   k = {k_low[key]:.2e} cm³/s (constant)")
            problematic.append(key)

print()
print(f"Found {len(problematic)} potentially problematic C2 reactions")
print()

# Check C2 + H specifically
print("="*80)
print("CRITICAL REACTION: C2 + H → CH + C")
print("="*80)
print()

c2_h_keys = [k for k in k_low.keys() if 'C2' in k and ('_H_' in k or '_H_CH' in k) and 'C2H' not in k]
for key in c2_h_keys:
    if 'CH_C' in key or 'C_H2' in key:
        ratio = k_high[key] / k_low[key] if k_low[key] > 0 else 1.0
        print(f"Reaction: {key}")
        print(f"  k(300K)  = {k_low[key]:.2e} cm³/s")
        print(f"  k(1500K) = {k_high[key]:.2e} cm³/s")
        print(f"  Ratio: {ratio:.3f}")
        if abs(ratio - 1.0) < 0.01:
            print(f"  ⚠ WARNING: Rate is constant (should it be?)")
        print()

# Check CH reactions
print("="*80)
print("KEY REACTIONS: CH CHEMISTRY")
print("="*80)
print()

ch_reactions = [k for k in k_low.keys() if ('CH_' in k or '_CH_' in k) and 'CH2' not in k and 'CH3' not in k and 'CH4' not in k and 'CH5' not in k]
print(f"Found {len(ch_reactions)} reactions involving CH")
print()

# Focus on CH + H (destroys CH)
print("CH + H reactions:")
print("-" * 80)
ch_h_keys = [k for k in ch_reactions if '_H_' in k or 'CH_H' in k or 'H_CH' in k]
for key in ch_h_keys[:10]:  # Show first 10
    ratio = k_high[key] / k_low[key] if k_low[key] > 0 else 1.0
    status = "constant" if abs(ratio - 1.0) < 0.01 else f"T-dep (×{ratio:.2f})"
    print(f"{key:50s} k(300K)={k_low[key]:.2e}  [{status}]")

print()

# Check H + C2H2 (we know this one is wrong)
print("="*80)
print("CONFIRMED PROBLEMATIC: H + C2H2 → C2")
print("="*80)
print()

c2h2_h_keys = [k for k in k_low.keys() if 'C2H2' in k and '_H_' in k and 'C2_H2' in k]
for key in c2h2_h_keys:
    ratio = k_high[key] / k_low[key] if k_low[key] > 0 else 1.0
    print(f"Reaction: {key}")
    print(f"  k(300K)  = {k_low[key]:.2e} cm³/s")
    print(f"  k(1500K) = {k_high[key]:.2e} cm³/s")
    print(f"  Ratio: {ratio:.3f}")

    # Calculate what it should be with Baulch
    k_baulch_300 = 1.67e-14 * (300**1.64) * np.exp(-15250/300)
    k_baulch_1500 = 1.67e-14 * (1500**1.64) * np.exp(-15250/1500)
    print(f"\n  Baulch et al. (2005):")
    print(f"  k(300K)  = {k_baulch_300:.2e} cm³/s (should be)")
    print(f"  k(1500K) = {k_baulch_1500:.2e} cm³/s (should be)")
    print(f"  Baulch ratio: {k_baulch_1500/k_baulch_300:.2e}")
    print()

# Check other potentially problematic reactions
print("="*80)
print("OTHER REACTIONS TO CHECK")
print("="*80)
print()

print("1. Radical-radical recombination (should be T-independent or weak T-dep):")
print("-" * 80)
recomb_keys = ['CH_CH2_C2H2_H_cm3', 'CH2_CH2_C2H2_H2_cm3', 'CH_CH_C2H2_cm3']
for pattern in recomb_keys:
    matches = [k for k in k_low.keys() if pattern in k]
    for key in matches[:3]:
        ratio = k_high[key] / k_low[key] if k_low[key] > 0 else 1.0
        status = "✓ constant" if abs(ratio - 1.0) < 0.01 else f"⚠ T-dep (×{ratio:.2f})"
        print(f"{key:50s} [{status}]")

print()
print("2. Atom abstraction reactions (should have activation energy):")
print("-" * 80)
abstraction_patterns = ['H_CH4', 'CH_CH4', 'C_CH4']
for pattern in abstraction_patterns:
    matches = [k for k in k_low.keys() if pattern in k]
    for key in matches[:3]:
        ratio = k_high[key] / k_low[key] if k_low[key] > 0 else 1.0
        if ratio > 1.5:
            status = f"✓ T-dep (×{ratio:.1f})"
        elif abs(ratio - 1.0) < 0.01:
            status = "⚠ constant (should have barrier?)"
        else:
            status = f"? weak T-dep (×{ratio:.2f})"
        print(f"{key:50s} [{status}]")

print()
print("="*80)
print("SUMMARY OF CONCERNS")
print("="*80)
print()
print("Reactions that may need review:")
print()
print("1. H + C2H2 → C2 + H2 + H")
print("   - CONFIRMED WRONG: Using constant k=1e-11, should be Baulch formula")
print("   - Already identified and corrected")
print()

print("2. C2 + H → CH + C")
print("   - Check if rate is energy-dependent (for hot C2)")
print("   - May need vibrational state-specific rate")
print()

print("3. Other C2 production reactions:")
for key in ['CH_CH2_C2', 'CH2_CH2_C2', 'C2H2_C_C2']:
    matches = [k for k in k_low.keys() if key in k]
    if matches:
        print(f"   - {matches[0]}: k = {k_low[matches[0]]:.2e} (verify with literature)")

print()
print("NEXT: Conduct literature review to verify these rates")
print("="*80)
