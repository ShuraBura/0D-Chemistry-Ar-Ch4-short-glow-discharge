#!/usr/bin/env python3
"""
Analyze how Te affects e + CH4 → CH3 rate and overall CH3 production.

User insight: e + CH4 is Te-dependent! Are we optimizing Te properly?
"""

import numpy as np
import json

# Load current result
with open('optimization_results_fixed_ionization/best_f22.8.json', 'r') as f:
    data = json.load(f)

current_Te = data['Te']
print('=' * 80)
print('Te DEPENDENCE ANALYSIS')
print('=' * 80)
print()
print(f'Current Te: {current_Te:.2f} eV')
print()

# Electron impact scaling function
def scale_electron_impact(k_ref, Te, Te_ref=1.0, E_threshold=None):
    if E_threshold is not None and E_threshold > 0:
        return k_ref * np.sqrt(Te/Te_ref) * np.exp(-E_threshold * (1/Te - 1/Te_ref))
    else:
        return k_ref * (Te/Te_ref)**0.7

# e + CH4 → CH3 + H rate
k_ref = 4.2e-11  # cm³/s at Te = 1.0 eV
E_threshold = 8.5  # eV

print('e + CH4 → CH3 + H + e (Electron-Impact Dissociation):')
print(f'  Base rate at Te=1.0 eV: {k_ref:.2e} cm³/s')
print(f'  Threshold energy: {E_threshold:.1f} eV')
print()

# Calculate rate at different Te values
Te_values = [0.5, 1.0, 1.24, 1.5, 2.0, 2.5, 3.0]
print('Rate vs Te:')
print('-' * 60)
print(f"{'Te (eV)':<10} {'k (cm³/s)':<15} {'Relative to Te=1.24':<20}")
print('-' * 60)

k_at_current = scale_electron_impact(k_ref, current_Te, E_threshold=E_threshold)

for Te in Te_values:
    k = scale_electron_impact(k_ref, Te, E_threshold=E_threshold)
    ratio = k / k_at_current
    marker = ' ← CURRENT' if abs(Te - current_Te) < 0.01 else ''
    print(f"{Te:<10.2f} {k:<15.2e} {ratio:<20.1f}× {marker}")

print()
print('=' * 80)
print('IMPACT ON CH3 PRODUCTION')
print('=' * 80)
print()

# From earlier analysis: CH3 production pathways
print('Current CH3 production breakdown:')
print('  1. Ar* + CH4 → CH3: 56.4%')
print('  2. ArPlus + CH4 → CH3: 25.5%')
print(f'  3. e + CH4 → CH3: 12.0% (k = {k_at_current:.2e} at Te={current_Te:.2f})')
print('  4. Other: 6.1%')
print()

# If we increased Te to 3.0 eV
Te_high = 3.0
k_high = scale_electron_impact(k_ref, Te_high, E_threshold=E_threshold)
boost = k_high / k_at_current

print(f'If we increased Te from {current_Te:.2f} eV to {Te_high:.1f} eV:')
print(f'  - e + CH4 → CH3 rate: {k_high:.2e} (×{boost:.0f})')
print(f'  - Currently contributes 12% of CH3 production')
print(f'  - Would contribute much more (×{boost:.0f} boost to this channel)')
print()

# Estimate total CH3 boost
# Assume other channels stay same, e + CH4 channel gets boosted
current_e_CH4_fraction = 0.12
other_fraction = 1 - current_e_CH4_fraction

# New total rate relative to current
new_total_relative = other_fraction + current_e_CH4_fraction * boost
print(f'Estimated total CH3 production boost: ×{new_total_relative:.1f}')
print(f'  (Assuming other channels unchanged)')
print()

# To get [CH3] 37× higher (what we need), we'd need production 37× higher
# Since production ∝ [CH3] at steady state
target_boost = 37
print(f'We need [CH3] ×{target_boost} to get C2H2 high enough.')
print(f'Boosting Te alone gives ×{new_total_relative:.1f} CH3 production.')
print()

if new_total_relative >= target_boost:
    print(f'✓ Te boost ALONE could achieve target!')
else:
    print(f'✗ Te boost helps but not enough (need ×{target_boost/new_total_relative:.1f} more)')
print()

print('=' * 80)
print('WHY IS Te SO LOW IN OPTIMIZATIONS?')
print('=' * 80)
print()

print(f'Optimizer bounds: Te ∈ [0.5, 3.0] eV')
print(f'Optimizer chose: Te = {current_Te:.2f} eV (near lower bound!)')
print()
print('Possible reasons:')
print('  1. Higher Te increases ionization (creates charge imbalance)')
print('  2. Higher Te increases electron-impact losses of other species')
print('  3. Objective function penalizes high Te indirectly')
print('  4. Charge balance constraint forces low Te')
print()

print('RECOMMENDATION:')
print('  - Investigate why optimizer prefers low Te')
print('  - Check if ionization scales faster than dissociation with Te')
print('  - If ionization is the problem, this explains EVERYTHING:')
print('    High Te → More e + CH4 → CH3 (GOOD!)')
print('    High Te → More e + Ar → Ar+ (BAD for charge balance)')
print('    Optimizer forced to low Te to maintain charge balance!')
