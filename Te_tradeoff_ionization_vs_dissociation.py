#!/usr/bin/env python3
"""
THE KEY INSIGHT: Why can't we use high Te to boost CH3?

Compare how ionization (e + Ar → Ar+) and dissociation (e + CH4 → CH3)
scale with Te.
"""

import numpy as np

def scale_rate(k_ref, Te, Te_ref=1.0, E_threshold=None):
    if E_threshold is not None and E_threshold > 0:
        return k_ref * np.sqrt(Te/Te_ref) * np.exp(-E_threshold * (1/Te - 1/Te_ref))
    return k_ref

# Base rates at Te = 1.0 eV
k_dissoc_ref = 4.2e-11  # e + CH4 → CH3
k_ion_ref = 8e-12       # e + Ar → Ar+

# Thresholds
E_dissoc = 8.5   # eV
E_ion = 15.76    # eV

print('=' * 80)
print('THE TE TRADEOFF: Dissociation vs Ionization')
print('=' * 80)
print()

print('e + CH4 → CH3 + H (dissociation):')
print(f'  E_threshold = {E_dissoc:.1f} eV')
print(f'  k_ref = {k_dissoc_ref:.2e} at Te=1.0 eV')
print()

print('e + Ar → Ar+ + 2e (ionization):')
print(f'  E_ion = {E_ion:.2f} eV')
print(f'  k_ref = {k_ion_ref:.2e} at Te=1.0 eV')
print()

Te_current = 1.24
k_dissoc_current = scale_rate(k_dissoc_ref, Te_current, E_threshold=E_dissoc)
k_ion_current = scale_rate(k_ion_ref, Te_current, E_threshold=E_ion)

print(f'At current Te = {Te_current:.2f} eV:')
print(f'  k(dissoc) = {k_dissoc_current:.2e}')
print(f'  k(ion) = {k_ion_current:.2e}')
print()

print('=' * 80)
print('SCALING WITH Te')
print('=' * 80)
print()

Te_values = [0.5, 1.0, 1.24, 1.5, 2.0, 2.5, 3.0]

print(f"{'Te':<6} {'k(dissoc)':<12} {'k(ion)':<12} {'Dissoc boost':<14} {'Ion boost':<14} {'Ion/Dissoc':<12}")
print('-' * 80)

for Te in Te_values:
    k_d = scale_rate(k_dissoc_ref, Te, E_threshold=E_dissoc)
    k_i = scale_rate(k_ion_ref, Te, E_threshold=E_ion)

    boost_d = k_d / k_dissoc_current
    boost_i = k_i / k_ion_current
    ratio = (k_i / k_d) / (k_ion_current / k_dissoc_current)  # Normalized

    marker = ' ← CURRENT' if abs(Te - Te_current) < 0.01 else ''
    print(f"{Te:<6.2f} {k_d:<12.2e} {k_i:<12.2e} {boost_d:<14.1f}× {boost_i:<14.1f}× {ratio:<12.2f}× {marker}")

print()
print('=' * 80)
print('THE PROBLEM')
print('=' * 80)
print()

Te_high = 3.0
k_d_high = scale_rate(k_dissoc_ref, Te_high, E_threshold=E_dissoc)
k_i_high = scale_rate(k_ion_ref, Te_high, E_threshold=E_ion)

boost_d = k_d_high / k_dissoc_current
boost_i = k_i_high / k_ion_current

print(f'If we increase Te from {Te_current:.2f} eV to {Te_high:.1f} eV:')
print(f'  - Dissociation (e + CH4 → CH3): ×{boost_d:.0f}')
print(f'  - Ionization (e + Ar → Ar+): ×{boost_i:.0f}')
print()

if boost_i > boost_d:
    print(f'✗ Ionization grows FASTER than dissociation!')
    print(f'   Ionization/Dissociation ratio increases by {boost_i/boost_d:.1f}×')
    print()
    print('This means:')
    print('  - High Te → More CH3 (GOOD!)')
    print('  - High Te → Even MORE ions (BAD for charge balance!)')
    print('  - Optimizer forced to choose low Te to avoid charge imbalance')
    print()
    print('This is THE FUNDAMENTAL CONSTRAINT!')
    print('  - Need high Te for CH3 production')
    print('  - But high Te ruins charge balance via ionization')
    print('  - Cannot have both!')

print()
print('=' * 80)
print('VERIFICATION')
print('=' * 80)
print()

print('From previous analysis:')
print('  - High Ar* → charge imbalance (Ar* + CH4 → CH4+ + e)')
print('  - High Te → charge imbalance (e + Ar → Ar+ + 2e)')
print()
print('BOTH pathways to high CH3 are blocked by charge balance!')
print()
print('The experiment must be using something else:')
print('  1. Pulsed discharge (time-averaged Te is low, peak Te is high)?')
print('  2. Spatial separation (high Te in core, CH3 accumulates in edge)?')
print('  3. Different mechanism entirely?')
