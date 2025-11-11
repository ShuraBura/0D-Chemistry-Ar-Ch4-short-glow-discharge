#!/usr/bin/env python3
"""
Analyze why high E-field doesn't fully compensate for high Te ionization.

Ion production rate: R_ion ∝ ne × n_Ar × k_ion(Te)
                     where k_ion ~ exp(-15.76/Te)

Ion loss rate: R_loss ∝ n_ion × k_drift
               where k_drift = mobility × E / L

At steady state: R_ion = R_loss
So: ne × n_Ar × k_ion(Te) = n_ion × mobility × E / L

Solving for n_ion:
n_ion = (ne × n_Ar × k_ion(Te) × L) / (mobility × E)

And Ni/Ne = n_ion / ne = (n_Ar × k_ion(Te) × L) / (mobility × E)
"""

import numpy as np

def scale_ionization(k_ref, Te, Te_ref=1.0, E_ion=15.76):
    return k_ref * np.sqrt(Te/Te_ref) * np.exp(-E_ion * (1/Te - 1/Te_ref))

# Constants
k_ion_ref = 8e-12  # cm³/s at Te = 1.0 eV
E_ion = 15.76  # eV
n_Ar = 1e16  # cm⁻³ (rough estimate)
L = 0.45  # m = 45 cm
mobility = 3057  # cm²/(V·s) for Ar+

print('=' * 80)
print('WHY HIGH E-FIELD DOESN\'T FULLY SOLVE THE PROBLEM')
print('=' * 80)
print()

print('At steady state: Ion production = Ion loss')
print('  Production: R_ion = ne × n_Ar × k_ion(Te)')
print('  Loss: R_loss = n_ion × (mobility × E / L)')
print()
print('Solving for Ni/Ne ratio:')
print('  Ni/Ne = (n_Ar × k_ion(Te) × L) / (mobility × E)')
print()
print('KEY INSIGHT:')
print('  - k_ion(Te) scales exponentially: ~ exp(-15.76/Te)')
print('  - E scales linearly in denominator')
print('  - Exponential >> Linear!')
print()

print('=' * 80)
print('NUMERICAL EXAMPLE')
print('=' * 80)
print()

# Compare different scenarios
Te_low = 1.24  # Current
E_low = 132.6  # Current

Te_high = 2.5  # High
E_high = 300   # High

k_ion_low = scale_ionization(k_ion_ref, Te_low)
k_ion_high = scale_ionization(k_ion_ref, Te_high)

# Ni/Ne ratio (proportional to)
factor_low = (k_ion_low * L) / (mobility * E_low)
factor_high = (k_ion_high * L) / (mobility * E_high)

ratio_increase = factor_high / factor_low

print(f'Scenario 1: Te = {Te_low:.2f} eV, E = {E_low:.1f} V/cm')
print(f'  k_ion = {k_ion_low:.2e}')
print(f'  Ni/Ne factor: {factor_low:.2e}')
print()

print(f'Scenario 2: Te = {Te_high:.2f} eV, E = {E_high:.1f} V/cm')
print(f'  k_ion = {k_ion_high:.2e}')
print(f'  Ni/Ne factor: {factor_high:.2e}')
print()

print(f'Ratio of Ni/Ne factors: {ratio_increase:.1f}×')
print()

if ratio_increase > 1:
    print(f'✗ Even with E increased {E_high/E_low:.2f}×, Ni/Ne is {ratio_increase:.1f}× HIGHER!')
    print()
    print('Why?')
    print(f'  - k_ion increased {k_ion_high/k_ion_low:.0f}× (exponential with Te)')
    print(f'  - E increased {E_high/E_low:.2f}× (linear scaling)')
    print(f'  - Net effect: {(k_ion_high/k_ion_low)/(E_high/E_low):.0f}× more ions!')

print()
print('=' * 80)
print('CONCLUSION')
print('=' * 80)
print()
print('User\'s insight about E-field helping is CORRECT!')
print('  - Increasing E from 132.6 to 300 V/cm reduces Ni by 2.26×')
print()
print('But the exponential scaling of ionization with Te DOMINATES:')
print('  - Increasing Te from 1.24 to 2.5 eV boosts ionization by 859×!')
print('  - E-field boost of 2.26× cannot compensate for 859× more production')
print()
print('The fundamental constraint remains:')
print('  Exponential (ionization) >> Linear (E-field removal)')
print()
print('Would need MASSIVE E-field (thousands of V/cm?) to compensate,')
print('which would be unphysical for this discharge!')
