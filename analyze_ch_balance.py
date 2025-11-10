#!/usr/bin/env python3
"""
Analyze CH production and loss balance at different Te values
"""

import numpy as np
import json
from define_rates import define_rates
from build_reactions import build_reactions

# Load best results - high Te vs low Te
with open('optimization_results_charge_balanced/best_f2962.3.json', 'r') as f:
    high_te = json.load(f)

with open('optimization_results_charge_balanced/best_f990.2.json', 'r') as f:
    low_te = json.load(f)

def analyze_ch_balance(data, label):
    print(f"\n{'='*80}")
    print(f"{label}")
    print(f"{'='*80}")

    Te = data['Te']
    Ne = data['Ne']
    E = data['E_field']
    Tgas = 400.0  # K

    # Get densities
    densities = data['all_densities']
    CH = densities['CH']
    CH4 = densities['CH4']
    CH4Plus = densities['CH4Plus']
    CHPlus = densities.get('CHPlus', 0)

    print(f"\nPlasma parameters:")
    print(f"  Te = {Te:.2f} eV")
    print(f"  Ne = {Ne:.2e} cm^-3")
    print(f"  E = {E:.1f} V/cm")
    print(f"  Ni/Ne = {data['Ni_over_Ne']:.2f}")

    print(f"\nCH density: {CH:.2e} cm^-3 ({CH/1e9:.2f}× target)")

    # Get rate constants - define_rates takes a params dict
    params = {
        'E_field': E,
        'L_discharge': 0.05,  # 5 cm
        'mobilities': {},  # Not needed for rate calculation
        'Te': Te,
        'Tg': Tgas
    }
    k = define_rates(params)

    print(f"\n{'='*80}")
    print("CH PRODUCTION PATHWAYS:")
    print(f"{'='*80}")

    total_production = 0

    # 1. e + CH4 → CH + H2 + H (dissociative ionization/excitation)
    if 'e_CH4_CH_H2_H_vib_cm3_1_3' in k:
        rate = k['e_CH4_CH_H2_H_vib_cm3_1_3']
        flux = rate * Ne * CH4
        total_production += flux
        print(f"\n1. e + CH4 → CH + H2 + H (dissociative)")
        print(f"   Rate: {rate:.2e} cm³/s")
        print(f"   Flux: {flux:.2e} cm^-3/s")
        print(f"   [Ne]={Ne:.2e}, [CH4]={CH4:.2e}")

    # 2. e + CH4Plus → CH + products (recombination)
    if 'e_CH4Plus_CH_H2_H_cm3_6_11' in k:
        rate = k['e_CH4Plus_CH_H2_H_cm3_6_11']
        flux = rate * Ne * CH4Plus
        total_production += flux
        print(f"\n2. e + CH4Plus → CH + H2 + H (recombination)")
        print(f"   Rate: {rate:.2e} cm³/s")
        print(f"   Flux: {flux:.2e} cm^-3/s")
        print(f"   [Ne]={Ne:.2e}, [CH4Plus]={CH4Plus:.2e}")

    # 3. CHPlus + e → CH (our new reaction!)
    if 'CHPlus_e_CH_cm3_6_32' in k and CHPlus > 0:
        rate = k['CHPlus_e_CH_cm3_6_32']
        flux = rate * Ne * CHPlus
        total_production += flux
        print(f"\n3. CHPlus + e → CH (NEW recombination)")
        print(f"   Rate: {rate:.2e} cm³/s")
        print(f"   Flux: {flux:.2e} cm^-3/s")
        print(f"   [Ne]={Ne:.2e}, [CHPlus]={CHPlus:.2e}")

    print(f"\n{'='*80}")
    print("CH LOSS PATHWAYS:")
    print(f"{'='*80}")

    total_loss = 0

    # 1. CH wall sticking
    if 'stick_CH_9_3' in k:
        rate = k['stick_CH_9_3']
        flux = rate * CH
        total_loss += flux
        print(f"\n1. CH → wall (sticking)")
        print(f"   Rate: {rate:.2f} s^-1")
        print(f"   Flux: {flux:.2e} cm^-3/s")

    # 2. CH ambipolar/drift loss
    if 'loss_CH_11_9' in k:
        rate = k['loss_CH_11_9']
        flux = rate * CH
        total_loss += flux
        print(f"\n2. CH ambipolar loss")
        print(f"   Rate: {rate:.2f} s^-1")
        print(f"   Flux: {flux:.2e} cm^-3/s")

    # 3. CH ionization (our new reaction!)
    if 'e_CH_CHPlus_2e_cm3_2_10' in k:
        rate = k['e_CH_CHPlus_2e_cm3_2_10']
        flux = rate * Ne * CH
        total_loss += flux
        print(f"\n3. e + CH → CHPlus + 2e (NEW ionization)")
        print(f"   Rate: {rate:.2e} cm³/s")
        print(f"   Flux: {flux:.2e} cm^-3/s")
        print(f"   [Ne]={Ne:.2e}, [CH]={CH:.2e}")

    print(f"\n{'='*80}")
    print(f"TOTALS:")
    print(f"  Production: {total_production:.2e} cm^-3/s")
    print(f"  Loss:       {total_loss:.2e} cm^-3/s")
    print(f"  Net:        {(total_production - total_loss):.2e} cm^-3/s")
    print(f"  P/L ratio:  {total_production/total_loss:.2f}")
    print(f"{'='*80}")

# Analyze both cases
analyze_ch_balance(low_te, "LOW Te CASE (CH=7.91×, Ni/Ne=0.16 - BAD CHARGE BALANCE)")
analyze_ch_balance(high_te, "HIGH Te CASE (CH=13.15×, Ni/Ne=1.59 - BETTER)")

print(f"\n\n{'='*80}")
print("SUMMARY: Why is CH higher at high Te?")
print(f"{'='*80}")
print("""
The key insight:
1. At low Te, dissociative ionization (e + CH4 → CH + ...) is LESS efficient
   because the cross-section drops at low electron energy
2. However, low Te also means low ionization, so Ni/Ne is terrible (0.16 vs 2-7)
3. At high Te, dissociative ionization becomes MUCH more efficient
4. The CH ionization we added (e + CH → CHPlus) helps, but doesn't compensate
   for the increased CH production

Solution: We need to either:
A. Reduce CH production at high Te (limit dissociative ionization)
B. Enhance CH loss at high Te (more CH chemistry or faster ionization)
C. Find alternative chemistry that produces less CH while maintaining charge balance
""")
