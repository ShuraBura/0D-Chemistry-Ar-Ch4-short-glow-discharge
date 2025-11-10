#!/usr/bin/env python3
"""
Compare CH chemistry before and after drift correction
"""

import json

print("=" * 80)
print("CH CHEMISTRY IMPROVEMENT ANALYSIS")
print("=" * 80)
print()

# Load old result (wrong drift 3.2e17)
with open('optimization_results_charge_balanced/best_f2695.4.json', 'r') as f:
    old = json.load(f)

# Load new result (corrected drift 1.06e17)
with open('optimization_results_charge_balanced/best_f61.8.json', 'r') as f:
    new = json.load(f)

print("COMPARISON: Old (drift=3.2e17) vs New (drift=1.06e17)")
print("-" * 80)
print()

print("PLASMA PARAMETERS:")
print(f"{'Parameter':<15} {'Old (3.2e17)':<20} {'New (1.06e17)':<20} {'Change':<15}")
print("-" * 70)
print(f"{'Te (eV)':<15} {old['Te']:<20.3f} {new['Te']:<20.3f} {new['Te']/old['Te']:<15.2f}×")
print(f"{'Ne (cm⁻³)':<15} {old['Ne']:<20.2e} {new['Ne']:<20.2e} {new['Ne']/old['Ne']:<15.2f}×")
print(f"{'E (V/cm)':<15} {old['E_field']:<20.1f} {new['E_field']:<20.1f} {new['E_field']/old['E_field']:<15.2f}×")
print(f"{'Ni/Ne':<15} {old['Ni_over_Ne']:<20.3f} {new['Ni_over_Ne']:<20.3f}")
print()

print("SPECIES DENSITIES (target comparison):")
print(f"{'Species':<15} {'Old':<20} {'New':<20} {'Improvement':<15}")
print("-" * 70)

H_target = 3.82e14  # Will update to 2.52e14
CH_target = 1.0e9
C2_target = 5.6e11

H_old_ratio = old['all_densities']['H'] / H_target
H_new_ratio = new['all_densities']['H'] / H_target
CH_old_ratio = old['all_densities']['CH'] / CH_target
CH_new_ratio = new['all_densities']['CH'] / CH_target
C2_old_ratio = old['all_densities']['C2'] / C2_target
C2_new_ratio = new['all_densities']['C2'] / C2_target

print(f"{'H (ratio)':<15} {H_old_ratio:<20.2f}× {H_new_ratio:<20.2f}× {(H_old_ratio-H_new_ratio)/H_old_ratio*100:+.1f}%")
print(f"{'CH (ratio)':<15} {CH_old_ratio:<20.2f}× {CH_new_ratio:<20.2f}× {(CH_old_ratio-CH_new_ratio)/CH_old_ratio*100:+.1f}%")
print(f"{'C2 (ratio)':<15} {C2_old_ratio:<20.2f}× {C2_new_ratio:<20.2f}× {(C2_old_ratio-C2_new_ratio)/C2_old_ratio*100:+.1f}%")
print()

print("OBJECTIVE FUNCTION:")
print(f"  Old: f(x) = 2695.4")
print(f"  New: f(x) = 61.8")
print(f"  Improvement: 44× better!")
print()

print("=" * 80)
print("KEY CH CHEMISTRY RATES")
print("=" * 80)
print()

# Extract key rates
print("CH PRODUCTION PATHWAYS:")
print(f"{'Reaction':<35} {'Old rate':<15} {'New rate':<15} {'Change':<10}")
print("-" * 75)

# Approximate from balance analysis:
# Old: CH2+H→CH+H2 = 77.2% of 8.445e14 = 6.52e14
# New: CH2+H→CH+H2 = 30.9% of 3.137e14 = 9.69e13

print(f"{'CH2+H→CH+H2':<35} {'6.5e14 (77%)':<15} {'9.7e13 (31%)':<15} {'6.7× less':<10}")
print(f"{'e+CH4→CH+H2+H & e+CH4→CH+H+H2':<35} {'1.9e14 (23%)':<15} {'2.2e14 (69%)':<15} {'1.1× more':<10}")
print()

print("CH CONSUMPTION (dominant pathway):")
print(f"{'CH+CH4→C2H4+H':<35} {'~3.3e15 (86%)':<15} {'7.3e14 (88%)':<15} {'4.5× less':<10}")
print()

print("=" * 80)
print("WHY DID CH IMPROVE?")
print("=" * 80)
print()
print("With reduced H drift (3.2e17 → 1.06e17 = 3× less):")
print()
print("1. H density decreased 3× (4.2e13 → 1.4e13)")
print("   → Less H available for CH2+H→CH+H2")
print()
print("2. CH2+H→CH+H2 rate dropped 6.7×")
print("   → This was the DOMINANT CH source (77% → 31%)")
print()
print("3. Total CH production dropped 2.7×")
print("   → CH density dropped 4.5× (1.22e10 → 2.70e9)")
print()
print("4. CH now only 2.7× target (was 12.2×)")
print("   → **Major improvement!**")
print()

print("=" * 80)
print("KEY INSIGHT")
print("=" * 80)
print()
print("The excessive H drift was creating too much H, which drove the")
print("CH2+H→CH+H2 reaction, overproducing CH. By reducing the drift")
print("to physically realistic levels, the CH chemistry is now much")
print("better balanced.")
print()
print("GOAL CHEMISTRY TO MAINTAIN:")
print("  - e+CH4 dissociation as primary CH source (~70%)")
print("  - CH2+H→CH+H2 as secondary source (~30%)")
print("  - CH+CH4→C2H4+H as dominant sink (~88%)")
