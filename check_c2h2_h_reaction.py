#!/usr/bin/env python3
"""
Check the thermochemistry and rate constant for H + C2H2 → C2 + H2 + H
"""

print("="*80)
print("ANALYSIS OF H + C2H2 → C2 + H2 + H REACTION")
print("="*80)

print("\n1. REACTION TYPE:")
print("   This is NOT an abstraction reaction (which would be H + C2H2 → C2H + H2)")
print("   This is breaking C2H2 into C2 + removing both H atoms")
print()

print("2. THERMOCHEMISTRY (Standard Heats of Formation at 298K):")
print("   C2H2 (acetylene):  ΔHf° = +227 kJ/mol  (stable, triple bond)")
print("   C2 (diatomic):     ΔHf° = +830 kJ/mol  (very unstable radical)")
print("   H2:                ΔHf° = 0 kJ/mol")
print("   H:                 ΔHf° = +218 kJ/mol")
print()

# Calculate reaction enthalpy
DH_C2H2 = 227  # kJ/mol
DH_C2 = 830
DH_H2 = 0
DH_H = 218

# For H + C2H2 → C2 + H2 + H
DH_rxn = (DH_C2 + DH_H2 + DH_H) - (DH_H + DH_C2H2)
print(f"   ΔHrxn = (C2 + H2 + H) - (H + C2H2)")
print(f"         = ({DH_C2} + {DH_H2} + {DH_H}) - ({DH_H} + {DH_C2H2})")
print(f"         = {DH_rxn:+.0f} kJ/mol")
print()
print(f"   >>> This reaction is HIGHLY ENDOTHERMIC ({DH_rxn:+.0f} kJ/mol) <<<")
print()

print("3. ACTIVATION ENERGY:")
print("   For endothermic reactions: Ea ≥ ΔHrxn")
print(f"   So Ea ≥ {DH_rxn} kJ/mol = {DH_rxn/8.314*1000:.1f} K")
print()
print("   At Tgas = 570 K, this barrier is INSURMOUNTABLE")
print("   This reaction should essentially NOT occur at 570K")
print()

print("4. RATE CONSTANT IN MODEL:")
print("   k = 1.0e-11 cm³/s (constant, no temperature dependence)")
print("   Reference: Baulch et al. (2005)")
print()
print("   >>> This appears to be WRONG! <<<")
print()

print("5. EXPECTED RATE WITH ARRHENIUS FORM:")
print("   If k = A * exp(-Ea/RT), with Ea ~ 603 kJ/mol:")
import numpy as np

T = 570  # K
R = 8.314  # J/mol/K
Ea = 603e3  # J/mol
A = 1e-11  # Assume pre-exponential

k_expected = A * np.exp(-Ea / (R * T))
print(f"   k(570K) ~ {A:.1e} * exp(-{Ea/1e3:.0f}/(8.314*570))")
print(f"           ~ {k_expected:.2e} cm³/s")
print()
print(f"   This is {1e-11/k_expected:.1e} times SMALLER than the model value!")
print()

print("="*80)
print("CONCLUSION")
print("="*80)
print()
print("The reaction H + C2H2 → C2 + H2 + H:")
print("  1. Is HIGHLY endothermic (+603 kJ/mol)")
print("  2. Should have massive activation barrier (>600 kJ/mol)")
print("  3. Should be negligible at Tgas=570K")
print("  4. But accounts for 95% of C2 production in the model!")
print()
print(">>> THIS RATE CONSTANT IS ALMOST CERTAINLY WRONG <<<")
print()
print("Possible explanations:")
print("  A) This is actually a different reaction (misidentified)")
print("  B) The rate is for high-temperature combustion (>1500K), not 570K")
print("  C) There's a typo in the reaction stoichiometry")
print("  D) This represents some other pathway (wall-catalyzed?)")
print()
print("RECOMMENDATION: Check Baulch et al. (2005) directly to verify this rate")
print("                and see at what temperature range it's valid.")
