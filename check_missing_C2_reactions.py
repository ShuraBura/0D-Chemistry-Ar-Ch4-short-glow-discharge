#!/usr/bin/env python3
"""
Check for potentially missing C2-producing reactions
"""

print("=" * 80)
print("CURRENT C2-PRODUCING REACTIONS IN NETWORK")
print("=" * 80)
print()

current_reactions = [
    ("CH + CH → C2 + H2", "CH_CH_C2_H2_cm3_5_4", "Dominant pathway"),
    ("CH + CH → C2 + H2", "CH_CH_C2_H2_cm3_7_44", "Duplicate?"),
    ("C + CH → C2 + H", "C_CH_C2_H_cm3_7_4", "Secondary"),
    ("CH + C → C2 + H", "CH_C_C2_H_cm3_7_9", "Secondary"),
    ("CH + C → C2 + H", "CH_C_C2_H_cm3_7_43", "Duplicate?"),
    ("CH + C → C2 + H2", "CH_C_C2_H2_cm3_7_24", "Alternative"),
    ("CH2 + CH → C2 + H2 + H", "CH2_CH_C2_H2_H_cm3_7_26", "Minor"),
    ("CH2 + CH2 → C2 + 2H2", "CH2_CH2_C2_H2_H2_cm3_7_58", "Minor"),
    ("C + C + M → C2 + M", "C_C_M_C2_M_cm6_7_64", "Three-body"),
]

for rxn, tag, note in current_reactions:
    print(f"{rxn:<30} [{tag}]")
    print(f"  Note: {note}")
print()

print("=" * 80)
print("POTENTIALLY MISSING C2-PRODUCING REACTIONS")
print("=" * 80)
print()

missing_reactions = [
    ("Electron-impact dissociation of C2Hx", [
        "e + C2H2 → C2 + H + H",
        "e + C2H2 → C2 + H2",
        "e + C2H4 → C2 + 2H2",
        "e + C2H6 → C2 + 3H2",
    ]),

    ("Ion recombination", [
        "C2H+ + e → C2 + H",
        "C2H2+ + e → C2 + H2",
        "C2H3+ + e → C2 + H2 + H",
    ]),

    ("Excited state reactions", [
        "Ar* + C2H2 → C2 + H2 + Ar",
        "Ar* + C2H4 → C2 + 2H2 + Ar",
    ]),

    ("Higher-order neutral reactions", [
        "C2H + H → C2 + H2",
        "C2H2 + C → C2 + CH2",
    ]),
]

for category, reactions in missing_reactions:
    print(f"\n{category}:")
    for rxn in reactions:
        print(f"  {rxn}")

print()
print("=" * 80)
print("ANALYSIS")
print("=" * 80)
print()

print("Current network has:")
print("  ✓ Neutral-neutral reactions (CH+CH, C+CH, etc.)")
print("  ✓ Three-body recombination (C+C+M)")
print()
print("Potentially missing:")
print("  ? Electron-impact dissociation of C2H2 → C2")
print("  ? Ion recombination C2H+ + e → C2")
print("  ? Excited state reactions Ar* + C2H2 → C2")
print()

print("RECOMMENDATION:")
print()
print("1. CHECK: Are e+C2H2→C2+H2 and similar reactions in the network?")
print("   These could be VERY FAST at our Ne ~ 2e9 cm⁻³")
print()
print("2. CHECK: Ion recombination C2H+ + e → C2 + H")
print("   This is typically fast (~1e-7 cm³/s)")
print()
print("3. If missing, ADD these reactions to build_reactions.py")
print()
print("4. If present but not tunable, FORCE them into optimizer")
