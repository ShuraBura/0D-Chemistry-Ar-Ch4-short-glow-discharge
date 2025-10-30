"""
Comprehensive analysis of reaction rate temperature dependencies
and missing chemistry in the Ar-CH4 plasma model
"""

print("=" * 100)
print("COMPREHENSIVE CHEMISTRY ANALYSIS: Temperature Dependencies and Missing Pathways")
print("=" * 100)

analysis = """

## CURRENT STATUS:

✓ Group 1 (Electron-Impact Excitation): Te-dependent
✓ Group 2 (Electron-Impact Ionization): Te-dependent
✓ Group 6 (Dissociative Recombination): Te-dependent
✓ Group 11 (Diffusion Losses): Tg-dependent (D ~ Tg^1.75)

## GROUPS THAT NEED TEMPERATURE DEPENDENCE:

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. GROUP 7: NEUTRAL-NEUTRAL REACTIONS (CRITICAL!)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Currently: 65 reactions, ALL CONSTANT at 400K

Problem: These are thermally-activated radical reactions!
- Many have activation barriers (Ea ~ 0-5 kcal/mol)
- Rates vary significantly with Tg: k(T) = A × T^n × exp(-Ea/(R×T))
- For Tg = 300-500K, rates can change by factors of 2-10×

Examples that MUST be Tg-dependent:
  • H + CH4 → CH3 + H2         (Ea ~ 10 kcal/mol, strong Tg dependence)
  • CH3 + CH3 → C2H6           (recombination, Tg^(-0.5) dependence)
  • CH3 + H → CH2 + H2         (Ea ~ 15 kcal/mol)
  • CH + H2 → CH2 + H          (Ea ~ 1.3 kcal/mol)

Impact: If you vary Tg = 300-500K and rates are frozen at 400K:
  - Low Tg (300K): Will overestimate radical reactions (too fast)
  - High Tg (500K): Will underestimate radical reactions (too slow)

PRIORITY: HIGH - Essential for accurate Tg studies

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
2. GROUP 8: TERMOLECULAR RECOMBINATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Currently: 3 reactions, constant
  • H + H + M → H2 + M
  • CH3 + CH3 + M → C2H6 + M
  • CH3 + H + M → CH4 + M

These scale as: k ~ T^(-n) where n ~ 1-3 (third body efficiency decreases with T)
Typical: k ~ T^(-1.5) to T^(-2.5)

Impact: Moderate (factor of 2-3× for 300-500K)
PRIORITY: MEDIUM

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
3. GROUPS 3-4: Ar* REACTIONS AND PENNING IONIZATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Currently: ~40 reactions, constant

Weak Tg dependence: k ~ sqrt(Tg) (thermal velocity dependence)
  - Ar* + CH4 → products
  - Ar* + H2 → products
  - Penning ionization: Ar* + M → M+ + Ar + e

Impact: Small (factor of ~1.3× for 300-500K)
PRIORITY: LOW - Can ignore for now

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
4. GROUP 5: ION-NEUTRAL REACTIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Currently: 12 reactions, constant

Most are Langevin-type (collision-limited), very fast
Some weak Tg dependence: k ~ Tg^(-0.16) (from collision theory)

Examples:
  • Ar+ + CH4 → CH3+ + H + Ar   (collision-limited, ~1e-9 cm³/s)
  • CH3+ + CH4 → C2H5+ + H2

Impact: Very small (factor of ~1.1× for 300-500K)
PRIORITY: LOW - These are already near collision rate

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
5. GROUP 9: WALL STICKING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Currently: Wall loss rates, constant

Should scale with thermal velocity: k_stick ~ sqrt(Tg)

Impact: Small to moderate (factor of ~1.3× for 300-500K)
PRIORITY: LOW-MEDIUM - Important for surface chemistry studies


━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
MISSING CHEMISTRY PATHWAYS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

A. H3+ FORMATION AND REACTIONS (IMPORTANT!)
   Missing:
   • H2+ + H2 → H3+ + H         (very fast, ~1e-9 cm³/s)
   • CH5+ → CH3+ + H2           (should this be reversible?)
   • H3+ + CH4 → CH5+ + H2      (proton transfer)

   Issue: H3+ appears in mobilities/drift but no formation pathway!

B. EXCITED STATE QUENCHING
   Present: H2*, CH4*, C2H2*, C3H4* mentioned
   Missing explicit quenching reactions:
   • H2* + M → H2 + M           (vibrational relaxation)
   • CH4* + M → CH4 + M

C. NEGATIVE ION DETACHMENT (Already present)
   Present:
   • H- + various positive ions
   • CH3- + various positive ions
   Good coverage ✓

D. H2+ AND CH2+ CHEMISTRY
   Issue: CH2+ and H2+ appear in ion list but limited formation/destruction
   • e + H2 → H2+ + 2e          (missing?)
   • H2+ + CH4 → CH4+ + H2      (missing?)

E. LARGER HYDROCARBON IONS (C3+, C4+ species)
   Currently limited C3/C4 ion chemistry
   For long-residence plasmas, might form:
   • C3Hx+ species
   • C4Hx+ species

F. ATTACHMENT/DETACHMENT BALANCE
   Present:
   • e + CH4 → CH3- + H         (attachment)
   • CH3- + positive ions       (mutual neutralization)

   Missing:
   • CH3- + M → CH3 + e + M     (collisional detachment)
   • H- + M → H + e + M         (collisional detachment)


━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
RECOMMENDATIONS (Priority Order):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. **CRITICAL - Group 7 Tg-dependence**
   Add Arrhenius forms for key radical reactions
   Estimate impact: 2-10× rate changes for Tg = 300-500K
   Effort: Medium (need to find Ea and n for ~20 key reactions)

2. **HIGH - H3+ chemistry**
   Add formation pathway: H2+ + H2 → H3+ + H
   Add key H3+ reactions
   Effort: Low (2-3 reactions)

3. **MEDIUM - Group 8 Tg-dependence**
   Add T^(-n) scaling for termolecular recombination
   Effort: Low (3 reactions)

4. **LOW - H2+ and CH2+ chemistry**
   Add missing formation/destruction pathways
   Effort: Low (3-5 reactions)

5. **OPTIONAL - Group 9 wall sticking Tg^0.5**
   Only if studying surface chemistry effects
   Effort: Low

6. **OPTIONAL - Excited state quenching**
   Only if tracking excited state densities explicitly
   Effort: Low-Medium

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SPECIFIC QUESTIONS FOR YOU:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. Do you plan to vary Tg in your simulations?
   - If YES → Must implement Group 7 Tg-dependence
   - If NO → Can skip for now

2. Do you see significant H3+ density in your results?
   - If YES → Need to add H3+ formation chemistry
   - If NO but it appears in species list → Should add anyway

3. What is your Tg range? (300K? 400K? 500K?)
   - This determines how critical Tg-dependence is

4. Are you matching experimental data at specific conditions?
   - If YES → Current rates (calibrated at Te~1eV, Tg~400K) are fine
   - If doing parameter sweeps → Need temperature dependence

5. Residence time / pressure regime?
   - Short (<1 ms) → Current chemistry sufficient
   - Long (>10 ms) → May need more C3+/C4+ chemistry
"""

print(analysis)

print("=" * 100)
print("NEXT STEPS:")
print("=" * 100)
print("""
Based on your answers, I can implement:
1. Tg-dependent neutral-neutral rates (Group 7) with Arrhenius forms
2. H3+ formation chemistry
3. Temperature scaling for termolecular reactions
4. Any other missing pathways you identify

What would you like to prioritize?
""")
