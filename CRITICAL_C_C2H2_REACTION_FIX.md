# CRITICAL FIX: C + C2H2 Reaction Products

**Date:** 2025-11-17
**Issue:** Wrong products for C + C2H2 reaction causing artificial C2 inflation

---

## Problem Identified

The model had:
```
C + C2H2 → C2 + CH2  (k = 2.0×10⁻¹⁰ cm³/s)
```

This reaction was:
- Contributing **33% of total C2 production**
- Labeled as having Baulch (2005) citation for rate constant
- **BUT WRONG PRODUCTS!**

---

## Literature Evidence

Web search for "C + C2H2" reactions found:

**Major channel:** C(³P) + C2H2 → **C3H2** (cyclopropenylidene)
**Subsequent step:** C3H2 → **C3H + H**

References cite:
- Ab initio calculations for C + HCCH → C3H2 → C3H + H pathway
- Variational RRKM theory confirming this mechanism
- Experimental studies of carbon atom reactions with acetylene

**NO evidence found for C + C2H2 → C2 + CH2 channel!**

---

## Thermochemical Check

**C + C2H2 → C2 + CH2:**
- ΔH = ΔHf(C2) + ΔHf(CH2) - ΔHf(C) - ΔHf(C2H2)
- ΔH = 830 + 391 - 717 - 227 = **+277 kJ/mol** (endothermic!)
- At T = 570 K: Would require high activation barrier

**C + C2H2 → C3H + H:**
- ΔH = ΔHf(C3H) + ΔHf(H) - ΔHf(C) - ΔHf(C2H2)
- ΔH ~ -50 kJ/mol (exothermic, estimate)
- Barrierless radical addition expected

**Thermochemistry supports C3H + H, not C2 + CH2!**

---

## Fix Applied

Changed build_reactions.py line 220:

**Before:**
```python
push(sto('C2H2', 1, 'C', 1), sto('C2', 1, 'CH2', 1), k['C2H2_C_C2_CH2_cm3_7_19'], ...)
```

**After:**
```python
push(sto('C2H2', 1, 'C', 1), sto('C3H', 1, 'H', 1), k['C2H2_C_C2_CH2_cm3_7_19'], ...)
```

**Comment added:** "FIXED: Literature shows C + C2H2 → C3H + H, NOT C2 + CH2!"

---

## Impact on Results

| Condition | C2 (cm⁻³) | vs Target | Notes |
|-----------|-----------|-----------|-------|
| Before fix | 4.43×10⁸ | 1263× low | With wrong products |
| After fix | 2.98×10⁸ | 1879× low | With correct products |

**C2 decreased by 33%** (as expected, since we removed 33% of production)

**Gap increased from 1263× to 1879×**, but this is the CORRECT chemistry!

---

## Remaining C2 Production Pathways

After fix, C2 production is now dominated by:

1. **CH + C → C2 + H/H2** (73% of production)
   - Rate constants: k ~ 1-1.2×10⁻¹⁰ cm³/s
   - **Status: NO LITERATURE CITATIONS!**
   - **Need verification:** Are these rates correct?

2. **e + C2H2 → C2 + H2** (12.7% of production)
   - Rate constant: k = 4.55×10⁻¹⁰ at Te=1.3 eV
   - E_threshold = 9.0 eV
   - Status: Has citation but threshold may be wrong (should be ~6.3 eV thermochemically)

3. **C2H + H → C2 + H2** (11.5% of production)
   - Rate constant: k = 1×10⁻¹⁰ cm³/s
   - Status: Uncertain citation

---

## Critical Questions

### 1. Is CH + C → C2 + H a real reaction?

**Evidence needed:**
- Literature search for CH + C kinetics
- Thermochemical analysis (is it exothermic?)
- Typical rates for radical-radical recombination

**Current rates:** k ~ 1-1.2×10⁻¹⁰ cm³/s
**Gas kinetic limit:** k ~ 3×10⁻¹⁰ cm³/s at 570 K

→ Rates are ~30% of collision frequency (plausible for barrierless)

**Thermochemistry:**
- CH + C → C2 + H:  ΔH = ΔHf(C2) + ΔHf(H) - ΔHf(CH) - ΔHf(C)
- ΔH = 830 + 218 - 595 - 717 = **-264 kJ/mol** (exothermic!)

→ Reaction is thermochemically favorable!

But: **Still need literature confirmation of rate constant!**

### 2. What is the actual C2 production mechanism?

If current model chemistry is correct (with fix), then:
- Need 1879× more C2 production
- OR 1879× less C2 destruction

**Possible explanations:**

A. **CH + C rate is underestimated by 1000-10000×**
   - Tests showed C2 plateaus at ~3.7×10¹¹ even at 20000× multiplier
   - Precursors [C] and [CH] deplete as rate increases
   - Cannot explain full 1879× gap alone

B. **Missing C2 production pathway**
   - Surface production of C2?
   - Different electron impact channel?
   - Ion-neutral chemistry?

C. **Wall sticking still too high**
   - Current: γ(C2) = 0.001 (based on thermal radical analogy)
   - If γ(C2) = 0.0001: only gives ~10× improvement
   - Cannot explain 1879× gap

D. **Model conditions incorrect**
   - Wrong Te, ne, or gas temperature?
   - Wrong pressure or flow rate?
   - But H and CH match well...

---

## Conclusion

1. **C + C2H2 → C2 + CH2 was WRONG** (now fixed to → C3H + H)

2. **C2 gap increased to 1879×** with correct chemistry

3. **Dominant pathway now CH + C → C2** (73%), but NO CITATIONS

4. **Need urgent literature verification:**
   - CH + C → C2 rate constant
   - e + C2H2 → C2 threshold energy
   - Any missing C2 production pathways

5. **Model self-consistency remains good:**
   - H matches within 1.5×
   - CH matches within 6×
   - Only C2 is way off (species-specific issue)

---

## Next Steps

**Priority 1:** Literature search for CH + C → C2 kinetics
**Priority 2:** Verify e + C2H2 → C2 threshold (9 eV vs 6.3 eV)
**Priority 3:** Search for alternative C2 production mechanisms
**Priority 4:** Consider ion-neutral chemistry (C₂⁺ + e → C2?)

---
