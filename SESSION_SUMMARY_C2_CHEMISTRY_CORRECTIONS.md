# Session Summary: C2 Chemistry Corrections from Baulch (2005)

**Date:** 2025-11-14
**Branch:** `claude/continue-next-section-01MTt2i7BNfZfUd7FMqPSADa`

---

## What We Accomplished

### 1. Identified Critical Chemistry Errors

We verified the model's C2 reactions against **Baulch et al. (2005)** - the gold standard for combustion kinetics - and found **major errors**:

| Reaction | Model Status | Baulch Status | Impact |
|----------|--------------|---------------|--------|
| **C2 + H → CH + C** | k = 9.6e-11 | ❌ NOT IN BAULCH (endothermic) | **93% of C2 destruction** |
| **CH + CH2 → C2 + H2 + H** | k = 1.2e-10 | ❌ NOT IN BAULCH ("garbage") | **21% of C2 production** |
| **CH2 + CH2 → C2 + 2H2** | k = 1.0e-11 | ❌ WRONG PRODUCTS | **16% of C2 production** |
| **H + C2H2 → C2 + H2 + H** | k = 1.0e-11 (const) | **7.5 billion × too high!** | **10% of C2 production** |
| **CH + CH → C2 + H2** | k = 1.0e-10 | ❌ WRONG PRODUCTS (should be C2H2) | Production pathway |
| **C2H2 + C → C2 + CH2** | k = 1.0e-10 | ✓ Should be 2.0e-10 | **18% of C2 production** |
| **CH + H → C + H2** | k = 1.2e-10 | ✓ Should be 2.0e-10 | CH destruction |

### Summary:
- **Top 3 C2 production pathways were WRONG** (21% + 18% + 16% = 55%)
- **Main C2 destruction pathway doesn't exist** at 570K (93% of destruction!)
- Several rates were off by factors of 2, or in one case, **7.5 billion!**

---

## 2. Implemented Corrections

### File: `define_rates.py`

**Disabled invalid reactions:**
```python
# Line 255-256: C2 + H → CH + C
k['C2_H_CH_C_cm3_7_6'] = 0.0  # Was 93% of C2 destruction - endothermic!

# Line 289-290: CH + CH2 → C2 + H2 + H
k['CH2_CH_C2_H2_H_cm3_7_26'] = 0.0  # Was 21% of C2 production - doesn't exist!
```

**Fixed rate constants:**
```python
# Line 252: CH + H → C + H2
k['CH_H_C_H2_cm3_7_3'] = 2.0e-10  # Was 1.2e-10 (Baulch 2005)

# Line 269: C2H2 + C → C2 + CH2
k['C2H2_C_C2_CH2_cm3_7_19'] = 2.0e-10  # Was 1.0e-10 (Baulch 2005)

# Line 317-319: H + C2H2 → C2 + H2 + H (Arrhenius form)
k['C2H2_H_C2_H2_H_cm3_7_50'] = 1.67e-14 * Tgas**1.64 * np.exp(-15250/Tgas)
# Was constant 1.0e-11 (7.5 billion × too high at 570K!)
```

### File: `build_reactions.py`

**Fixed reaction stoichiometry:**
```python
# Line 258-259: CH2 + CH2 → C2H2 + H2 (NOT C2 + H2 + H2!)
push(sto('CH2', 2), sto('C2H2', 1, 'H2', 1), ...)  # Corrected products

# Line 141-142: CH + CH → C2H2 (NOT C2 + H2!)
push(sto('CH', 2), sto('C2H2', 1), ...)  # Corrected products

# Line 245-246: CH + CH → C2H2 (NOT C2 + H2!) - duplicate in neutral section
push(sto('CH', 2), sto('C2H2', 1), ...)  # Corrected products
```

---

## 3. Test Results After Corrections

**Test run with corrected chemistry:**

```
Conditions: P = 400 mTorr, T = 570 K, ne = 2.3×10⁹ cm⁻³, Te = 1.3 eV

RESULTS:
  H  = 1.57×10¹⁴ cm⁻³  (target: 2.30×10¹⁴, error: 31.8%) ✓
  CH = 8.00×10⁹ cm⁻³   (target: 1.34×10⁹, error: 497%) ✗
  C2 = 9.98×10⁷ cm⁻³   (target: 5.60×10¹¹, error: 100%) ✗

Comparison:
  Old C2 (wrong rates): 2.75×10⁸ cm⁻³
  New C2 (corrected):   9.98×10⁷ cm⁻³
  Change: 0.4× (DECREASED!)

RMS Error: 293.4%
```

### Unexpected Result: C2 Decreased!

Despite removing the dominant (93%) fake C2 destruction pathway, **C2 actually decreased by 2.7×**.

**Why?**
1. ❌ Removed C2 + H → CH + C (should increase C2) ✓
2. ❌ Removed CH + CH2 → C2 (decreases C2 production)
3. ❌ Removed CH + CH → C2 + H2 (now produces C2H2 instead)
4. ❌ Removed CH2 + CH2 → C2 (now produces C2H2 instead)
5. ❌ Fixed H + C2H2 → C2 rate (~99.9% slower)
6. ✓ Increased C2H2 + C → C2 rate (2× faster)

**Net effect:** Loss of fake production pathways outweighed removal of fake destruction!

---

## 4. Key Implications

### The Model's C2 Chemistry Was Fundamentally Broken

1. **C2 was artificially suppressed** by a non-existent endothermic destruction pathway (C2 + H)
2. **C2 was artificially inflated** by multiple non-existent production pathways
3. **The two wrongs partially canceled out** - giving deceptively reasonable C2 values

### After Corrections:

- **C2 is now 5600× too low** (9.98×10⁷ vs. target 5.60×10¹¹ cm⁻³)
- **CH is now 6× too high** (8.00×10⁹ vs. target 1.34×10⁹ cm⁻³)
- **H is 32% too low** (reasonable)

### Real C2 Formation Pathways at 570K:

Based on Baulch (2005), the VALID C2 formation pathways are:

✓ **Electron impact:**
- e + C2H2 → C2 + H2 + e

✓ **Ion chemistry:**
- C2H+ + e → C2 + H

✓ **C-atom reactions:**
- C + CH → C2 + H
- C2H2 + C → C2 + CH2 (rate corrected)
- C + C2H3 → C2 + CH3

✓ **Radical abstraction:**
- C2H + H → C2 + H2

✓ **Three-body recombination (minor):**
- C + C + M → C2 + M

❌ **NO LONGER VALID:**
- CH + CH → C2 + H2 (produces C2H2 instead!)
- CH + CH2 → C2 (doesn't exist)
- CH2 + CH2 → C2 (produces C2H2 instead!)
- H + C2H2 → C2 (rate was 7.5×10⁹ too high)

---

## 5. What's Missing?

The corrected model still cannot match experimental C2 (off by 5600×).

**Possible explanations:**

1. **Missing C2 production pathway(s)** at low temperature
   - Vibrationally-excited species reactions?
   - State-specific chemistry (C2 electronic states)?

2. **Electron impact rates too low**
   - e + C2H2 → C2 may be more efficient than thought

3. **Ion chemistry underestimated**
   - More C2H+ formation?
   - Other ion pathways to C2?

4. **Spatial effects** (0D limitation)
   - C2 formed in plasma core, measured at edge?
   - Transport/diffusion creates apparent mismatch?

5. **Temperature effects**
   - Local Tgas > 570K in discharge?
   - Vibrational temperatures Tv >> Tgas?

---

## 6. Files Modified

### Core Chemistry Files:
- ✅ `define_rates.py` - 7 rate constant corrections
- ✅ `build_reactions.py` - 3 stoichiometry corrections

### Documentation Created:
- ✅ `CRITICAL_C2_CHEMISTRY_ERRORS.md` - Detailed error analysis
- ✅ `C2_REACTIONS_ANALYSIS.md` - Complete C2 reaction catalog
- ✅ `RATE_AND_LITERATURE_REVIEW.md` - Literature verification (from previous session)

### Analysis Scripts:
- ✅ `test_corrected_chemistry.py` - Quick test of corrections
- ✅ `recalculate_c2_corrected.py` - Detailed pathway analysis (partial)

---

## 7. Next Steps / Recommendations

### Immediate Priorities:

1. **Run detailed C2 pathway analysis** with corrected chemistry
   - Identify dominant production/destruction routes
   - Quantify each pathway's contribution

2. **Search for missing C2 formation pathways**
   - Literature review for low-T C2 chemistry
   - Check plasma-specific reactions
   - Look for state-specific rates

3. **Verify electron impact rates**
   - e + C2H2 → C2 efficiency
   - Check Te dependence

4. **Consider spatial/transport effects**
   - Is 0D model adequate?
   - Need 1D radial model?

### Long-term:

- Implement vibrational state chemistry (C2(v), CH(v))
- Add electronic state chemistry (C2 a³Πu vs X¹Σ+g)
- Consider 1D spatial model
- Benchmark against other Ar/CH4 plasma studies

---

## 8. Key Learnings

1. **Trust the literature** - Baulch (2005) is the gold standard for a reason
2. **Verify all reactions** - Even "confirmed" reactions had wrong products
3. **Two wrongs don't make a right** - Fake production + fake destruction = misleading results
4. **Chemistry matters more than fitting** - Correcting fundamental errors is more important than matching targets with wrong chemistry

---

## Summary

We've completed a major cleanup of the C2 chemistry, removing non-physical reactions and correcting rates against Baulch (2005). The model now has **chemically correct C2 pathways**, but reveals that C2 is severely under-produced compared to experiments. This sets the stage for identifying what's truly missing - whether it's:
- Missing chemical pathways
- Underestimated reaction rates
- Physics beyond 0D approximation
- State-specific chemistry effects

**The foundation is now solid. Time to build the right chemistry on top of it!**

---

**Status:** ✅ All chemistry corrections committed and tested
**Branch:** `claude/continue-next-section-01MTt2i7BNfZfUd7FMqPSADa`
**Ready for:** Push and continue investigation

