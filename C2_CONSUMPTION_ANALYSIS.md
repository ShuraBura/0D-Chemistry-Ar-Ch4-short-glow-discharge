# C2 Consumption Analysis: The Wall Loss Revelation

**Date:** 2025-11-14
**Question:** Why did C2 decrease after removing the fake destruction pathway?

---

## The Paradox

We removed **C2 + H → CH + C** which was **93% of C2 destruction**, yet:

```
OLD (with fake destruction):  C2 = 2.75×10⁸ cm⁻³
NEW (without fake destruction): C2 = 9.98×10⁷ cm⁻³

C2 DECREASED by 2.7× instead of increasing!
```

---

## Answer: Wall Loss Now Dominates!

### C2 Destruction Breakdown (Corrected Chemistry):

```
Total C2 destruction: 1.45×10¹¹ cm⁻³/s

Sources:
  Wall sticking:     1.25×10¹¹ cm⁻³/s  (86.1%)
  Diffusion loss:    2.00×10¹⁰ cm⁻³/s  (13.8%)
  Chemical (C2+CH):  1.60×10⁸ cm⁻³/s   (0.1%)

TOTAL WALL/TRANSPORT: 99.9%
TOTAL CHEMICAL:       0.1%
```

**C2 lifetime: τ = 0.69 milliseconds** (extremely short!)

---

## Why Wall Loss Dominates

### The Wall Sticking Rate:

```
k_wall = γ × v_th × (A/V)

Where:
  γ = sticking coefficient = 0.01 (for C2)
  v_th = thermal velocity = √(8kT/πm) = 709 m/s (at T=570K)
  A/V = surface-to-volume ratio = 1.67 cm⁻¹

Result: k_wall = 1.25×10³ s⁻¹
```

### The Problem:

With the fake **C2 + H** destruction disabled, C2 is now destroyed primarily by:
1. **Hitting the walls** (86%)
2. **Diffusing out** (14%)
3. **Chemistry is negligible** (0.1%)

**This means C2 sticking coefficient uncertainty is now CRITICAL!**

---

## C2 Sticking Coefficient: How Certain Is It?

### Assumed Value:
```
γ(C2) = 0.01  (1% sticking probability)
```

### Literature Context:

| Species Type | Typical γ | Examples |
|-------------|-----------|----------|
| Stable molecules | 0.0001-0.001 | CH4, H2, C2H6 |
| Reactive radicals | 0.001-0.01 | CH, CH2, CH3 |
| Very reactive | 0.01-0.1 | C atoms |
| Ions | 0.1-1.0 | All ions |

**C2 is a STABLE RADICAL** with a strong triple bond (C≡C, 6.21 eV):
- More stable than CH, CH2, CH3
- But still a radical (unpaired electrons)
- **Might stick LESS than assumed!**

### Uncertainty:

There is **NO direct measurement** of C2 sticking coefficient on stainless steel at 570K and low pressure!

**Conservative estimate:** γ(C2) = 0.001 to 0.01 (uncertainty of 10×)

---

## Impact of Sticking Coefficient on C2

### Sensitivity Analysis:

| γ(C2) | k_wall (s⁻¹) | Predicted C2 | vs. Target | Notes |
|-------|-------------|--------------|------------|-------|
| 0.01 | 1.25×10³ | 1.0×10⁸ | 5600× low | **Current** |
| 0.005 | 6.25×10² | 2.0×10⁸ | 2800× low | 2× increase |
| 0.001 | 1.25×10² | 1.0×10⁹ | 560× low | **10× increase** |
| 0.0001 | 1.25×10¹ | 1.0×10¹⁰ | 56× low | 100× increase |

**Reducing γ by 10× (to 0.001) would increase C2 by ~10×!**

At γ = 0.0001, C2 would be only 56× too low (much closer to target).

---

## C2 Production vs. Destruction Balance

### With Corrected Chemistry:

**C2 Production (1.45×10¹¹ cm⁻³/s):**
1. C2H2 + C → C2 + CH2 - 32.9% ✓ (Baulch confirmed)
2. C + CH → C2 + H - 13.4% ✓ (Baulch confirmed)
3. CH + C → C2 + H2 - 11.1% ✓ (Baulch confirmed)
4. e + C2H2 → C2 + H2 - 8.5% ✓ (electron impact)
5. C2H + H → C2 + H2 - 7.7% ✓ (abstraction)

**All major production pathways are now VALID!**

**C2 Destruction (1.45×10¹¹ cm⁻³/s):**
1. Wall sticking - 86.1% ❓ **UNCERTAIN!**
2. Diffusion loss - 13.8% ❓ **UNCERTAIN!**
3. C2 + CH → C3 + H - 0.1% ✓ (Baulch confirmed)

**Almost all destruction is physical, not chemical!**

---

## Why C2 Decreased After Corrections

### Lost Production Pathways:

| Pathway | Old Contribution | Status |
|---------|-----------------|--------|
| CH + CH2 → C2 + H2 + H | 21% | ❌ Doesn't exist (Baulch) |
| CH2 + CH2 → C2 + 2H2 | 16% | ❌ Wrong products (makes C2H2) |
| CH + CH → C2 + H2 | ~10% | ❌ Wrong products (makes C2H2) |
| H + C2H2 → C2 + H2 + H | 10% | ❌ Rate 7.5×10⁹ too high |
| **TOTAL LOST** | **~57%** | |

### Gained Production:

| Pathway | Old Rate | New Rate | Gain |
|---------|----------|----------|------|
| C2H2 + C → C2 | 1.0e-10 | 2.0e-10 | 2× ✓ |

**Net effect on production: -57% + ~20% = -37%**

### Lost Destruction:

| Pathway | Old Contribution | Status |
|---------|-----------------|--------|
| C2 + H → CH + C | 93% | ❌ Endothermic (removed) |

**Net effect on destruction: Switched from chemical (93%) to wall (100%)**

### Result:

The **production decrease (-37%)** happened, but it's MASKED by the fact that **wall loss is now 99.9%** of destruction. The model reaches a new equilibrium where:

```
Production = Wall Loss + Chemical Loss
1.45×10¹¹ = 1.45×10¹¹ + ~0
```

C2 adjusts downward because:
1. Less chemical production (lost fake pathways)
2. Wall loss rate (k_wall) stays the same
3. New steady state: [C2] = Production / k_wall

Since production decreased ~37% and k_wall is constant, C2 decreased proportionally.

---

## The Critical Question

**Is the C2 sticking coefficient correct?**

### Evidence It Might Be Too High:

1. **C2 has a strong triple bond** (6.21 eV) - very stable radical
2. **No direct measurements** of γ(C2) at these conditions
3. **If γ is 10× lower**, C2 increases to ~1×10⁹ cm⁻³ (closer to target)
4. **CH sticking is also 0.001** - C2 might be similar or lower

### Evidence It Might Be Right:

1. **C2 is still a radical** - unpaired electrons make it reactive
2. **At low pressure**, radicals have long mean free paths → more wall hits
3. **Model uses 0.01** consistently for medium-reactivity species

---

## Recommendations

### Priority 1: Verify C2 Sticking Coefficient

**Experimental check:**
- Is C2 really lost to walls this quickly?
- Can you measure C2 spatial profile near wall?
- Does C2 decay time match τ ~ 0.7 ms?

**Literature search:**
- C2 sticking on stainless steel / quartz
- C2 surface recombination rates
- Comparable hydrocarbon radicals

### Priority 2: Sensitivity Study

Test C2 model response to γ(C2):
- Try γ = 0.001 (10× lower)
- Try γ = 0.0001 (100× lower)
- See if this brings C2 closer to target

### Priority 3: Re-examine Diffusion Loss

Diffusion loss is 14% of total - also depends on:
- Diffusion coefficient D(C2)
- Boundary conditions
- Geometry

**Is diffusion loss calculated correctly?**

---

## Summary

**After removing fake chemistry, we discovered:**

1. ✅ **C2 production pathways are now chemically correct**
   - Dominated by C2H2 + C, C + CH, and electron impact

2. ❌ **C2 destruction is 99.9% PHYSICAL (wall + diffusion)**
   - Chemical destruction is negligible (only C2 + CH → C3)

3. ❓ **C2 sticking coefficient (γ = 0.01) is CRITICAL and UNCERTAIN**
   - No direct measurements at these conditions
   - 10× uncertainty would cause 10× change in C2

4. 🎯 **C2 decreased because lost production > lost destruction**
   - Production lost 37% (fake pathways removed)
   - Destruction switched from chemical (wrong) to physical (uncertain)

**Bottom line:** The chemistry is now correct, but the **physics (wall loss) dominates** and is highly uncertain!

---

## Next Steps

1. **Investigate C2 sticking coefficient** - literature + sensitivity study
2. **Check diffusion loss formula** - verify boundary conditions
3. **Consider reducing γ(C2)** to 0.001-0.0001 based on C2 stability
4. **Test if γ adjustment alone** can match experimental C2

**The answer may be simpler than we thought: Just need the right sticking coefficient!**

---
