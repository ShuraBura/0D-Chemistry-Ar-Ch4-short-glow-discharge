# Chemistry Analysis - Root Cause of CH Problem

## Executive Summary

**Problem:** CH is 46x too high (4.59×10¹⁰ vs target 1.0×10⁹ cm⁻³)

**Root Cause:** **C2H2 → C2 → CH feedback loop**

The diagnostic reveals that the CH problem is NOT primarily from electron-impact on CH4, but from a **downstream chemistry feedback loop** involving C2H2.

---

## The Feedback Loop

```
High C2H2 (1.04×10¹³ cm⁻³)
     ↓ (91% of C2 production)
C2H2 + H → C2 + H2 + H  [Rate: 3.57×10¹⁵ cm⁻³ s⁻¹]
     ↓
High C2 (7.31×10¹¹ cm⁻³, 5.6x too high)
     ↓ (59% of CH production, 62% of C2 loss)
C2 + H → CH + C  [Rate: 2.42×10¹⁵ cm⁻³ s⁻¹]
     ↓
High CH (4.59×10¹⁰ cm⁻³, 46x too high)
```

**Key insight:** The problem cascades from **C2H2** down to CH.

---

## Detailed Chemistry Breakdown

### H Chemistry (0.67x of target - need +49%)

**Loss dominated by ONE reaction:**
- **H + CH4 → CH3 + H2**: 92% of all H loss (2.996×10¹⁷ cm⁻³ s⁻¹)

**Production spread across many reactions:**
- CH2 + CH3 reactions: 22%
- C + CH3 reactions: 19%
- CH + CH4 reactions: 15-8%
- Electron-impact: 5%

**Problem:** H is consumed 60x faster than produced because one reaction dominates loss.

**Solution direction:** Can't reduce H + CH4 → CH3 + H2 much (it's important for overall chemistry). Need to increase H production instead.

---

### CH Chemistry (46x too high - THE PROBLEM!)

**Production dominated by C2-related reactions (94% total!):**
1. **C2 + H → CH + C**: 59% (2.42×10¹⁵ cm⁻³ s⁻¹) ⚠️
2. CH2 + H → CH + H2: 21%
3. C + H → CH: 14%
4. **Electron-impact (e + CH4):** Only 3.5%! (Not the main source!)

**Loss spread across multiple CH4 reactions:**
1. CH + CH4 → C2H4 + H: 19%
2. CH + CH4 → CH2 + CH3: 16%
3. CH wall loss: 11%
4. CH wall sticking: 7%
5. CH + CH3 reactions: 5-4%

**Critical finding:**
- Electron-impact only contributes 3.5% to CH production
- **C2 + H → CH + C contributes 59%** - this is the real culprit!

**Implication:** Reducing electron density (Ne) won't solve the problem because electron-impact is a minor CH source. The problem is the C2 → CH conversion.

---

### C2 Chemistry (5.6x too high)

**Production dominated by ONE reaction:**
- **C2H2 + H → C2 + H2 + H**: 91% (3.57×10¹⁵ cm⁻³ s⁻¹) ⚠️
- All other reactions: <6% each

**Loss:**
- **C2 + H → CH + C**: 62% (feeds CH!)
- C2 wall sticking: 23%
- C2 wall loss: 15%

**Problem:** C2 is produced almost entirely from C2H2, then immediately converts to CH.

---

### C2H2 Chemistry (NOT a target, but the root cause!)

**Current density:** 1.04×10¹³ cm⁻³ (very high!)

**Top C2H2 production pathways (from diagnostic):**
1. CH + CH3 → C2H2 + H2
2. CH3 + CH3 → C2H2 + H2 + H2
3. Various other hydrocarbon reactions

**C2H2 loss:**
- Electron-impact: e + C2H2 → C2 + H2 (creates the C2!)
- Wall losses
- Reactions with radicals

**Problem:** C2H2 is a stable intermediate that accumulates, then feeds into C2 production.

---

## Why Previous Optimizations Didn't Work

### 2-Parameter Optimization (Ne + E only)
- Reduced Ne by 50% (3.3×10⁹ → 1.66×10⁹)
- Only reduced CH by 4% (46x → 44x)
- **Why it failed:** Electron-impact is only 3.5% of CH production!

### 32-Parameter Optimization
- Achieved 97% improvement in objective function
- Stopped at step 7 (unknown final densities)
- Was tuning rates, but without targeting C2H2 specifically

---

## Recommended Solutions

### Option 1: Reduce C2H2 → C2 Conversion (Direct)
**Target reaction:** C2H2 + H → C2 + H2 + H

**Current rate:** Unknown (need to check literature)

**Action:**
- Check if this rate is at literature minimum or can be reduced
- This reaction is 91% of C2 production - reducing it would directly impact the loop

**Pros:** Directly attacks root cause
**Cons:** Might be well-constrained in literature

### Option 2: Increase C2 Wall Losses (Indirect)
**Current:** C2 wall loss is 15% of C2 destruction

**Action:**
- Increase C2 wall loss coefficient (currently in range [1e-4, 2e3])
- Prevent C2 from accumulating and converting to CH

**Pros:** Wall losses have huge tuning range (20 million x!)
**Cons:** Physical parameter, should match experimental geometry

### Option 3: Reduce C2 + H → CH Conversion (Most Direct!)
**Target reaction:** C2 + H → CH + C

**Current rate:** 2.42×10¹⁵ cm⁻³ s⁻¹ (59% of CH production!)

**Action:**
- Check this rate against literature
- If it can be reduced, this directly cuts off CH production

**Pros:** This is THE dominant CH production pathway
**Cons:** Need to check literature constraints

### Option 4: Increase CH + CH → C2 Conversion (Use CH to make C2)
**Target reaction:** CH + CH → C2 + H2

**Current rate:** Only 0.01% of C2 production!

**Potential:** This reaction could consume excess CH while producing useful C2

**Action:**
- Check if this rate can be increased (was 2.16×10⁻¹⁰, above lit max of 1.8×10⁻¹⁰)
- Might help but CH + CH reaction rate is limited by CH² dependence

### Option 5: Multi-Species Optimization
**Targets:** H, CH, C2, **AND C2H2**

**Action:**
- Add C2H2 as a 4th target species
- Set target C2H2 < 1×10¹³ cm⁻³ (reduce by factor of 1.04)
- This forces optimizer to address root cause

**Pros:** Systematic, addresses cascade effect
**Cons:** Need to know what C2H2 should be

---

## Immediate Next Steps

1. **Check literature rates for key reactions:**
   - C2H2 + H → C2 + H2 + H (91% of C2 production)
   - C2 + H → CH + C (59% of CH production)
   - CH + CH → C2 + H2 (could consume CH)

2. **Run diagnostic with optimized parameters from step 4** (f(x)=1097):
   - Extract parameters from optimization
   - See what changed
   - Check if C2H2 → C2 → CH loop was broken

3. **Consider adding C2H2 as optimization target:**
   - Current: 1.04×10¹³ cm⁻³
   - Reasonable target: ~5×10¹² cm⁻³ (factor of 2 reduction)

---

## Conclusion

**The CH problem is NOT primarily from electron-impact chemistry.**

It's a **downstream cascade** from high C2H2:
- C2H2 dissociates to make C2 (91% of C2 comes from this)
- C2 reacts with H to make CH (59% of CH comes from this)
- CH accumulates because these production rates overwhelm wall losses

**To fix CH, we must either:**
1. Reduce C2H2 production/accumulation (upstream fix)
2. Reduce C2H2 → C2 conversion (block the cascade)
3. Reduce C2 → CH conversion (block final step)
4. Increase CH wall losses beyond what we've tried

The 2-parameter optimization failed because it targeted Ne, but electron-impact is only 3.5% of CH production. The real source is the C2 feedback loop.
