# CH Density Problem Analysis Report

## Executive Summary

**Problem:** CH density is **53x too high** (5.3×10¹⁰ cm⁻³ vs target 1.0×10⁹ cm⁻³)

**Critical Finding:** The MATLAB code was already heavily tuned outside literature bounds to try to reduce CH, but it's still 53x too high. This indicates a **fundamental chemistry issue**, not just parameter tuning.

---

## Current Status vs Targets

| Species | Target | Current | Status |
|---------|--------|---------|--------|
| H       | 5.18×10¹³ | 3.48×10¹³ | 67% of target (need +49%) |
| **CH**  | **1.0×10⁹** | **5.3×10¹⁰** | **53x TOO HIGH** |
| C2      | 1.3×10¹¹ | 7.59×10¹¹ | 5.8x too high |

---

## Critical Findings

### 1. CH_CH_C2_H2_cm3_5_4: CH + CH → C2 + H2

**THE KEY REACTION** - Currently **OUTSIDE literature bounds!**

```
Reaction:      CH + CH → C2 + H2
Current value: 2.16×10⁻¹⁰ cm³/s
Literature:    [1.20×10⁻¹⁰, 1.80×10⁻¹⁰] (Baulch et al. 2005)
VIOLATION:     20% ABOVE maximum!
```

**Why this matters:**
- This reaction is **quadratic in [CH]** → Very powerful for high CH densities
- Removes 2 CH molecules, produces 1 C2
- Already pushed 20% above literature to consume more CH
- This is the **primary C2 source** (explains why C2 is also high)

### 2. Twenty-Nine Rates Outside Literature Bounds

The MATLAB code has **29 rates violating literature constraints**, including:

**Electron-ion recombination (25-50% above max):**
- ArPlus_e_Ar: 25% above
- CH3Plus_e_CH3: 25% above
- CH5Plus_e_CH4_H: 25% above
- HMinus_ArPlus_H_Ar: **50% above**

**Ion-neutral reactions (20% violations):**
- CH3Plus_CH4_CH5Plus_CH2: 20% below
- ArPlus_CH4_Ar_CH4Plus: 20% below

### 3. Fifty-Seven CH Rates Already at Min/Max Boundaries

**CH Production Rates (at MINIMUM to reduce CH production):**
```
e_CH4_CH_H_H2_cm3_1_11:        AT MIN (3.0x range available)
ArStar_CH2_CH_H_cm3_3_4:       AT MIN (1.5x range)
CH2_H_CH_H2_cm3_7_1:           AT MIN (2.2x range)
```

**CH Loss Rates (at MAXIMUM to increase CH loss):**
```
loss_CH_11_9 (wall loss):      AT MAX (10x range!) ***CRITICAL***
stick_CH_9_3 (wall stick):     AT MAX (5x range)
CH_H_C_H2_cm3_7_3:            AT MAX (1.5x range)
```

**Interpretation:** The code was already aggressively tuned to minimize CH production and maximize CH loss, yet CH is STILL 53x too high.

---

## Key CH Reactions Analysis

### CH Production Pathways

| Reaction | Rate Constant | Range | Status | Source |
|----------|---------------|-------|--------|--------|
| e + CH4 → CH + H2 + H | 3.0×10⁻¹¹ | [2×10⁻¹¹, 1×10⁻¹⁰] | Middle | Janev & Reiter (2002) |
| e + CH4 → CH + H + H2 | 2.0×10⁻¹¹ | [2×10⁻¹¹, 6×10⁻¹¹] | **MIN** | Janev & Reiter (2002) |
| Ar* + CH2 → CH + H + Ar | 8.0×10⁻¹¹ | [8×10⁻¹¹, 1.2×10⁻¹⁰] | **MIN** | Phelps (1999) |
| CH2 + H → CH + H2 | 1.0×10⁻¹¹ | [1×10⁻¹¹, 2.25×10⁻¹¹] | **MIN** | Baulch et al. (2005) |

**Primary CH source:** Electron impact on CH4 → Depends on [e] and [CH4]

### CH Loss Pathways

| Reaction | Rate Constant | Range | Status | Source |
|----------|---------------|-------|--------|--------|
| **CH + CH → C2 + H2** | **2.16×10⁻¹⁰** | **[1.2×10⁻¹⁰, 1.8×10⁻¹⁰]** | **20% ABOVE MAX!** | **Baulch et al. (2005)** |
| CH + CH3 → C2H2 + H2 | 1.0×10⁻¹⁰ | [8×10⁻¹¹, 1.2×10⁻¹⁰] | Middle | Baulch et al. (2005) |
| CH + CH3 → C2H3 + H | 8.0×10⁻¹¹ | [8×10⁻¹¹, 1.2×10⁻¹⁰] | **MIN** | Baulch et al. (2005) |
| CH + H → C + H2 | 1.2×10⁻¹⁰ | [8×10⁻¹¹, 1.2×10⁻¹⁰] | **MAX** | Baulch et al. (2005) |
| CH → wall | 1.0×10⁴ s⁻¹ | [1×10³, 1×10⁴] | **MAX (10x range!)** | Estimated |
| CH wall sticking | 6.25×10³ s⁻¹ | [1.25×10³, 6.25×10³] | **MAX (5x range)** | Jauberteau et al. (1998) |

---

## Why CH is So High: Root Cause Analysis

### Hypothesis 1: Electron Density Too High

If [e] is higher than experiment, then:
- e + CH4 → CH production is artificially high
- This would also affect ionization rates

**Need to check:** What is [e] in simulation vs experiment?

### Hypothesis 2: CH Wall Loss Model Incorrect

Current model uses:
```
loss_CH_11_9 = 1×10⁴ s⁻¹ (at MAX of [1×10³, 1×10⁴])
```

But this is a **crude approximation** for wall loss. The real wall loss should depend on:
- Diffusion to walls (geometry-dependent)
- Sticking probability
- Surface reactions

**Problem:** The loss rate is linear in [CH], but at high [CH] densities, walls might saturate.

### Hypothesis 3: Missing CH Loss Pathway

Are there CH loss reactions missing from the chemistry? Common ones:
- CH + O → CO + H (but no oxygen in system)
- CH + O2 → HCO + O (but no oxygen)
- CH + other radicals?

**Need to check:** Is the CH chemistry complete for Ar/CH4 systems?

### Hypothesis 4: CH + CH → C2 Rate Literature Value Wrong

The Baulch et al. (2005) range [1.2×10⁻¹⁰, 1.8×10⁻¹⁰] might be:
- For a different temperature (current sim: 400K)
- For a different pressure regime
- Updated by more recent studies

**Need to:** Look up original Baulch et al. (2005) paper and check for newer studies.

---

## Strategy Options

### Option A: Bring All Rates Within Literature Bounds (Conservative)

**Action:**
1. Fix CH_CH_C2_H2_cm3_5_4 from 2.16×10⁻¹⁰ → 1.8×10⁻¹⁰ (literature max)
2. Fix all 29 violations to be within bounds
3. Run simulation and see impact

**Expected Impact:**
- Reducing CH + CH → C2 by 1.2x will **increase** CH (bad!)
- But might make C2 more reasonable (currently 5.8x too high)
- Uncertain overall effect

**Pros:** Scientifically rigorous, defensible
**Cons:** Might make CH problem worse

### Option B: Investigate Electron Density (Physical)

**Action:**
1. Check [e] in simulation output
2. Compare to experimental measurements (if available)
3. If [e] is too high, investigate why:
   - Ionization rates too high?
   - Recombination rates too low?
   - Boundary conditions wrong?

**Expected Impact:**
- Reducing [e] would reduce all electron-impact reactions
- Would directly reduce CH production from e + CH4 → CH
- Would also affect H production and ion chemistry

**Pros:** Addresses root cause if [e] is the problem
**Cons:** Might break other species agreements

### Option C: Re-examine Wall Loss Model (Geometry)

**Action:**
1. Check experimental geometry (L, R values)
2. Calculate proper diffusion-limited wall loss
3. Compare to current loss_CH_11_9 = 1×10⁴ s⁻¹
4. Check if wall loss should be higher

**Formula for diffusion loss:**
```
ν_wall = D_CH * (2.405/R)² + D_CH * (π/L)²
```

Where D_CH is CH diffusion coefficient in Ar at 400K.

**Expected Impact:**
- If calculated loss is higher, CH density would drop
- But this only helps if diffusion model is more accurate

**Pros:** Physically motivated correction
**Cons:** Requires additional parameters (D_CH)

### Option D: Constrained Multi-Parameter Optimization (Brute Force)

**Action:**
1. Fix the 29 out-of-bounds rates to be within literature
2. Select 10-20 key rates with largest ranges
3. Run global optimization (differential evolution) with:
   - Objective: Match H, CH, C2 targets
   - Constraints: All rates within literature bounds
   - Weighted error (CH weighted 10x)

**Expected Impact:**
- Will find best possible solution within literature constraints
- Might reveal if targets are achievable or contradictory

**Pros:** Systematic, explores full parameter space
**Cons:** Computationally expensive, might fail if targets are impossible

### Option E: Literature Survey for CH + CH → C2 + H2 (Rigorous)

**Action:**
1. Look up Baulch et al. (2005) original paper
2. Check temperature dependence: k(400K) vs k(298K)
3. Search for newer studies (post-2005)
4. Check if other groups have measured this for Ar/CH4 plasmas

**Expected Impact:**
- Might find that range is different at 400K
- Might find newer, more accurate measurements
- Could justify using a different value

**Pros:** Most scientifically rigorous
**Cons:** Time-consuming, might not change anything

---

## Recommended Strategy

I recommend a **hybrid approach** combining multiple options:

### Phase 1: Quick Checks (Physical)

1. **Check electron density** in simulation output
   - Compare [e] vs experimental measurement
   - If too high, investigate ionization/recombination balance

2. **Verify geometry parameters** (L, R)
   - Ensure they match experimental setup
   - Recalculate wall loss rates if needed

3. **Check all species densities** (not just H, CH, C2)
   - Look for other discrepancies
   - Might reveal systematic issues

### Phase 2: Literature Verification (Rigorous)

4. **Look up CH + CH → C2 + H2 rate** in original sources
   - Temperature dependence
   - Pressure dependence
   - Any updates since 2005?

5. **Verify other key rates** against literature:
   - e + CH4 → CH + H2 + H
   - CH wall loss mechanisms

### Phase 3: Constrained Optimization (If Needed)

6. **Bring all rates within literature bounds** first

7. **Run constrained optimization** with:
   - 10-20 rates with largest ranges
   - All rates constrained to literature
   - Weighted objective (CH errors weighted 10x)

---

## Questions for You

Before proceeding, I'd like your input on:

1. **Do you have experimental [e] measurements?**
   - If electron density is too high in sim, that's likely the root cause

2. **What is the experimental geometry?**
   - Length L, radius R
   - This affects all wall loss calculations

3. **Priority: Scientific rigor vs fitting targets?**
   - Option A: Stay strictly within literature bounds (might not hit targets)
   - Option D: Allow optimization within bounds (best fit possible)

4. **Have you checked other species besides H, CH, C2?**
   - e.g., CH3, CH4, C2H2, ions
   - Might reveal systematic issues

5. **Do you want me to look up literature right now?**
   - I can search for CH + CH → C2 + H2 studies
   - Check for temperature/pressure dependencies

---

## Next Steps

Based on your answers, I'll proceed with the most appropriate strategy. My current recommendation:

1. **Immediate:** Check [e] density from last simulation output
2. **Quick:** Verify geometry and recalculate wall losses
3. **If still stuck:** Literature survey for CH + CH rate
4. **Final approach:** Constrained multi-parameter optimization

What would you like me to focus on first?
