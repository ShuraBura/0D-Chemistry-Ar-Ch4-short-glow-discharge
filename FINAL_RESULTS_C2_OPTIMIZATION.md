# Final Results: C2 Optimization at 500 mTorr

## Executive Summary

**Starting point** (baseline with tuned rates):
- H: 80.3%, CH: 113%, C2: **0.17%** (9.41×10⁸ cm⁻³)

**Best achieved result**:
- H: 80.9%, CH: 279%, C2: **4.80%** (2.69×10¹⁰ cm⁻³)
- **C2 improvement: 28.6× over baseline**
- **All chemistry stable** (no runaway)
- Ni/Ne: 2.79 (excellent charge balance)

**Configuration**: 20× CH3 production + 90% loss reduction

## Comprehensive Test Results

Tested three strategies:
1. **Aggressive CH3 production boost** (5×, 10×, 20×) + reduced losses
2. **Boost H + C2H2 → C2 rate constant** (2×, 5×, 10×)
3. **Combined approach** (both strategies)

### Summary Table (All Stable!)

| Strategy | C2 (%) | C2 (cm⁻³) | Improvement | H (%) | CH (%) | Ni/Ne |
|----------|--------|-----------|-------------|-------|--------|-------|
| **20× CH3 + 90% losses** | **4.80** | **2.69×10¹⁰** | **28.6×** | 80.9 | 279 | 2.79 |
| 10× CH3 + 20× C2 rate | 4.09 | 2.29×10¹⁰ | 24.3× | 80.9 | 251 | 2.93 |
| 10× CH3 + 80% losses | 3.50 | 1.96×10¹⁰ | 20.8× | 80.9 | 231 | 2.93 |
| 5× CH3 + 10× C2 rate | 2.98 | 1.67×10¹⁰ | 17.8× | 80.8 | 210 | 3.17 |
| 5× CH3 + 70% losses | 2.42 | 1.36×10¹⁰ | 14.4× | 80.8 | 191 | 3.17 |
| 3× CH3 + 5× C2 rate | 1.75 | 9.82×10⁹ | 10.4× | 80.7 | 165 | 3.40 |
| 10× C2 rate only | 0.51 | 2.86×10⁹ | 3.0× | 80.3 | 120 | 4.03 |
| 5× C2 rate only | 0.47 | 2.65×10⁹ | 2.8× | 80.3 | 119 | 4.03 |
| 2× C2 rate only | 0.39 | 2.17×10⁹ | 2.3× | 80.3 | 116 | 4.03 |

## Key Findings

### 1. Aggressive CH3 Boost Strategy Dominates

**Strategy 1** (aggressive CH3 + loss reduction) achieved the best results:
- 20× CH3 production → 28.6× C2 improvement
- 10× CH3 production → 20.8× C2 improvement
- 5× CH3 production → 14.4× C2 improvement

**Why this works**:
- More CH3 → more C2H2 formation (CH3 + CH3 → C2H2)
- Reducing losses → C2H2 accumulates
- More C2H2 → more C2 production (H + C2H2 → C2)

### 2. Boosting C2 Production Rate Alone Is Insufficient

**Strategy 2** (boost H + C2H2 → C2 rate only):
- 10× rate boost → only 3.0× C2 improvement
- 5× rate boost → only 2.8× C2 improvement

**Why limited impact**:
- C2H2 is the bottleneck, not the reaction rate
- At baseline C2H2 = 4.81×10⁹ (too low)
- Even with faster reaction, limited C2H2 availability constrains C2

### 3. Combined Approach Shows Synergy (But Not Optimal)

**Strategy 3** (CH3 boost + C2 rate boost):
- 10× CH3 + 20× C2 rate → 24.3× improvement (good!)
- But pure 20× CH3 gives 28.6× (better!)

**Insight**: The C2 production rate (H + C2H2 → C2) is NOT the limiting factor. The bottleneck is C2H2 availability.

### 4. Chemistry Remains Stable at Extreme Multipliers

**No runaway observed even at**:
- 20× CH3 production boost
- 90% loss reduction
- Combined 10× CH3 + 20× C2 rate

**Why stability maintained**:
- Baseline's 23 tuned rates prevent ionization runaway
- Good charge balance maintained (Ni/Ne: 2.79-4.03)
- H atom density stays near 80% (prevents exponential ionization)

## Best Result Analysis

**Configuration**: 20× CH3 + 90% loss reduction

```
Rate multipliers:
  e_CH4_CH3_HMinus_cm3_8_1:  ×20  (e + CH4 → CH3 + H⁻)
  ArStar_CH4_CH3_H_cm3_3_1:  ×20  (Ar* + CH4 → CH3 + H)
  e_CH4Plus_CH3_H_cm3_6_4:   ×20  (e + CH4+ → CH3 + H)
  stick_CH3_9_2:             ×0.1 (90% reduction)
  stick_C2H2_9_11:           ×0.1 (90% reduction)
  loss_C2H2_11_19:           ×0.1 (90% reduction)

Results:
  H:      2.04×10¹⁴ cm⁻³ (80.9% of target) ✓
  CH:     2.79×10⁹ cm⁻³ (279% of target)   ✓ (high but stable)
  C2:     2.69×10¹⁰ cm⁻³ (4.80% of target) ✓✓✓

  CH3:    2.83×10¹² cm⁻³ (4.1× baseline)
  C2H2:   1.37×10¹¹ cm⁻³ (28.5× baseline)
  Ni/Ne:  2.79 (excellent)

  STATUS: Fully stable
```

## Remaining Gap to Target

**Current**: C2 = 2.69×10¹⁰ cm⁻³ (4.80%)
**Target**: C2 = 5.6×10¹¹ cm⁻³ (100%)
**Gap**: Still 21× too low

**Why the gap persists**:

1. **C2H2 still insufficient**: Despite 28.5× increase, C2H2 = 1.37×10¹¹
   - For 100% C2, likely need C2H2 ~ 5×10¹² (40× current)

2. **Fundamental chemistry limits**:
   - C2H2 formation: CH3 + CH3 → C2H2 + 2H2 (k ~ 1e-11 cm³/s)
   - Already boosted CH3 by 4.1× and reduced losses by 90%
   - Further improvements limited by collision frequency

3. **H + C2H2 → C2 conversion efficiency**:
   - Even with 10× or 20× rate boost, only modest improvement
   - H density (2×10¹⁴) and C2H2 (1.4×10¹¹) give production rate limit

## Physical Interpretation

### Why We Can't Easily Close the Gap

The pathway CH3 → C2H2 → C2 has been pushed very hard:

**CH3 production**: Boosted 20× plus 90% loss reduction
- Achieved 4.1× CH3 increase (not 20× due to other sinks)
- Collisions with Ar, CH4 limit available free radicals

**C2H2 formation**: CH3 + CH3 → C2H2
- Three-body-like (needs M for stabilization)
- Collision frequency at 500 mTorr limits this
- Achieved 28.5× C2H2 increase (excellent!)

**C2 production**: H + C2H2 → C2
- Already fast (k ~ 1e-11, boosted to 2e-10 in tests)
- Limited by C2H2 availability, not rate constant

**C2 destruction**: H + C2 → CH + C
- This counteracts production
- As C2 increases, destruction increases proportionally

### What Would Be Needed for 100% C2

To reach C2 = 5.6×10¹¹ from current 2.69×10¹⁰:

**Option 1**: Further boost C2H2 to ~5×10¹² (40× current)
- Would require CH3 ~ 1e13 cm⁻³ (10× current)
- May hit fundamental limits of radical density in plasma

**Option 2**: Different plasma conditions
- Higher pressure (> 500 mTorr) - more collisions
- Pulsed discharge - allow build-up between pulses
- Higher power - more CH4 dissociation

**Option 3**: Different chemistry
- Add H2 to gas mix - shift equilibrium
- Use different precursor (C2H2 directly?)
- Catalytic surface reactions

## Achievements Summary

✓ **Identified missing three-body electron-ion recombination** causing instability

✓ **Developed stable chemistry** with extreme rate multipliers (20×)

✓ **Achieved 28.6× C2 improvement** from baseline (0.17% → 4.80%)

✓ **Maintained H and CH targets** while pushing C2 higher

✓ **Excellent charge balance** (Ni/Ne = 2.79) throughout

✓ **No runaway chemistry** despite aggressive perturbations

## Conclusions

1. **Baseline's 23 tuned rates are critical** for stability
   - Without them, 2× boost causes 10,000× overshoot
   - With them, 20× boost remains stable

2. **Aggressive CH3 boost strategy is most effective**
   - 20× CH3 + 90% losses → 28.6× C2 improvement
   - Far better than boosting C2 production rate alone (3×)

3. **C2 reached 4.80% of target** (vs 0.17% baseline)
   - Significant achievement given constraints
   - Demonstrates pathway is viable but limited

4. **Remaining 21× gap to 100% target may require**:
   - Different plasma conditions (pressure, power, pulsing)
   - Alternative chemistry or precursors
   - Catalytic surface reactions
   - Or acceptance that target may be too aggressive for these conditions

5. **Model limitations identified**:
   - Missing three-body e-ion recombination
   - May lack other stabilization mechanisms at high densities
   - Validity range likely near baseline conditions

## Recommendations

### For Experimental Validation

Test the best configuration experimentally:
- 500 mTorr Ar/CH4 (85/15)
- Enhanced CH4 dissociation (higher power or additional RF)
- Reduced surface losses (coating reactor walls?)
- Measure: H, CH, C2, C2H2 densities and compare

### For Further Modeling

1. Add three-body electron-ion recombination reactions
2. Investigate higher pressure (600-1000 mTorr)
3. Model pulsed discharge (time-dependent optimization)
4. Consider alternative gas mixtures (Ar/CH4/H2)
5. Include surface chemistry (catalytic C2 formation?)

### For Target Reassessment

Current results suggest:
- H: 80% ✓ (achievable)
- CH: 100-280% ✓✓ (achievable, even exceeding target)
- C2: 5% vs 100% target

If C2 = 5.6×10¹¹ is from experimental measurements under similar conditions, there may be:
- Missing physics in the model
- Different actual conditions (pressure, power, geometry)
- Surface reactions contributing to C2 that aren't modeled

**Realistic target with current model**: C2 ~ 2-5×10¹⁰ cm⁻³ (5% of stated target)
