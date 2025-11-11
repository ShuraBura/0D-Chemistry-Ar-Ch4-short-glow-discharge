# Fundamental Constraint Analysis: C2H2 vs Charge Balance

## Summary

We have discovered a **fundamental trade-off** in the Ar/CH4 plasma model:

**High C2H2 (>1e12 cm⁻³) is incompatible with stable charge balance (Ni/Ne = 2-7).**

## Mechanism

The coupling mechanism is:

```
Ar* metastables + CH4 → TWO pathways:
   1. Ar* + CH4 → CH3 + H + Ar    (dissociation - produces CH3 for C2H2)
   2. Ar* + CH4 → CH4+ + e + Ar   (Penning ionization - creates ions)
```

You **cannot selectively boost** pathway (1) without also boosting pathway (2).

Both reactions use the same Ar* population, so:
- High Ar* → High CH3 → High C2H2 ✓
- High Ar* → High Penning ionization → High Ni/Ne ✗

## Optimization Results

### Approach 1: Fixed Ionization Rates (Correct Physics)
- **Ionization**: Fixed at literature values (not tunable)
- **Result**: Ni/Ne = 2.95 ✓, C2H2 = 3.54e+09 (0.35% of 1e12 target)
- **Ar***: 2.58e+08 cm⁻³ (limited by charge balance constraint)

### Approach 2: Boosted Metastables (Attempted Fix) - WITH BUG
- **Bug**: Used wrong rate name ('e_Ar_ArStar_cm3_1_4' instead of 'e_Ar_ArStar_cm3_1_7')
- **Result**: Optimizer couldn't tune Ar* excitation, converged to similar result as Approach 1
- **Lesson**: Always verify rate constant names match the model!

### Approach 3: Boosted Metastables (FIXED) - Revealing Trade-off

#### Best results explored by optimizer:

| Eval | Ar* (cm⁻³) | C2H2 (cm⁻³) | Ni/Ne | Status |
|------|-----------|------------|-------|--------|
| 2 | 1.62e+12 | **6.54e+13** (6.5× target!) | 1322 | ✗ Extreme charge imbalance |
| 15 | 2.62e+11 | 6.44e+12 (6.4× target) | 90 | ✗ Still unstable |
| 32 | 7.15e+10 | 1.30e+12 (1.3× target) | 50.9 | ✗ Unstable |
| 53 | 5.81e+10 | 2.51e+12 (2.5× target) | 36.3 | ✗ Unstable |
| 86 | 2.88e+10 | 7.76e+10 | 24.3 | ✗ Still too high |
| **128** | **4.17e+09** | **2.69e+10** | **2.94** | **✓ STABLE!** |

#### Convergence trajectory:
- **Ni/Ne**: 1322 → 90 → 50.9 → 36.3 → 24.3 → **2.94** ✓
- **C2H2**: 6.54e+13 → 6.44e+12 → 1.30e+12 → 2.51e+12 → 7.76e+10 → **2.69e+10**

The optimizer clearly shows:
- To enforce Ni/Ne < 7, it must reduce Ar* to ~4×10⁹
- At this Ar* level, C2H2 can only reach ~3×10¹⁰ (0.03× target)
- To reach C2H2 > 1×10¹² requires Ar* > 5×10¹⁰
- But Ar* > 5×10¹⁰ → Ni/Ne > 30 (unstable)

### Approach 4: Unphysical Ionization (What Previous Optimizer Did)
- **Method**: Boosted ionization rates 10× beyond literature values
- **Result**: Ni/Ne = 215, C2H2 = 4.09e+12 (4× target!)
- **Problem**: Violates physics, creates charge imbalance
- **Why it worked**: High ionization → High Ar* → High CH3 → High C2H2
- **Why it's wrong**: Real experiments cannot sustain Ni/Ne = 215

## The Constraint Visualized

```
Ar* Level vs Outcomes:

Ar* ~ 1e12:  C2H2 = 6e13 ✓ (6× target)   Ni/Ne = 1322 ✗ (200× too high)
Ar* ~ 1e11:  C2H2 = 6e12 ✓ (6× target)   Ni/Ne = 90   ✗ (13× too high)
Ar* ~ 5e10:  C2H2 = 2e12 ✓ (2× target)   Ni/Ne = 36   ✗ (5× too high)
Ar* ~ 3e10:  C2H2 = 8e10                 Ni/Ne = 24   ✗ (3.4× too high)
Ar* ~ 4e09:  C2H2 = 3e10                 Ni/Ne = 2.94 ✓ STABLE

Constraint:
  To get Ni/Ne < 7:  Need Ar* < 1e10
  To get C2H2 > 1e12: Need Ar* > 5e10

  THESE ARE INCOMPATIBLE!
```

## Implications

1. **The model cannot simultaneously achieve:**
   - C2H2 ≥ 1×10¹² cm⁻³
   - Ni/Ne in range [2, 7] (stable plasma)

2. **Three possible explanations:**

   a) **The model is missing physics** (most likely)
      - Perhaps missing a C2H2 production pathway that doesn't require high Ar*
      - Perhaps missing a way to produce CH3 without Ar* (direct e + CH4 → CH3?)
      - Perhaps Penning ionization rate is overestimated
      - Perhaps there's a CH3 production mechanism we're not including

   b) **The experimental targets are inconsistent**
      - Experiment achieves C2H2 > 1e12 with stable plasma
      - But our model says this requires unstable plasma
      - This suggests the model is incomplete

   c) **Parameter space issue**
      - Perhaps there's a narrow window of conditions we haven't explored
      - But optimizer explored 128+ combinations and found clear trend
      - Seems unlikely there's a "hidden" solution

3. **Best achievable with current model:**
   - **With stable plasma (Ni/Ne ~ 3)**: C2H2 ~ 3×10¹⁰ cm⁻³
   - **With C2H2 ≥ 1×10¹²**: Ni/Ne > 30 (unstable)

## Recommendations

1. **Review the model for missing reactions:**
   - Direct electron-impact dissociation: e + CH4 → CH3 + H + e
   - Alternative CH3 sources that don't require Ar*
   - Check if Penning dissociation/ionization branching ratios are correct
   - Verify that Ar* + CH4 channels are correctly implemented

2. **Re-examine experimental conditions:**
   - What is the actual Ar* density in the experiment?
   - How is charge balance maintained at high C2H2?
   - Are there pulsed or time-varying effects not captured in 0D steady-state model?

3. **Consider alternative model structures:**
   - Time-dependent (non-steady-state) behavior
   - Spatial gradients (move to 1D model?)
   - Different gas mixtures or additives

## Files

- `optimize_fixed_ionization_rates.py` - Approach 1 (stable but low C2H2)
- `optimize_boosted_metastables.py` - Approach 3 (reveals trade-off)
- `why_c2h2_needs_unstable_plasma.py` - Analysis showing Ar* is bottleneck
- `analyze_boosted_metastables_result.py` - Comparison of all approaches
- `optimization_results_fixed_ionization/best_f22.8.json` - Stable result
- `optimization_results_boosted_metastables/best_f45.9.json` - Best stable result with boosting

## Conclusion

The optimization campaign has been extremely valuable - it revealed a **fundamental physical constraint** in the current model rather than just a tuning problem. The model appears to be missing some physics that allows the experiment to achieve high C2H2 with stable charge balance.

**Next step**: Investigate missing reaction pathways or reconsider model assumptions.
