# Critical Finding: Chemistry Model Instability

## Summary

Attempting to improve C2 from the baseline (16.8% → 100% of target) reveals a **fundamental instability** in the chemistry model. Even modest 2× boosts to reaction rates cause **10,000× overshoots** in densities, rendering the chemistry completely unphysical.

## The Multi-Stage Strategy

### Stage 1: Analyze H Impact on C2

**Goal**: Determine if reducing H would help C2 (since H + C2 → CH accounts for 72.8% of C2 loss)

**Result**: **H HELPS C2, doesn't hurt it!**

Key finding from `analyze_H_vs_C2_tradeoff.py`:
```
Rate constants:
  k(H + C2H2 → C2): 5.00e-10 cm³/s  (C2 production)
  k(H + C2 → CH):   3.53e-11 cm³/s  (C2 destruction)
  Ratio: 14.2× faster production than destruction

If we boost C2H2 from 4.81e9 to 5e10 (only 10× increase):
  At baseline H: C2 reaches 92% of target!
```

**Conclusion**: Don't reduce H. Focus on boosting C2H2.

### Stage 2: Boost C2H2 Production

**Goal**: Increase C2H2 from 4.81×10⁹ to ~5×10¹⁰ (10× increase)

**Strategy**: Since C2H2 ∝ [CH3]², boost CH3 production:
- Current CH3: 6.85×10¹¹
- Target CH3: ~2×10¹² (3× increase → 9× C2H2 increase)

**CH3 production pathways**:
1. e + CH4 → CH3 + H⁻  (rate: 7.42e-13 cm³/s in baseline)
2. Ar* + CH4 → CH3 + H  (rate: 6.92e-10 cm³/s in baseline)
3. e + CH4⁺ → CH3 + H   (rate: 6.44e-07 cm³/s in baseline)

### Stage 3: Systematic Testing

Created `test_C2H2_boost_systematic.py` to test specific multipliers starting from baseline.

## The Critical Problem

**EVERY test shows massive instability:**

### Test 1: Boost e + CH4 → CH3 by 2×
```
Input:  2× multiplier on single rate
Output:
  H:      83,777% of target  (1,048× baseline!)
  CH:     1,196% of target   (11.9× baseline)
  C2:     5,248% of target   (312× baseline!)
  CH3:    9,498× baseline    (completely unphysical)
  Ni/Ne:  3,700,000          (vs 3.12 baseline!)
```

### Test 2: Boost Ar* + CH4 → CH3 by 2×
```
Output:
  H:      92,787% of target
  CH:     1,331% of target
  C2:     5,393% of target
  CH3:    10,239× baseline
  Ni/Ne:  3,359,000
```

### Test 3-7: All Other Combinations
Similar behavior - any perturbation triggers runaway chemistry with 1000×-10,000× overshoots.

## Analysis of the Instability

### What's Happening?

1. **Tiny rate = 7.42e-13**: The baseline `e + CH4 → CH3` rate is extraordinarily small
2. **No negative feedback**: When CH3 increases, there's no proportional loss mechanism
3. **Exponential growth**: CH3 → C2H2 → more reactions → more CH3 (positive feedback loop)
4. **Runaway ionization**: More chemistry → more ionization → Ni/Ne goes to millions

### Why the Baseline Works

The baseline achieves H=79.9%, CH=100.7%, C2=16.8% with:
```
Te: 1.31 eV  (very low - limits ionization)
Ne: 1.22e8   (carefully tuned)
E:  250 V/cm
```

**The baseline sits at a knife-edge stability point.**

Any change triggers instability because:
- Rates are finely balanced
- Missing or inadequate loss mechanisms
- Positive feedback loops dominate

### Evidence This is Fundamental

1. **Global optimizer failed** (`optimize_from_baseline.py`):
   - With 23 tunable rates, optimizer made C2 WORSE (0.16% vs 16.8%)
   - Easier to match H/CH by increasing C2 losses than fixing chemistry

2. **Focused optimizer failed** (`optimize_C2H2_boost.py`):
   - With 9 rates (CH3/C2H2 only), produced:
     - H: 324,562% of target
     - CH3: 4.71e15 (6,878× baseline)
     - Ni/Ne: 5.6 million

3. **Systematic tests failed** (`test_C2H2_boost_systematic.py`):
   - Even 2× single rate boost → 10,000× overshoot
   - ALL 7 test scenarios show same problem
   - No way to "gently" nudge chemistry higher

## Possible Explanations

### Hypothesis 1: Artificially Tuned Rates
The baseline rate constants may not be purely physical but include artificial constraints/tuning to achieve stability. This would explain:
- Why rates are at such specific values
- Why any change breaks stability
- Why `e_CH4_CH3_HMinus` = 7.42e-13 is so small

### Hypothesis 2: Missing Physics
The model may be missing critical loss/quenching mechanisms that provide negative feedback at higher densities:
- Three-body recombination
- Collisional quenching
- Plasma-wall interactions
- Self-absorption of radiation

### Hypothesis 3: Model Validity Boundary
The chemistry model may only be valid near certain plasma conditions:
- Te ~ 1.3 eV
- Ne ~ 1e8 cm⁻³
- Low ionization fraction

Outside this regime, the ODE system becomes stiff/unstable.

## Implications

### What We've Proven

✓ **H is not the problem** - higher H helps C2 production
✓ **C2H2 is the bottleneck** - need 1000× increase
✓ **Baseline is stable** - H=79.9%, CH=100.7%, C2=16.8%
✗ **Cannot push higher** - any perturbation → instability

### The Hard Truth

**The baseline (C2 = 16.8%) may be the BEST achievable result with this chemistry model.**

Pushing C2 higher requires:
1. Boosting CH3 (to increase C2H2)
2. But boosting CH3 triggers runaway chemistry
3. No "gentle" path exists

### Options Going Forward

1. **Accept baseline as optimal**
   - H: 79.9% ✓
   - CH: 100.7% ✓✓✓
   - C2: 16.8% (not ideal but stable)
   - Ni/Ne: 3.12 ✓

2. **Revise chemistry model**
   - Add missing loss mechanisms
   - Include three-body reactions
   - Model plasma-wall sheaths properly
   - Validate at higher densities

3. **Reconsider targets**
   - Are all three targets achievable simultaneously?
   - Can C2 target be reduced?
   - Is there experimental evidence for target densities?

4. **Try different plasma conditions**
   - Higher pressure (> 500 mTorr)
   - Pulsed discharge (time-varying conditions)
   - Different gas mixture (vary Ar/CH4 ratio)

## Recommendation

Before continuing optimization, we need to:

1. **Validate the chemistry model** at conditions away from baseline
2. **Understand why baseline rates are what they are** - physical or fitted?
3. **Identify missing physics** that could provide stability at higher densities
4. **Set realistic expectations** - can all 3 targets be achieved simultaneously?

The current model appears to hit a **fundamental limit at C2 ~ 17% of target**. Breaking through requires understanding the root cause of the instability.
