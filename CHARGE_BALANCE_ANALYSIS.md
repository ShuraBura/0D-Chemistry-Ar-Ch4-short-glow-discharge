# CHARGE BALANCE CRISIS - ROOT CAUSE ANALYSIS

## Executive Summary

V5 optimization achieved H=0.40x, CH=5.49x, C2=0.37x but has a **CATASTROPHIC charge balance violation**:
- **Ion/electron ratio: 0.524x** (target: 4.5x) 
- **88.4% deficit in positive ions**
- Missing 5.69×10⁹ cm⁻³ of positive ions

## Root Cause: Missing Ion-Electron Recombination Chemistry

### 1. Electron Production/Loss Imbalance

```
Electron production (ionization): 5.861×10¹⁴ cm⁻³/s
Electron loss (all mechanisms):   1.287×10¹³ cm⁻³/s
NET IMBALANCE:                    5.732×10¹⁴ cm⁻³/s (97.8%!)
```

**Why**: Electrons are produced 45x faster than they're consumed!

### 2. Ion Loss Mechanisms

Top ion loss pathways:
- **CH5+**: 98.6% lost via drift (ambipolar diffusion)
- **CH3+**: 40.9% → C2H5+ + H2, 39.3% → CH5+ + CH2, 19.4% drift
- **ArPlus**: 40.9% → CH4+, 26.3% → CH3+, 23.4% → ArHPlus + CH3

**Problem**: Ions are being lost rapidly via:
1. **Drift to walls** (mobility × E_field / L_discharge)
2. **Ion-neutral reactions** (creating other ions)
3. **Wall sticking** (minor)

But there's **NO compensating electron-ion recombination** for major ions!

### 3. Missing Dissociative Recombination Reactions

**Ions WITHOUT dissociative recombination in model**:
1. **CH4+** (313 million cm⁻³, 38.4% of all ions!) ← **CRITICAL MISSING!**
2. CH2+ (62k cm⁻³, trace)
3. CHPlus (negative density, artifact)
4. H3Plus (negative density, artifact)

**Currently tuned recombination**: Only 1 out of 25 reactions!
- Only C2HPlus + e → C2 + H (minor ion, 0.01% of total)

**NOT being tuned** (major ions):
- CH5+ + e → products (26.7% of ions)
- CH3+ + e → products (5.5% of ions)  
- C2H5+ + e → products (11.5% of ions)
- ArPlus + e → Ar (11.0% of ions)

### 4. Current Optimization Strategy Gaps

V5 optimization (106 tunable rates) includes:
- ✓ Wall sticking: 28 rates
- ✓ Volume loss: 25 rates
- ✓ Electron-impact: 18 rates (but NO ionization - fixed!)
- ✗ **Dissociative recombination: 1 rate** (should be 25!)
- ✗ **Ion mobilities: 0 rates** (drift = mobility × E / L, NOT tunable)

## Ion Chemistry Breakdown (V5 Result)

| Ion | Density (cm⁻³) | % of Total | Top Loss Mechanism | Has Recomb? |
|-----|----------------|------------|-------------------|-------------|
| CH4+ | 3.13×10⁸ | 38.4% | Drift: 98%+ | **NO** ❌ |
| CH5+ | 2.18×10⁸ | 26.7% | Drift: 98.6% | YES (not tuned) |
| C2H5+ | 9.42×10⁷ | 11.5% | Drift: unknown | YES (not tuned) |
| ArPlus | 8.99×10⁷ | 11.0% | Ion-neutral: 90.6% | YES (not tuned) |
| ArHPlus | 5.68×10⁷ | 7.0% | Ion-neutral: 100% | NO |
| CH3+ | 4.47×10⁷ | 5.5% | Ion-neutral: 80.6%, Drift: 19.4% | YES (not tuned) |

**Total positive ions**: 8.17×10⁸ cm⁻³  
**Electrons**: 1.43×10⁹ cm⁻³  
**Negative ions**: 6.74×10⁷ cm⁻³ (HMinus: 99.8%)

## Recommendations

### IMMEDIATE ACTIONS

1. **Add CH4+ dissociative recombination reaction**
   - Literature: CH4+ + e → CH3 + H (or other channels)
   - Typical rate: ~10⁻⁷ cm³/s (need to verify)
   - **This is THE MOST CRITICAL missing reaction**

2. **Enable tuning for ALL dissociative recombination reactions**
   - Current: 1/25 tuned
   - Recommended: At least tune for major ions (CH5+, CH3+, C2H5+, ArPlus)
   - These directly control electron consumption!

3. **Investigate drift/mobility physics**
   - Current: drift = mobility × E_field / L_discharge (NOT tunable)
   - Question: Are mobility values correct?
   - Consider: Should mobilities be tunable within reasonable bounds?

4. **Increase recombination weights in objective function**
   - Current charge penalty: 75x (moderate)
   - Consider: 200-500x to force charge balance

### SECONDARY ACTIONS

5. **Add ArHPlus dissociative recombination**
   - Currently missing, 7% of ions
   - ArHPlus + e → Ar + H

6. **Review electron attachment reactions**
   - Currently: e + CH3 → CH3Minus (2 reactions total)
   - HMinus production is strong (67 million cm⁻³)
   - These consume electrons but create negative ions

7. **Consider 3-body recombination**
   - For high-density regions
   - Ion+ + e + M → neutral + M

## Physics Questions

1. **Why is drift so dominant?**
   - E_field = 62.6 V/cm in V5
   - Is this realistic for the discharge?
   - Typical short glow discharge: 10-100 V/cm

2. **What's the actual electron temperature?**
   - Affects ionization rates
   - Affects recombination rates
   - Currently using fixed rates

3. **Is charge neutrality required?**
   - Sheath region: ion/e ~4.5x (OK to have deficit)
   - But 0.52x is VERY far from quasi-neutrality
   - May need to reconsider target

## Expected Impact of Fixes

If CH4+ recombination is added and tuned:
- Electron loss will increase ~40% (from CH4+ recombination)
- Ion/electron balance will improve significantly
- May achieve ion/e ratio closer to 2-3x (still below 4.5x due to sheath)
- Overall objective function should improve

## Literature to Review

1. CH4+ dissociative recombination coefficients
2. Ion mobility values in Ar/CH4 mixtures
3. Sheath physics in short glow discharges
4. Charge balance in similar 0D models
