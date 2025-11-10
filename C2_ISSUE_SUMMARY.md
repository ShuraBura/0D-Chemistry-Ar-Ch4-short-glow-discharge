# C2 Production Issue - Root Cause Analysis

## Problem Summary
C2 density collapses to near-zero at steady state (2.65e7 cm⁻³ vs target 5.60e11, only 0.005% of target).

## Investigation Results

### 1. CRITICAL FIX Status ✓
The "CRITICAL FIX: Force-include C2-producing reactions" **IS WORKING**:
- All 6 C2-producing reactions are correctly force-included in optimizer
- Reactions are in the tunable set and will be optimized

### 2. ROOT CAUSE IDENTIFIED ⚠️
**C2 wall loss dominates production by 32×:**

```
C2 Production:  7.10e+10 cm⁻³/s
C2 Wall Loss:   2.28e+12 cm⁻³/s
Ratio:          0.031 (only 3.1%)
```

**C2 Production Sources (with default rates):**
1. e + C2H2 → C2 + H2:     6.87e+10 cm⁻³/s (97% of production!)
2. CH + CH → C2 + H2:      8.30e+08 cm⁻³/s
3. C + CH → C2 + H:        1.51e+09 cm⁻³/s
4. C2H+ + e → C2 + H:      2.38e+07 cm⁻³/s

**C2 Wall Loss:**
- loss_C2_11_3: k = 200 s⁻¹ (default)
- In best result: k = 671 s⁻¹ (optimizer made it WORSE!)

## Why C2 Collapses at Steady State

At transient state (from best_f49.5.json):
- C2 = 1.14e+10 cm⁻³ (still has some C2)

At true steady state:
- C2 → 2.65e+07 cm⁻³ (collapsed by 430×!)
- C2 production → 0 because C2H2 and other precursors are consumed

## Solution Options

### Option A: Reduce C2 Wall Loss (Recommended)
1. **Lower loss_C2_11_3** bound from 200-5000 to 10-100 s⁻¹
2. **Lower stick_C2_9_9** (C2 sticking coefficient)
3. This allows C2 to accumulate to target levels

### Option B: Boost C2 Production
1. **Increase e + C2H2 → C2 + H2** rate (currently contributes 97%)
2. **Increase CH + CH → C2 + H2** rate
3. Need to boost by 32× to match current wall loss

### Option C: Fix Both (Best)
1. Reduce C2 wall loss by 10×
2. Boost C2 production by 3×
3. This gives margin for steady state

## Recommended Next Steps

1. **Modify rate_database_complete.py:**
   ```python
   # Reduce C2 wall loss
   'loss_C2_11_3': RateInfo(min=10.0, max=100.0, ...)  # was 200-5000
   'stick_C2_9_9': RateInfo(min=100.0, max=1000.0, ...) # was 3049
   ```

2. **Run optimizer with fix:**
   ```bash
   python3 optimize_charge_balanced_500mTorr_NO_TIMEOUT.py
   ```

3. **Monitor C2 production/loss balance** in results

## Key Insight
The CRITICAL FIX correctly force-includes C2 reactions, but the fundamental issue is **wall loss overwhelms production**. Fixing the rate bounds will allow the optimizer to find physically reasonable C2 densities.
