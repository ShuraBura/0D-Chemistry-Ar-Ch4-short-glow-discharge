# Optimization Status Report
## Session: 2025-11-10

### Executive Summary

Successfully restarted work after code execution failure. Root cause of C2 collapse identified and fixes implemented. Multiple optimization runs in progress.

---

## Critical Discoveries

### Discovery 1: Unbalanced Objective Function Weights
**Problem:** Original optimizer had:
- H weight: 20.0
- CH weight: 20.0
- **C2 weight: 3.0** ← Only 15% of H/CH importance!

**Impact:** Optimizer sacrificed C2 to improve H and CH, leading to maximum C2 wall loss rates

**Fix:** Created `optimize_charge_balanced_500mTorr_FIXED_WEIGHTS.py` with C2 weight = 20.0

---

### Discovery 2: Transient vs Steady State Mismatch
**Problem:** Optimizer integrated to t=100s, which is NOT steady state!

**Evidence from best_f109.6.json:**

| Time | H (cm⁻³) | CH (cm⁻³) | C2 (cm⁻³) | Status |
|------|----------|-----------|-----------|--------|
| t=100s | 1.06e+13 (4%) | 2.90e+09 (290%) | 1.43e+10 (3%) | Transient |
| t=1000s | 3.66e+14 (145%) | 7.33e+09 (733%) | 1.19e+08 (0.02%) | **TRUE Steady State** |

**Impact:**
- Optimizer sees H at 4% of target → tries to increase H
- True steady state has H at 145% of target → overshoots!
- **34× discrepancy in H density** between t=100s and t=1000s
- C2 collapses by 120× from t=100s to t=1000s

**Fix:** Created `optimize_charge_balanced_500mTorr_STEADY_STATE.py` with integration to t=500s

---

### Discovery 3: Force-Included C2 Reactions Are Working
**Status:** ✓ VERIFIED

All 6 C2-producing reactions ARE in the tunable set:
- `CH_CH_C2_H2_cm3_5_4`: CH + CH → C2 + H2
- `e_C2H2_C2_H2_cm3_1_16`: e + C2H2 → C2 + H2 (electron-impact!)
- `C2HPlus_e_C2_H_cm3_6_18`: C2H+ + e → C2 + H (ion recombination)
- `C_CH_C2_H_cm3_7_4`: C + CH → C2 + H
- `CH_C_C2_H_cm3_7_9`: CH + C → C2 + H
- `C2H_H_C2_H2_cm3_7_47`: C2H + H → C2 + H2

The CRITICAL FIX from commit 19a53fe is working correctly.

---

### Discovery 4: C2 Wall Loss Dominates Production

**At best result (true steady state):**
- C2 production: ~3.6e+11 cm⁻³/s (mainly e + C2H2 → C2)
- C2 wall loss: ~5.4e+11 cm⁻³/s
- **Production/Loss ratio: 0.67**

**Optimizer choices:**
- C2 volumetric loss: 1.06e+03 s⁻¹ (allowed: 0.0001-2000)
- C2 wall sticking: 3.46e+03 s⁻¹ (allowed: 1250-6250)
- Both near maximum → optimizer prefers high C2 loss

---

## Current Status

### Three Optimizers Created:

1. **optimize_charge_balanced_500mTorr_NO_TIMEOUT.py**
   - Integration: t=100s
   - C2 weight: 20.0 ✓
   - Force-includes C2 reactions ✓
   - **Issue:** Optimizes for transient state

2. **optimize_charge_balanced_500mTorr_FIXED_WEIGHTS.py**
   - Integration: t=100s
   - C2 weight: 20.0 ✓
   - Force-includes C2 reactions ✓
   - **Status:** Stopped after discovering transient issue
   - **Results:** 3 results saved, best f(x) = 154.2

3. **optimize_charge_balanced_500mTorr_STEADY_STATE.py** ← CURRENTLY RUNNING
   - Integration: **t=500s** ✓
   - C2 weight: 20.0 ✓
   - Force-includes C2 reactions ✓
   - **Status:** Running in background
   - **Current best:** f(x) = 194.9 after ~10 evaluations
   - **Performance:** ~8-10 evaluations/min (5× slower than t=100s version)

---

## Analysis of Current Results

### Steady-State Optimizer (t=500s)
**Best result: f(x) = 194.9**
- H: 1.07e+13 (4.2% of target)
- CH: 3.72e+09 (372% of target)
- C2: 2.19e+10 (3.9% of target)
- Ni/Ne: 1.55

**Observations:**
- H still at ~1e+13 even at t=500s for current parameter sets
- May need t>500s for some parameters to reach true steady state
- Or H evolution dynamics vary significantly with parameters

### Key Chemistry Insights

**Precursor densities in best results:**
- C2H2: ~1.6e+12 cm⁻³ (excellent for e + C2H2 → C2)
- CH: ~3-4e+09 cm⁻³ (good for CH + CH → C2)
- e: ~2.1e+09 cm⁻³ (good electron density)
- C: ~2e+10 cm⁻³ (good for C + CH → C2)

**Reaction rates in best result:**
- CH + CH → C2: k = 1.62e-10 cm³/s
- e + C2H2 → C2: k = 6.83e-11 cm³/s (dominant pathway)
- C2H+ + e → C2: k = 2.19e-07 cm³/s (fast but low C2H+ density)

---

## Remaining Issues

### Issue 1: CH Overshoot
At true steady state, CH reaches 7.33e+09 (733% of target).

**Impact on objective function:**
- CH error dominates: (6.33)² = 40.1
- H error: (0.45)² = 0.20
- C2 error: (-1.0)² = 1.00
- Total species error: 20×(40.1 + 0.20 + 1.00) = 826

**Why this matters:** Optimizer will try to reduce CH, which may drive carbon into other species and away from C2.

### Issue 2: Integration Time Uncertainty
- t=100s: definitely transient (H increases 34× from t=100s to t=1000s)
- t=500s: may still be transient for some parameter sets
- t=1000s: appears to be steady state

**Tradeoff:**
- Longer integration → more accurate steady state
- Longer integration → slower optimization (t=1000s would be 10× slower)

### Issue 3: C2 Wall Loss
Even with equal weights and force-included reactions, optimizer chooses high C2 wall loss rates.

**Possible solutions:**
1. Reduce C2 wall loss upper bounds (stick_C2_9_9: 6250 → 2000)
2. Add penalty term for high C2 loss rates
3. Constrain C2 wall loss to measured values if available

---

## Recommendations

### Option A: Let Current Optimizer Run Longer
- Steady-state optimizer (t=500s) is running
- Give it several hours to explore parameter space
- Monitor for convergence

**Pros:**
- Already running, no additional setup
- May find good solution with current constraints

**Cons:**
- Slow (5× longer than t=100s)
- t=500s may still not be true steady state for all parameters

### Option B: Increase Integration Time to t=1000s
- Guarantees true steady state
- Eliminates transient vs steady state issue

**Pros:**
- Optimizes for correct steady state
- Eliminates uncertainty

**Cons:**
- 10× slower than t=100s version
- May take days to converge

### Option C: Constrain C2 Wall Loss + Use t=500s
- Modify steady-state optimizer to reduce C2 loss bounds
- Keep t=500s integration

**Changes:**
```python
# In rate_database_complete.py
'stick_C2_9_9': RateInfo(min=500.0, max=2000.0, ...)  # was 1250-6250
```

**Pros:**
- Prevents optimizer from maximizing C2 loss
- Still uses t=500s (reasonable speed)
- Forces C2 to accumulate

**Cons:**
- May be overconstraining if wall loss is truly high
- Need experimental data to validate bounds

### Option D: Use Adaptive Integration Time
- Start with t=100s for initial exploration
- Switch to t=500s or t=1000s for final refinement
- Check convergence: if |dy/dt| / |y| < 1e-10, consider converged

**Pros:**
- Fast initial exploration
- Accurate final optimization
- Detects when steady state is reached

**Cons:**
- More complex implementation
- Requires modifying optimizer

---

## Diagnostic Scripts Created

1. **test_C2_production_fix.py**
   - Verifies force-included reactions are in tunable set ✓
   - Analyzes C2 production vs loss balance
   - Calculates production/loss ratio

2. **check_steady_state_best.py**
   - Extends integration to t=1000s for any result file
   - Compares optimizer-reported vs true steady state
   - Verifies convergence (dN/dt ≈ 0)

3. **monitor_optimization.py**
   - Tracks optimization progress
   - Shows top 5 results
   - Analyzes C2-producing reaction rates and precursor densities
   - Checks for convergence

---

## Files Created This Session

- `C2_ISSUE_SUMMARY.md`: Root cause analysis
- `test_C2_production_fix.py`: Verification and diagnostics
- `optimize_charge_balanced_500mTorr_FIXED_WEIGHTS.py`: Equal weights fix
- `optimize_charge_balanced_500mTorr_STEADY_STATE.py`: Steady-state integration
- `check_steady_state_best.py`: Steady-state verification
- `monitor_optimization.py`: Progress monitoring
- `OPTIMIZATION_STATUS_REPORT.md`: This file

---

## Next Steps

**Immediate (< 1 hour):**
1. Monitor steady-state optimizer progress
2. Check if results improve with t=500s integration
3. Verify C2 production mechanisms are working

**Short-term (1-4 hours):**
1. Let steady-state optimizer run for several hours
2. Analyze top results for true steady state behavior
3. Determine if t=500s is sufficient or need t=1000s

**If current optimizer doesn't solve C2 issue:**
1. Implement Option C (constrain C2 wall loss bounds)
2. Or implement Option D (adaptive integration time)
3. May need experimental data to validate C2 loss rates

---

## Key Metrics to Monitor

1. **Objective function value** - should decrease over time
2. **Species errors at TRUE steady state** (not t=100s values!)
   - H: target 2.52e+14
   - CH: target 1.0e+9
   - C2: target 5.6e+11
3. **C2 production/loss ratio** - should be ≥ 1.0 for steady state
4. **Convergence spread** - top 10 results should cluster

---

## Success Criteria

**Minimum acceptable:**
- H: 2.0-3.0e+14 (80-120% of target)
- CH: 0.5-2.0e+9 (50-200% of target)
- C2: 4.0-7.0e+11 (70-125% of target)
- Ni/Ne: 2-7
- All species at true steady state (dN/dt < 1e-10 × N)

**Ideal:**
- H: 2.4-2.6e+14 (95-105% of target)
- CH: 0.8-1.2e+9 (80-120% of target)
- C2: 5.0-6.2e+11 (90-110% of target)
- Ni/Ne: 3-5
- Objective function < 50

---

*Report generated: 2025-11-10*
*Steady-state optimizer status: RUNNING*
*Estimated time to completion: 4-8 hours*
