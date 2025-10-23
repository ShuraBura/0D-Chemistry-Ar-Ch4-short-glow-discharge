# Plasma Chemistry Optimization - Status Report

**Date:** 2025-10-23
**Status:** Optimization Running (Started 14:56, ~2 min elapsed)

---

## Problem Statement

Match experimental target densities for H, CH, C2 at sheath edge:
- **H:** 5.18×10¹³ cm⁻³
- **CH:** 1.0×10⁹ cm⁻³
- **C2:** 1.3×10¹¹ cm⁻³

All rate constants must remain within literature-validated ranges.

---

## Progress Summary

### Phase 1: Analysis ✓ COMPLETE

**Key Finding:** Original simulation had fundamental issues:
1. Electron density 3x too high (1.0×10¹⁰ vs experimental 3.3×10⁹)
2. 29 reaction rates outside literature bounds
3. CH+CH→C2+H2 rate 20% above literature maximum

**Actions Taken:**
- Built complete rate database with 225 rates, literature ranges, citations
- Analyzed all CH production/loss pathways
- Identified 57 rates already at min/max boundaries
- Created comprehensive CH chemistry report

### Phase 2: Corrections ✓ COMPLETE

**Fixed Issues:**
1. Set Ne = 3.3×10⁹ cm⁻³ (experimental value)
2. Corrected all 29 out-of-bounds rates to literature limits
3. Established baseline simulation

**Baseline Results (with corrections):**
```
Species | Target     | Baseline   | Ratio  | Status
--------|------------|------------|--------|--------
H       | 5.18×10¹³  | 3.45×10¹³  | 0.67x  | Need +49%
CH      | 1.0×10⁹    | 4.59×10¹⁰  | 46x    | TOO HIGH!
C2      | 1.3×10¹¹   | 7.31×10¹¹  | 5.6x   | Too high
```

**Interpretation:**
- Corrections helped (CH: 53x → 46x)
- But CH still severely elevated
- **Conclusion:** Need multi-parameter optimization

### Phase 3: Optimization ⏳ IN PROGRESS

**Started:** 14:56 (2 minutes ago)
**Expected duration:** 10-30 minutes
**CPU usage:** 95.8% (actively computing)

**Optimization Strategy:**
- **Parameters:** 31 total
  - 30 reaction rates (selected by impact & range)
  - 1 E field [10, 200] V/cm
- **All rates constrained within literature bounds**
- **Method:** Differential evolution
  - Maxiter: 100
  - Population: 10
  - Strategy: best1bin
- **Objective:** Weighted least squares
  - CH weighted 20x (hardest, 46x off)
  - C2 weighted 3x (5.6x off)
  - H weighted 1x (0.67x, closest)

**Key Tunable Rates:**
```
1. loss_C2_11_3        : C2 wall loss    [1e-4, 2e3]    20 million x range!
2. loss_C2H2Star_11_25 : C2H2* loss      [100, 1e4]     100x range
3. loss_CH_11_9        : CH wall loss    [1e3, 1e4]     10x range ⚠ CRITICAL
4. e_CH4_CH_cm3_1_3    : CH production   [2e-11, 1e-10] 5x range
5. stick_CH_9_3        : CH sticking     [1.25e3, 6.25e3] 5x range
... and 25 more
```

---

## Expected Outcomes

### Best Case
- All three species within ±20% of targets
- Demonstrates targets are achievable within literature constraints
- Provides optimized rate set for future simulations

### Likely Case
- CH improves significantly (46x → 5-10x)
- H and C2 within factor of 2
- Some trade-offs required between species
- Identifies which constraints are limiting

### Worst Case
- Targets mutually incompatible within literature bounds
- Would indicate:
  - Chemistry model incomplete
  - Experimental measurement uncertainty
  - Need for additional reactions or species

---

## Files Generated

### Analysis Tools
- `rate_database_complete.py` - All 225 rates with ranges & citations
- `check_rate_bounds.py` - Validates rates vs literature
- `analyze_ch_chemistry.py` - CH pathway analysis
- `CH_ANALYSIS_REPORT.md` - Comprehensive strategy report

### Optimization Tools
- `correct_rates_to_literature.py` - Brings rates within bounds
- `run_baseline_corrected.py` - Corrected baseline simulation
- `optimize_constrained.py` - Multi-parameter optimization (RUNNING)

### Results (pending)
- `optimized_parameters.txt` - Will contain final optimized rates & E field
- Optimization output - Final densities and convergence info

---

## Next Steps

1. **Wait for optimization completion** (~8-28 minutes remaining)
2. **Analyze results:**
   - Check final H, CH, C2 densities
   - Verify all rates within literature bounds
   - Examine which rates changed most
3. **If successful:**
   - Document optimized parameter set
   - Run validation simulation
   - Compare to experimental data
4. **If unsuccessful:**
   - Analyze why targets incompatible
   - Consider relaxing constraints or examining chemistry model
   - May need sensitivity analysis to identify missing reactions

---

## Literature Sources Referenced

19 unique sources including:
- Baulch et al. (2005) - Hydrocarbon combustion kinetics
- Janev & Reiter (2002) - Electron impact dissociation
- UMIST (2012) - Ion-electron recombination
- Phelps (1999) - Ar* reactions
- Morgan (1992) - Electron-CH4 reactions
- Anicich (2003) - Ion-molecule reactions
- Jauberteau et al. (1998) - Wall sticking coefficients

All rates constrained to peer-reviewed experimental or theoretical values.

---

**Status:** Optimization actively running, results expected within 30 minutes.
