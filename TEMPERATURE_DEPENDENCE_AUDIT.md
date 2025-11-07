# Temperature Dependence Audit Report

## Summary
**CRITICAL FINDING**: Some simulation scripts use hardcoded rate constants without temperature dependence, while others correctly use temperature-dependent rates. This inconsistency may affect results.

## Rate Definition Files

### 1. `define_rates.py` (❌ NO TEMPERATURE DEPENDENCE)
- **ALL rates are hardcoded constants**
- Does NOT use `Te` (electron temperature) parameter
- Does NOT use `Tgas` (gas temperature) parameter
- Example: `k['e_CH4_CH3_H_cm3_1_1'] = 4.2e-11` (constant)

**Used by:**
- `run_optimized.py`
- `optimize_v7b_checkpointed.py`
- `test_quick.py`
- `main.py`
- `benchmark.py`
- And many other files (92 total)

### 2. `define_rates_tunable.py` (✅ HAS TEMPERATURE DEPENDENCE)
- **Properly implements Te and Tgas dependence**
- Uses scaling functions for different reaction types

**Temperature-dependent reactions:**

1. **Electron-impact dissociation** (Group 1):
   - Uses: `scale_electron_impact(k_ref, Te, E_threshold)`
   - Scaling: `k ~ sqrt(Te) * exp(-E_threshold*(1/Te - 1/Te_ref))`
   - Example: `k['e_CH4_CH3_H_cm3_1_1'] = scale_electron_impact(4.2e-11, Te, E_threshold=8.5)`

2. **Electron-impact ionization** (Group 2):
   - Uses: `scale_ionization(k_ref, Te, E_ion)`
   - Scaling: `k ~ sqrt(Te) * exp(-E_ion*(1/Te - 1/Te_ref))`
   - Example: `k['e_CH4_CH3Plus_H_cm3_2_1'] = scale_ionization(1e-11, Te, E_ion=12.6)`

3. **Dissociative recombination** (Group 6):
   - Uses: `scale_recombination(k_ref, Te, alpha)`
   - Scaling: `k ~ Te^(-alpha)` where alpha ~ 0.5-1.0
   - Example: `k['ArPlus_e_Ar_cm3_6_1'] = scale_recombination(1.5e-7, Te, alpha=0.7)`

4. **Diffusion/loss rates** (Group 11):
   - Depend on Tgas through diffusion coefficient scaling
   - Scaling: `D ~ (Tgas/T_ref)^1.75 / P`
   - Affects all wall loss rates: `k_loss = D / L_diff²`

**Temperature-independent reactions** (as expected):
- Ar* reactions (Group 3) - thermal energies
- Penning ionization (Group 4) - excitation energy driven
- Ion-neutral reactions (Group 5) - thermal
- Neutral-neutral reactions (Group 7) - thermal (could have weak Tg dependence)
- Termolecular recombination (Group 8) - thermal

**Used by:**
- `comprehensive_ch_analysis.py` ✓
- `test_corrected_ch_ch4_rate.py` ✓
- `test_te_e_sweep_optimized.py` ✓
- Many other test/analysis scripts (62 total)

## Checkpoint Files
Checkpoint files (e.g., `checkpoint_f3407.json`) **DO** contain temperature parameters:
```json
{
  "params": {
    "Te": 1.4700762646387704,
    "Tgas": 300,
    ...
  }
}
```

## Problem Analysis

### Issue 1: Inconsistent Rate Definitions
**Severity: HIGH**

Some production scripts use `define_rates.py` which ignores the Te parameter even when it's present in params:
- `run_optimized.py` - Main simulation script
- `optimize_v7b_checkpointed.py` - Optimization script

This means:
- Electron-impact rates don't increase with higher Te
- Ionization rates don't increase with higher Te
- Recombination rates don't decrease with higher Te
- Loss rates don't scale with Tgas

### Issue 2: Test Script Bug
The temperature dependence test script (`test_te_dependence.py`) has a bug:
```python
KeyError: 'H2Plus'
```
Missing H2Plus mobility in test parameters.

## Recommendations

### Option 1: Fix define_rates.py (RECOMMENDED)
**Action:** Update `define_rates.py` to include Te and Tgas dependence like `define_rates_tunable.py`

**Pros:**
- Fixes all scripts using `define_rates`
- Maintains backward compatibility (can default Te=1.0, Tgas=400)
- Physically correct

**Cons:**
- Changes existing behavior
- Need to verify all scripts pass correct parameters

### Option 2: Switch all scripts to define_rates_tunable
**Action:** Change imports in production scripts from `define_rates` to `define_rates_tunable`

**Pros:**
- Uses existing tested code
- Clean separation of concerns

**Cons:**
- Need to update many files
- May break scripts that don't pass Te/Tgas parameters

### Option 3: Document and Accept Current Behavior
**Action:** Document that `define_rates.py` uses rates calibrated for Te=1 eV, Tgas=400 K

**Pros:**
- No code changes
- Existing results remain valid

**Cons:**
- Physically incorrect for varying temperatures
- Can't study Te or Tgas effects with these scripts

## Impact Assessment

### Critical Impact (High Priority)
If optimizations were run with varying Te but using `define_rates.py`:
- Results are **inconsistent** with physics
- Te variations had NO effect on electron-impact rates
- True optimal Te may be different

### Medium Impact
Analysis scripts using `define_rates_tunable.py`:
- Results are **physically correct**
- Temperature effects properly captured

### Low Impact
If all simulations used Te ≈ 1 eV:
- Results are approximately correct (rates calibrated for Te=1 eV)
- But still missing temperature scaling

## Verification Steps

1. Check if any optimization varies Te:
```bash
grep -l "Te.*=" optimize_*.py | head -5
```

2. Check recent best results for Te values:
```bash
ls optimization_results*/*.json | xargs -I {} sh -c 'echo {}; cat {} | grep "\"Te\"" | head -1'
```

3. Verify which rate definition is used in critical scripts:
```bash
grep "from define_rates" run_*.py optimize_*.py
```

## Recommended Fix

Add to the top of `define_rates.py`:
```python
# Extract temperature parameters (default to calibration values)
Te = params.get('Te', 1.0)  # eV
Tgas = params.get('Tgas', params.get('Tg', 400))  # K

# Add scaling functions from define_rates_tunable.py
def scale_electron_impact(k_ref, Te, Te_ref=1.0, E_threshold=None):
    if E_threshold is not None and E_threshold > 0:
        return k_ref * np.sqrt(Te/Te_ref) * np.exp(-E_threshold * (1/Te - 1/Te_ref))
    else:
        return k_ref * (Te/Te_ref)**0.7

# ... etc
```

Then update all electron-impact, ionization, and recombination rates to use these functions.

---

**Audit Date:** 2025-11-07  
**Auditor:** Claude Code  
**Status:** Requires user decision on mitigation strategy
