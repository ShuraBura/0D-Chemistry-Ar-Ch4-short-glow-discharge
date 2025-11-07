# Temperature Scaling Update for define_rates.py

## Date: 2025-11-07

## Summary
Updated `define_rates.py` to include proper temperature dependence for electron-impact, ionization, and recombination reactions. This addresses the critical issue identified in the temperature dependence audit.

## Changes Made

### 1. Added Temperature Scaling Functions

Three new internal functions added to `define_rates()`:

```python
def scale_electron_impact(k_ref, Te, Te_ref=1.0, E_threshold=None):
    """Scale electron-impact rates: k ~ sqrt(Te) * exp(-E_threshold/Te)"""

def scale_ionization(k_ref, Te, Te_ref=1.0, E_ion=12.0):
    """Scale ionization rates: k ~ sqrt(Te) * exp(-E_ion/Te)"""

def scale_recombination(k_ref, Te, Te_ref=1.0, alpha=0.7):
    """Scale recombination rates: k ~ Te^(-alpha)"""
```

### 2. Temperature-Dependent Rate Groups

**Group 1: Electron-Impact Dissociation (24 rates)**
- Now scale with Te using threshold energies
- Example: `k['e_CH4_CH3_H_cm3_1_1'] = scale_electron_impact(4.2e-11, Te, E_threshold=8.5)`

**Group 2: Ionization (9 rates)**
- Now scale with Te using ionization thresholds
- Example: `k['e_CH4_CH3Plus_H_cm3_2_1'] = scale_ionization(1e-11, Te, E_ion=12.6)`

**Group 6: Dissociative Recombination (21 rates)**
- Now scale as Te^(-alpha)
- Example: `k['ArPlus_e_Ar_cm3_6_1'] = scale_recombination(1.5e-7, Te, alpha=0.7)`
- Ion-ion recombination remains temperature-independent

### 3. Temperature-Independent Groups (Unchanged)

- **Group 3**: Ar* reactions (thermal energies)
- **Group 4**: Penning ionization (excitation energy driven)
- **Group 5**: Ion-neutral reactions (thermal)
- **Group 7**: Neutral-neutral reactions (thermal)
- **Group 8**: Termolecular recombination (thermal)
- **Groups 9-11**: Wall reactions and losses (system-dependent)

## Backward Compatibility

✅ **Fully backward compatible**

- Default Te = 1.0 eV if not specified
- Default Tgas = 400 K if not specified
- All existing scripts work without modification
- Results identical at Te=1.0 eV

```python
# Old usage (still works, defaults to Te=1.0 eV)
params = {'E_field': 50, 'L_discharge': 0.45, 'mobilities': {...}}
k = define_rates(params)

# New usage with custom Te
params = {'Te': 2.0, 'E_field': 50, 'L_discharge': 0.45, 'mobilities': {...}}
k = define_rates(params)
```

## Verification Results

### Temperature Dependence Tests

**At Te=2.0 eV vs Te=1.0 eV:**

| Rate Type | Example | Change | Expected | Status |
|-----------|---------|--------|----------|--------|
| Electron-impact | e_CH4_CH3_H | 99× increase | Increase | ✓ |
| Ionization | e_Ar_ArPlus | 3739× increase | Increase | ✓ |
| Recombination | ArPlus_e_Ar | 0.62× decrease | Decrease | ✓ |
| Neutral-neutral | CH_CH3_C2H4 | 0× change | Constant | ✓ |

### Literature Range Verification

**At Te=1.0 eV (reference condition):**
- Total rates: 256
- Rates with literature data: 227
- Within literature range: 180 (79%)
- Out of range: 47 (21%)

**Out-of-range rates are primarily:**
- Loss rates (system-specific, 19 rates)
- Stick reactions (system-specific, 2 rates)
- Ion-ion recombination (slightly above max, 12 rates)
- Some ion-neutral rates (slightly below min, 3 rates)

These are acceptable as they are system-dependent or have been tuned to experimental data.

## Impact on Previous Results

### Scripts Using define_rates.py (92 files)

**LOW IMPACT** if Te was never varied:
- Results remain valid if Te ≈ 1 eV was used
- Rates were calibrated for Te=1 eV

**HIGH IMPACT** if Te was varied in optimizations:
- Previous optimizations that varied Te may have incorrect electron-impact and recombination rates
- Need to re-verify optimal Te values
- True temperature effects were not captured

### Scripts Using define_rates_tunable.py (62 files)

**NO IMPACT**:
- Already had correct temperature scaling
- Results remain valid

## Physical Correctness

### Before Update
❌ Electron-impact rates constant regardless of Te
❌ Ionization rates constant regardless of Te
❌ Recombination rates constant regardless of Te
❌ Cannot study Te effects properly

### After Update
✓ Electron-impact rates increase with Te (exponential + power law)
✓ Ionization rates increase strongly with Te (threshold behavior)
✓ Recombination rates decrease with Te (power law)
✓ Can properly study Te effects

## Testing

### Test Script Provided
Run `verify_rates_against_literature.py` to:
- Check all rates against literature ranges
- Verify temperature scaling behavior
- Test at multiple Te values (0.5, 1.0, 2.0, 3.0 eV)

```bash
python3 verify_rates_against_literature.py
```

### Quick Test
```bash
python3 -c "
from define_rates import define_rates

params = {
    'E_field': 50,
    'L_discharge': 0.45,
    'mobilities': {...}
}

# Test default (Te=1.0)
k1 = define_rates(params)

# Test high Te
params['Te'] = 2.0
k2 = define_rates(params)

print(f'Electron-impact increase: {k2[\"e_CH4_CH3_H_cm3_1_1\"]/k1[\"e_CH4_CH3_H_cm3_1_1\"]:.1f}x')
print(f'Recombination decrease: {k2[\"ArPlus_e_Ar_cm3_6_1\"]/k1[\"ArPlus_e_Ar_cm3_6_1\"]:.2f}x')
"
```

## Recommendations

### For New Work
✅ Use updated `define_rates.py` with explicit Te parameter

### For Existing Results
1. **If Te was never varied**: Results are valid
2. **If Te was varied**: Consider re-running critical optimizations
3. **For publications**: Cite this update and verify key conclusions

### For Future Development
Consider adding Tgas-dependent neutral-neutral rates:
```python
# Future enhancement
k_neutral = k_ref * (Tgas/T_ref)**n  # Arrhenius-type
```

## Files Modified
1. `define_rates.py` - Added temperature scaling
2. `verify_rates_against_literature.py` - New verification script
3. `TEMPERATURE_DEPENDENCE_AUDIT.md` - Original audit report
4. `TEMPERATURE_SCALING_UPDATE.md` - This document

## References

- Morgan et al. (1992) - Electron-impact cross sections
- Janev & Reiter (2002) - Atomic and molecular processes
- Phelps (1999) - Argon cross sections
- Baulch et al. (2005) - Combustion chemistry
- UMIST (2012) - Recombination rate database

---
**Author**: Claude Code
**Date**: 2025-11-07
**Commit**: Temperature scaling update to define_rates.py
