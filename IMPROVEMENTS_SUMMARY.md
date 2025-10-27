# Plasma Model Improvements - Summary

## Date: 2025-10-27
## Based on: Complete Audit Recommendations

---

## Changes Implemented

### 1. Wall Loss Rates (FIXED) ✓
**Issue**: Wall loss rates were 5-10× too high for P=0.4 Torr, L=0.45 cm

**Changes**:
- Reduced all wall loss rates from 1e3-1e4 s⁻¹ to **100-800 s⁻¹**
- Now consistent with geometry: k_loss = γ × v_th / L_eff

**Examples**:
| Species | Old Rate (s⁻¹) | New Rate (s⁻¹) | Change |
|---------|----------------|----------------|--------|
| CH      | 1.0e4         | 8.0e2          | -92%   |
| H       | (via stick)   | (via stick)    | -      |
| CH3     | 1.2e3         | 3.0e2          | -75%   |
| C2      | 8.0e2         | 2.0e2          | -75%   |

---

### 2. Electron-Impact Rates Corrected (FIXED) ✓
**Source**: Janev & Reiter (2002)

**Rate Corrections**:
| Reaction | Old Rate | New Rate | Source |
|----------|----------|----------|--------|
| e + CH4 → CH3 + H + e | 6.0e-11 | **4.2e-11** | Janev-Reiter |
| e + CH4 → CH2 + H2 + e | 3.0e-11 | **1.1e-11** | Janev-Reiter |
| e + CH4 → CH + H2 + H + e | 3.0e-11 | **0.7e-11** | Janev-Reiter |
| e + C2H6 → C2H4 + H2 + e | 1.5e-11 | **7.0e-12** | Janev-Reiter |
| e + C2H4 → C2H2 + H2 + e | 1.8e-11 | **1.0e-11** | Updated |

---

### 3. Neutral-Neutral Rates Corrected (FIXED) ✓
**Source**: Baulch (2005)

**Rate Corrections**:
| Reaction | Old Rate | New Rate | Source |
|----------|----------|----------|--------|
| CH + CH4 → C2H4 + H | 1.2e-11 | **1.5e-10** | Baulch 2005 |
| CH + CH3 → C2H4 | 8.0e-11 | **1.5e-10** | Baulch 2005 |
| CH + CH → C2H2 + H | 1.0e-10 | **1.0e-10** | Confirmed |

---

### 4. Missing Electron-Impact Reactions (ADDED) ✓
**Source**: Janev-Reiter / Kushner

**New Reactions**:
```python
e + C2H2 → C2H + H + e          # k = 1.0e-11 cm³/s  (Critical for C2 production!)
e + C2H4 → C2H3 + H + e         # k = 8.0e-12 cm³/s
e + C2H6 → C2H5 + H + e         # k = 1.2e-11 cm³/s
```

---

### 5. Missing Neutral-Neutral Reactions (ADDED) ✓
**Source**: Baulch, Kushner

**New Reactions**:
```python
C + C + M → C2 + M              # k = 1.0e-32 cm⁶/s  (Three-body C2 formation)
H + C2H4 → C2H3 + H2            # k = 1.0e-11 cm³/s  (C2H3 production pathway)
```

---

### 6. Missing Three-Body Recombination (ADDED) ✓
**Source**: Standard plasma chemistry databases

**New Reactions**:
```python
H + H + M → H2 + M              # k = 1.0e-32 cm⁶/s  (Updated from 8e-33)
CH3 + H + M → CH4 + M           # k = 5.0e-31 cm⁶/s  (NEW - CH4 reformation)
```

---

## Impact Summary

### Reaction Network Statistics
- **Original reactions**: ~237
- **Updated reactions**: **244**
- **New reactions added**: 7
- **Rate constants corrected**: 12+

### Expected Model Improvements

1. **More accurate CH4 dissociation**
   - Reduced e-impact rates → Less CH4 depletion
   - Better agreement with experiments

2. **Enhanced C2 production pathway**
   - New e + C2H2 → C2H + H channel
   - C + C + M → C2 mechanism
   - Should help reach C2 target: 1.12e17 m⁻³

3. **Better H balance**
   - Reduced wall losses
   - CH3 + H recombination
   - Should help reach H target: 8.57e21 m⁻³

4. **Improved CH production**
   - Corrected CH + CH4 rate (12× increase!)
   - Should help reach CH target: 2.76e14 m⁻³

---

## Files Modified

1. **define_rates.py**
   - Updated 12 rate constants
   - Added 7 new rate definitions
   - Reduced wall loss rates (Group 11)
   - Updated comments with sources

2. **build_reactions.py**
   - Added 7 new reaction stoichiometries
   - All new reactions linked to rate constants

---

## Validation Status

✓ Modules import successfully
✓ Reaction network builds (244 reactions)
✓ Rate constants validated:
  - e + CH4 → CH3: 4.2e-11 cm³/s ✓
  - Wall loss CH: 800 s⁻¹ ✓
  - CH + CH4 → C2H4: 1.5e-10 cm³/s ✓

---

## Still TODO (From Audit)

### High Priority
1. **Implement Te-dependent rates** (requires BOLSIG+ or lookup tables)
   - Currently all rates assume Te = 1 eV
   - Should use <σv>(Te) from cross-section data
   - Impact: Critical for accurate electron chemistry

2. **Add missing species** (if needed)
   - CH2(singlet) - for C2H4 formation chemistry
   - Ar*(2p) - might already be represented by ArStar
   - H2(v=1-4) - vibrational states (optional)

### Medium Priority
3. **Validate against TALIF targets**
   - Run simulation to steady-state
   - Compare: H, CH, C2 densities
   - Tune parameters if needed

4. **Benchmark performance**
   - Compare simulation speed before/after changes
   - Check for numerical stability issues

---

## How to Use

```bash
# Run the updated simulation
python3 main.py

# The model now includes:
# - More accurate rate coefficients
# - Expanded reaction network
# - Physically correct wall losses
```

---

## References

- Janev & Reiter (2002-2004): CHy/C2Hy cross-sections
- Baulch et al. (2005): Neutral hydrocarbon chemistry
- Kushner models: J. Appl. Phys. 73, 2003; 88, 2226
- Phelps (1999): Ar cross-sections
- Anicich (2003): Ion-molecule reactions
- NIST Kinetics Database

---

## Next Steps

1. **Run full simulation** to 100s steady-state
2. **Compare with TALIF targets**:
   - H: 8.57e21 m⁻³ (target)
   - CH: 2.76e14 m⁻³ (target)
   - C2: 1.12e17 m⁻³ (target)
3. **Implement Te-dependence** if large discrepancies remain
4. **Consider parameter tuning**: Te, E-field, pressure sensitivity

---

**Status**: All immediate improvements from audit completed ✓
**Next milestone**: Validate against experimental data
