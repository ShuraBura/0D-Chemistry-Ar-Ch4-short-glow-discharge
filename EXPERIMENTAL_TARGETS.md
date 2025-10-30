# Experimental TALIF Targets - Region Analysis

## Experimental Data Summary

### Cathode Glow (CG) Region - Updated 2025-10-30 (Measured Spatial Averages)
| Species | Density (m⁻³) | Density (cm⁻³) | Notes |
|---------|---------------|----------------|-------|
| H       | 9.58e21       | **9.58e15**    | High - production region (measured avg 0-1mm) |
| CH      | 4.29e14       | **4.29e8**     | Low - transient (measured avg 0-1mm) |
| C₂      | 1.50e17       | **1.50e11**    | Moderate (measured avg 0-1mm) |

**Previous estimates (2025-10-29):** H: 8.58e15, CH: 4.6e8, C₂: 1.44e11 cm⁻³
**Agreement:** All within ±12% - excellent validation!

**Physics:**
- ❌ ne: Low, hard to measure
- ❌ Te: NOT well-defined, non-thermal (high-E tail)
- ✓ H production: Local chemistry dominates
- ⚠️ Diffusion: Matters, but chemistry sets ratios

**E-field**: Higher (near cathode)
**Pressure**: 0.4 Torr
**Suggested Te**: ~1-2 eV (effective for high-E tail)

---

### Cathode Sheath Boundary (CSB) Region - Updated 2025-10-29
| Species | Density (m⁻³) | Density (cm⁻³) | Notes |
|---------|---------------|----------------|-------|
| H       | 6.35e20       | **6.35e14**    | Lower - diffused from CG |
| CH      | 9.27e14       | **9.27e8**     | 2.0× higher than CG |
| C₂      | 5.56e17       | **5.56e11**    | 3.9× higher than CG |

**Physics:**
- ✓ ne: Easier to measure
- ✓ Te: Better defined, more thermal
- ❌ H production: Mostly from diffusion, not local
- ⚠️ Chemistry: C₂/CH chemistry active, but H is imported

**E-field**: Lower (away from cathode)
**Suggested Te**: ~0.5-1 eV (bulk plasma)

---

## Model Configuration Recommendations

### Option 1: Target CG (Cathode Glow) ✓ RECOMMENDED for 0-D

**Rationale**:
- Chemistry dominates over transport
- H is produced locally (not diffused)
- 0-D model is chemistry-focused
- Better suited for steady-state chemical balance

**Model Parameters**:
```python
params = {
    'P': 0.4,           # Torr
    'Tg': 400,          # K (or 300 K?)
    'ne': 5e9,          # cm⁻³ (low, but ionization active)
    'Te': 1.5,          # eV (effective Te for non-thermal EEDF)
    'E_field': 80,      # V/cm (higher near cathode)
    'L_discharge': 0.45,# cm (actual gap)
}
```

**Targets**:
- H:  8.58e15 cm⁻³
- CH: 4.6e8 cm⁻³
- C₂: 1.44e11 cm⁻³

**Key ratios**:
- H/CH = 1.87e7 (H dominates!)
- C₂/CH = 313 (C₂ >> CH)
- H/C₂ = 5.96e4

---

### Option 2: Target SE (Sheath Edge)

**Rationale**:
- Te easier to define
- ne easier to measure
- BUT: H is diffused, not chemically produced here

**Model Parameters**:
```python
params = {
    'P': 0.4,
    'Tg': 400,
    'ne': 2e10,         # cm⁻³ (higher)
    'Te': 0.8,          # eV (more thermal)
    'E_field': 50,      # V/cm (bulk)
    'L_discharge': 0.45,
}
```

**Targets**:
- H:  6.35e14 cm⁻³
- CH: 9.27e8 cm⁻³
- C₂: 5.56e11 cm⁻³

**Key ratios**:
- H/CH = 6.85e5 (still H-dominated, but less so)
- C₂/CH = 600 (C₂ >> CH)
- H/C₂ = 1.14e3

**Problem**: 0-D model will predict H production, but here H is mainly transported in

---

## Diffusion Length Considerations

### Cathode Glow (CG)
- **Diffusion length (measured 2025-10-30)**: L_diff_H = **0.057 cm** (fitted from exponential decay)
- **Previous estimate**: L_diff ~ 0.77 cm (from D/k_loss)
- **Interpretation**: MUCH shorter than estimated! Rapid wall loss dominates
- For CH, C₂: Use default L_diff = 0.1 cm (insufficient data for fitting)

### Sheath Edge (SE)
- Closer to bulk plasma
- H diffusion from CG dominates local production
- **Interpretation**: Need to include H flux as boundary condition

---

## Recommended Approach: CG with Te Considerations

### Why CG?
1. ✓ 0-D model is chemistry-based (not transport)
2. ✓ CG is production region (chemistry sets ratios)
3. ✓ Can use effective Te for non-thermal EEDF
4. ✓ Wall loss coefficients already updated

### How to Handle Non-Thermal Te?

The audit is correct: Te-dependence is critical, BUT for non-thermal EEDF in CG:

**Option A: Effective Te** (current approach)
- Use Te = 1.5-2 eV as "effective" temperature
- Weights high-E tail contributions
- Simple, reasonable for 0-D model

**Option B: Two-Temperature Model**
- Separate T_bulk (0.5 eV) and T_tail (3-5 eV)
- Weight electron-impact rates by tail fraction
- More complex, but more accurate

**Option C: BOLSIG+ with High-E tail**
- Use BOLSIG+ with E/N appropriate for cathode region
- Generates <σv> including tail effects
- Most accurate, requires external tool

### Recommended Starting Point

```python
# Focus on CG region with effective Te
params = {
    'P': 0.4,              # Torr
    'Tg': 400,             # K
    'ne': 5e9,             # cm⁻³ (low, CG)
    'Te': 1.5,             # eV (effective, accounts for tail)
    'E_field': 80,         # V/cm (cathode region)
    'L_discharge': 0.45,   # cm
}

# Targets (CG) - Updated 2025-10-30 (Measured Spatial Averages)
target_H  = 9.58e15  # cm⁻³
target_CH = 4.29e8   # cm⁻³
target_C2 = 1.50e11  # cm⁻³

# Diffusion lengths (measured)
L_diff_H  = 0.057    # cm (fitted from profile)
L_diff_CH = 0.1      # cm (default - insufficient data)
L_diff_C2 = 0.1      # cm (default - insufficient data)
```

---

## Next Steps

1. **Run simulation with CG parameters**
2. **Compare ratios**: H/CH, C₂/CH (more robust than absolute values)
3. **If needed, tune**:
   - Te (1.0 → 2.0 eV range)
   - E-field (50 → 100 V/cm)
   - Wall loss rates (already improved)
4. **Advanced**: Implement two-temperature EEDF if simple Te doesn't work

---

## Key Insight

The **ratio** C₂/CH = 408 (CG) or 752 (SE) tells us:
- C₂ is much more stable than CH
- CH is very reactive (short lifetime)
- Model must balance:
  - CH production (from e + CH4, CH + radicals)
  - CH loss (to C₂H₂, C₂H₄, walls)

This ratio is **chemistry-dominated**, so 0-D model should capture it well!

---

## Question for You

**Which region do you want to model?**
- **CG**: Better for 0-D chemistry, but Te is non-thermal
- **SE**: Easier Te/ne, but H is diffused not produced

I recommend **CG** because your model is chemistry-focused.
