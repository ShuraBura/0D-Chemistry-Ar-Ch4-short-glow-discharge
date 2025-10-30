# Spatial Averaging and H₂ Recombination - Critical Notes

## Spatial Averaging (CRITICAL CLARIFICATION!)

### Your Statement:
**"This is spatially averaged density for all species (in 0-1 mm from cathode)"**

### Implications:

This changes EVERYTHING about model interpretation!

#### What This Means:

TALIF measurements give **spatial average** over 0-1 mm region:
```
n_avg = (1/L) ∫ n(z) dz   from z=0 to z=1 mm
```

NOT a single-point measurement!

#### Why This Helps the 0-D Model:

**0-D model predicts**: Spatially uniform density (by definition)
**Reality**: Density varies with position z

If actual profile is:
- Peak at z=0.5 mm: n_peak = 1.5 × n_avg
- Falls off near cathode (z=0) and edge (z=1 mm)

Then:
- **0-D should match n_avg** (what TALIF measures)
- **NOT n_peak** (which would be higher)

**This is FAVORABLE for model validation!**

---

## H₂ Recombination (Already Included ✓)

### 1. Three-Body Recombination (Gas Phase)

**Already in model** (`define_rates.py`, line 206):
```python
k['H_H_M_H2_M_cm6_8_1'] = 1.0e-32  # cm⁶/s

Reaction: H + H + M → H₂ + M
```

At P = 0.4 Torr, n_total ~ 1e16 cm⁻³:
```
Rate = k × n_H² × n_M
     = 1e-32 × n_H² × 1e16
     = 1e-16 × n_H²  cm⁻³/s
```

For n_H = 1e15 cm⁻³:
```
Rate = 1e-16 × (1e15)² = 1e14 cm⁻³/s
```

This is **slow** compared to wall loss (4.5e17 cm⁻³/s), but not negligible!

### 2. Wall Recombination (Heterogeneous)

**Implemented via γ_H** in `define_rates_tunable.py`:

When H hits copper wall:
- Fraction **γ_H**: Recombines → H₂ (lost from H population)
- Fraction **(1-γ_H)**: Reflects back as H

The H₂ formed at wall returns to gas phase.

**Critical**: H₂ can then be dissociated back to H:
```
e + H₂ → H + H        (k = 6e-12 cm³/s)
ArStar + H₂ → H + H   (k = 6e-11 cm³/s)
```

This creates a **cycle**:
```
H → (wall) → H₂ → (dissociation) → 2H
```

**Net effect**: Effective loss rate is reduced, H density increases!

---

## Copper Surface Contamination

### Initial State: Clean Copper
**γ_H (clean Cu)**: 0.001-0.01 (very low, H reflects)

### During Measurement: Contaminated
**Your note**: "decontaminated during measurements"
(Assuming you meant **contaminated**?)

**Contaminants**:
- Hydrocarbon deposits from CH₄ decomposition
- Oxide layer from residual O₂
- Sputtered material re-deposition

**Effect on γ_H**:
- **Hydrocarbons**: May reduce γ_H further (H adsorbs but doesn't recombine)
- **Oxides**: Typically increase γ_H (more surface sites)
- **Sputtering**: Exposes fresh Cu (cleans surface, reduces γ_H)

**Typical range during operation**:
```
γ_H = 0.001 (fresh, sputtered)
    → 0.01 (steady-state, some contamination)
    → 0.05-0.1 (heavily contaminated, oxides)
```

**For sweep**: Test γ_H = 0.001 to 0.1 to cover range!

---

## E-field Gradient in Cathode Fall

### Why E-Field Should Be Tunable

Cathode fall voltage: 253 V over 4 mm
**Average**: E_avg = 632 V/cm

But E(z) is **highly non-uniform**:
```
E(z) ~ E_0 × exp(-z/λ)   (approximate)
```

**At different positions**:
- z = 0 (cathode surface): E ~ 1000-2000 V/cm (very high!)
- z = 0.5 mm (CG center): E ~ 600-800 V/cm
- z = 1 mm (CG edge): E ~ 400-500 V/cm
- z = 2 mm (bulk): E ~ 100-200 V/cm

**For spatially-averaged model** (0-1 mm):
Should use **effective E** that represents average over region:
```
E_eff ~ 400-800 V/cm
```

**Recommended sweep range**: E = 400, 500, 600, 800 V/cm

---

## H₂ Build-up and Dissociation Cycle

### Critical Feedback Loop

1. **H production**: e + CH₄ → CH₃ + H
2. **H loss to wall**: H → (wall, γ_H) → ½ H₂
3. **H₂ accumulation**: n_H₂ builds up
4. **H₂ dissociation**: e + H₂ → H + H
5. **Back to step 2**: Cycle repeats

**Steady-state H₂**:
```
n_H₂_ss = (H→H₂ rate) / (H₂→H rate)
```

If H₂ builds up to ~1e13 cm⁻³:
```
H₂ dissociation rate = 6e-12 × 1e9 × 1e13 = 6e10 cm⁻³/s

This is H production rate = 1.2e11 cm⁻³/s (doubled!)
```

**This could be critical for closing the gap!**

---

## Updated Sweep Parameters

### Comprehensive Sweep (Realistic)

```python
# Primary variables
ne_values = [3e8, 5e8, 8e8, 1e9]           # cm⁻³ (4 values, MAX 1e9)
Te_values = [3.0, 5.0, 7.0]                # eV (3 values, high)
E_field_values = [400, 500, 600, 800]      # V/cm (4 values, NEW!)
Tg_values = [570, 700]                     # K (2 values)
gamma_H_values = [0.001, 0.01, 0.05, 0.1]  # (4 values, contamination range)

# Fixed
L_diff = 0.1  # cm (from measurement)
```

**Total**: 4 × 3 × 4 × 2 × 4 = **384 combinations**

**Runtime**: ~30-60 minutes

### Reduced Sweep (Faster)

```python
ne_values = [5e8, 1e9]                    # cm⁻³ (2 values)
Te_values = [3.0, 5.0, 7.0]               # eV (3 values)
E_field_values = [400, 600, 800]          # V/cm (3 values)
Tg_values = [570, 700]                    # K (2 values)
gamma_H_values = [0.001, 0.01, 0.05]      # (3 values)
```

**Total**: 2 × 3 × 3 × 2 × 3 = **108 combinations**

**Runtime**: ~15-25 minutes

---

## Expected Results with Spatial Averaging

### Best Case:

Model matches all three species within factor 2-3:
- H:  within factor 2
- CH: within factor 2
- C₂: within factor 2

**Interpretation**:
- ✓ Chemistry validated
- ✓ Spatial averaging works
- ✓ 0-D model sufficient for CG

### Likely Case:

Model matches CH and C₂, H off by factor 5-10:
- CH: within factor 2 ✓
- C₂: within factor 2 ✓
- H:  factor 5-10 lower

**Interpretation**:
- ✓ Hydrocarbon chemistry validated
- ⚠️ H cycle needs refinement
- Check: H₂ build-up and dissociation
- Check: γ_H effective value

### Acceptable Case:

All species within factor 10:
- Within factor 10 is **very good** for plasma models!
- Absolute densities have uncertainties
- Ratios (H/CH, C₂/CH) more reliable

---

## Diagnostic Checks After Sweep

### 1. Check H₂ Density

If model predicts n_H₂ > 1e13 cm⁻³:
→ H₂ dissociation is major H source
→ This is GOOD (explains H density)

### 2. Check γ_H Sensitivity

If best fit has γ_H = 0.001:
→ Clean copper (H reflects)
→ Consistent with "decontaminated"

If best fit has γ_H = 0.05-0.1:
→ Contaminated surface
→ More H loss

### 3. Check E-field Sensitivity

If best fit has E = 400 V/cm:
→ Lower E in spatially-averaged CG
→ Less ionization

If best fit has E = 800 V/cm:
→ Higher E domina
tes in CG
→ More ionization

### 4. Check Ratios

**H/CH ratio**:
- Target: 3.1e7
- If matched: H production OK
- If off: Check H loss mechanisms

**C₂/CH ratio**:
- Target: 408
- If matched: C₂ chemistry OK
- If off: Check C₂ production/loss

---

## Implementation: Updated Sweep

I'll create two sweep scripts:
1. **Comprehensive** (384 runs, ~1 hour)
2. **Reduced** (108 runs, ~20 min)

Both include:
- ✓ E-field sweep (400-800 V/cm)
- ✓ Extended γ_H (0.001-0.1)
- ✓ H₂ recombination already included
- ✓ Spatial averaging interpretation

---

## Questions Answered:

1. **E-field tunable?** ✓ YES, now 400-800 V/cm
2. **H₂ measurements?** ✗ NO, but model includes H+H→H₂
3. **Cu contamination?** ~ γ_H = 0.001-0.1 (sweep covers range)
4. **Spatial averaging?** ✓ YES, 0-1 mm average (FAVORABLE!)

Ready to implement?
