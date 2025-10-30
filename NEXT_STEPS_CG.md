# Next Steps for CG Model Validation

## Summary of Current Status

### What We've Done ✓

1. **Implemented chemistry improvements from expert audit**:
   - Fixed electron-impact rates (Janev-Reiter 2002)
   - Corrected neutral-neutral rates (Baulch 2005)
   - Added missing reactions (7 new reactions)
   - Total: 244 reactions

2. **Incorporated your experimental data**:
   - Cathode fall: 253 V / 4 mm = **632.5 V/cm**
   - Diffusion length: **0-1 mm** (0-0.1 cm)
   - TALIF targets for H, CH, C₂

3. **Created tunable model**:
   - `define_rates_tunable.py`: Parameters scale with Te, ne, L_diff, E-field
   - Wall losses: k = D / L_diff²
   - Ion drift: k = μ × E / L

4. **Documentation**:
   - `CG_PHYSICS_ANALYSIS.md`: Physics constraints
   - `CG_TE_TREATMENT.md`: Te handling strategies
   - `EXPERIMENTAL_TARGETS.md`: CG vs SE analysis

---

## Critical Findings

### Your Measured L_diff = 0.1 cm Changes Everything!

**Wall loss rates scale as 1/L²**:

| Species | Old k_loss (s⁻¹) | New k_loss (s⁻¹) | Factor |
|---------|------------------|-------------------|--------|
| H       | 800              | **30,000**        | 37×    |
| CH      | 800              | **15,000**        | 19×    |
| C₂      | 200              | **10,000**        | 50×    |

**Ion drift rates with E = 800 V/cm**:

| Ion | μ (cm²/V/s) | k_drift (s⁻¹) | Lifetime |
|-----|-------------|---------------|----------|
| Ar⁺ | 3057        | **5.4×10⁶**   | 0.18 μs  |
| CH₃⁺| 4950        | **8.8×10⁶**   | 0.11 μs  |

**Ions are lost in < 1 μs!**

---

## The Challenge

With such fast losses, steady-state densities will be **low** unless:

1. **ne is high** (1e10-1e11 cm⁻³) → more production
2. **Te is high** (3-7 eV) → more ionization/dissociation
3. **Wall recombination** returns H as H₂ (partial recycling)
4. **Spatial effects** not captured by 0-D model

### Example Calculation

Target: n_H = 8.57e15 cm⁻³

**Production rate**:
```
R_prod = k × ne × n_CH4
       = 4.2e-11 × ne × 1.5e15  (cm⁻³/s)
```

**Loss rate**:
```
R_loss = k_loss × n_H
       = 3e4 × n_H  (cm⁻³/s)
```

**Steady-state** (R_prod = R_loss):
```
n_H = (4.2e-11 × ne × 1.5e15) / 3e4
    = 2.1e-6 × ne

For n_H = 8.57e15:  ne = 4.1e21 cm⁻³ (!)
```

**This is impossible!** Plasma density cannot exceed gas density (9.66e15 cm⁻³).

**Conclusion**: The model is missing something, OR parameters need careful tuning.

---

## Possible Solutions

### Option 1: Multi-Channel H Production

H is produced by many reactions:
```
e + CH4 → CH3 + H
e + CH4 → CH + H2 + H
e + H2 → H + H
ArStar + CH4 → products (some H)
CH + radicals → products + H
...
```

**Total production** may be 5-10× higher than single channel.

### Option 2: Partial Wall Return

If walls recombine H → H₂ with efficiency γ < 1:
- Some H sticks permanently
- Some H₂ returns and is re-dissociated
- **Effective** k_loss < D/L²

Maybe: k_loss_eff = γ × D/L²  where γ = 0.1-0.5

### Option 3: Spatial Averaging

TALIF measures **column-integrated** or **spatially-averaged** densities.

Your 0-D model represents **single point** in CG.

If CG has spatial gradients:
- Center (peak): n_H = 1e16
- Edge: n_H = 1e14
- **Average**: n_H = 8.57e15 ✓

0-D model should match **peak** or **average**?

### Option 4: E-Field Lower Than Assumed

If E-field in CG is actually **300-500 V/cm** (not 800):
- Ion drift reduced 2-3×
- Allows higher ion densities
- More chemistry possible

Cathode fall E-field **non-uniform**:
- Peak at cathode surface: 1000+ V/cm
- CG region (a few mm away): 300-500 V/cm?

---

## Recommended Next Steps

### Phase 1: Parameter Sweep (ESSENTIAL)

Run comprehensive sweep to find viable (Te, ne, L_diff, E) combinations:

```python
# Parameter ranges
Te_values = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0]  # eV
ne_values = [1e8, 5e8, 1e9, 5e9, 1e10, 5e10]    # cm⁻³
L_diff_values = [0.05, 0.1, 0.15, 0.2, 0.3]      # cm
E_field_values = [400, 600, 800, 1000]           # V/cm

# Optional: Wall return coefficient
wall_return_H = [0.0, 0.3, 0.5, 0.7]
```

**Total combinations**: 7 × 6 × 5 × 4 = **840 runs**
**Runtime**: ~30-60 minutes

**Goal**: Find parameter set that matches H, CH, C₂ within factor 2-5.

### Phase 2: Sensitivity Analysis

For best-fit parameters, test:
1. **Rate scaling**: Multiply all e-impact rates by 0.5-2×
2. **Wall return**: Test γ = 0-1 for H return fraction
3. **Multi-channel production**: Sum all H-producing reactions

### Phase 3: Physical Validation

Compare best-fit parameters with:
1. **ne measurement** (if available): Langmuir probe, interferometry?
2. **Te estimate**: From EEDF measurement or Ar line ratios?
3. **E-field profile**: Measure or calculate from Poisson equation?
4. **Spatial profiles**: TALIF scan vs position?

---

## Implementation Plan

### Step 1: Create Parameter Sweep Script

Modify `parameter_sweep_cg.py` to include:
- L_diff sweep
- E_field sweep
- Wall return coefficient
- Multi-dimensional output

```python
def run_comprehensive_sweep():
    """Multi-dimensional parameter sweep."""
    results = []

    for Te in Te_values:
        for ne in ne_values:
            for L_diff in L_diff_values:
                for E_field in E_field_values:
                    params = build_params(Te, ne, L_diff, E_field)
                    result = run_simulation(params)
                    results.append(result)

    # Find best fit
    best = min(results, key=lambda r: r['error'])
    return best, results
```

### Step 2: Run and Analyze

```bash
python3 comprehensive_sweep_cg.py
```

Expected output:
- Best (Te, ne, L_diff, E) combination
- Heatmaps showing parameter sensitivity
- Error decomposition (H, CH, C₂ individually)

### Step 3: Validate and Iterate

If best fit is **good** (error < factor 3):
→ Model validated! ✓

If best fit is **poor** (error > factor 10):
→ Need to add physics:
- Wall recycling
- Spatial transport
- Two-temperature EEDF

---

## Questions for Experimental Validation

To constrain the model further, please provide (if available):

1. **Electron density (ne)**:
   - Measurement method?
   - Value and uncertainty?
   - Where measured (CG center, edge)?

2. **Electron temperature (Te)**:
   - Is EEDF measured?
   - Ar line ratios for Te estimate?

3. **Gas temperature (Tg)**:
   - Rotational temperature from OH/CH/N₂?
   - Thermocouple measurement?

4. **E-field profile**:
   - Measured vs position?
   - Or: voltage probe measurements?

5. **Spatial profiles**:
   - TALIF scan vs position?
   - How is "CG" region defined?
   - Peak density or average?

6. **Wall material**:
   - Stainless steel, glass, polymer?
   - Surface condition (clean, coated, sputtered)?

7. **L_diff measurement**:
   - How measured (TALIF intensity decay)?
   - Uncertainty?

---

## Deliverables

### Code Ready to Run:
1. ✓ `define_rates_tunable.py` - Tunable rate module
2. ✓ `run_cg_optimized.py` - Single-point CG simulation
3. ⏳ `comprehensive_sweep_cg.py` - Multi-D parameter sweep (need to create)

### Documentation:
1. ✓ `CG_PHYSICS_ANALYSIS.md` - Physics constraints
2. ✓ `CG_TE_TREATMENT.md` - Te treatment options
3. ✓ `EXPERIMENTAL_TARGETS.md` - CG vs SE comparison
4. ✓ `IMPROVEMENTS_SUMMARY.md` - Chemistry improvements
5. ✓ `NEXT_STEPS_CG.md` - This file

---

## Timeline Estimate

| Task | Time | Priority |
|------|------|----------|
| Create comprehensive sweep | 1-2 hr | HIGH |
| Run 840 simulations | 1-2 hr | HIGH |
| Analyze results | 1 hr | HIGH |
| If needed: Add wall return | 2 hr | MEDIUM |
| If needed: Two-temperature Te | 3 hr | MEDIUM |
| Validate against experiments | Ongoing | HIGH |

**Total**: ~5-10 hours of computation + analysis

---

## Expected Outcomes

### Best Case:
Parameter sweep finds (Te, ne, L_diff, E) that matches TALIF ± factor 2
→ **Model validated!**

### Likely Case:
Parameter sweep finds partial match (some species good, others off)
→ Identify which physics is missing
→ Add wall recycling or multi-temperature EEDF

### Worst Case:
No parameter set matches data
→ 0-D model insufficient (need 1-D or 2-D transport)
→ Or: Reaction rates significantly wrong

---

## Shall I proceed with creating the comprehensive sweep tool?

This will be the key to understanding if your model can match the CG data!
