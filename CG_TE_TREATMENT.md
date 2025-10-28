# Electron Temperature (Te) Treatment in Cathode Glow (CG)

## The Challenge

In the **Cathode Glow (CG) region**, the electron energy distribution function (EEDF) is:
- ❌ **NOT Maxwellian** (not thermal)
- ✓ **Non-thermal** with high-energy tail
- ⚠️ Te is **poorly defined** or undefined

### Why Non-Thermal?

1. **High E-field** near cathode accelerates electrons
2. **Short mean free path** → electrons don't thermalize
3. **Ionization zone** → continuous injection of cold electrons
4. **Result**: EEDF has:
   - **Bulk electrons**: ~0.5-1 eV (thermal-like)
   - **High-E tail**: 3-10+ eV (from acceleration)

---

## Impact on Rate Coefficients

Electron-impact reaction rates depend on EEDF:

```
k = ∫ σ(E) * v(E) * f(E) dE
```

Where:
- σ(E) = cross-section (energy-dependent)
- v(E) = electron velocity
- f(E) = EEDF (electron energy distribution)

For **non-thermal EEDF**:
- High-E tail contributes disproportionately to rates
- Threshold reactions (ionization, dissociation) enhanced
- Simple "effective Te" approximation may not work well

---

## Our Approaches

### 1. Effective Te (Simple, Current)

**Method**: Use single Te value that captures high-E tail contribution

**Range**: Te = 1.5 - 3 eV (higher than bulk ~0.5 eV)

**Pros**:
- ✓ Simple to implement
- ✓ Computationally fast
- ✓ Reasonable for 0-D model

**Cons**:
- ❌ Not physically rigorous
- ❌ May not capture all threshold effects
- ❌ Requires tuning against data

**Implementation**:
```python
params['Te'] = 1.5  # eV (effective)
# Use constant rates from define_rates()
```

---

### 2. Parameter Sweep (Pragmatic)

**Method**: Test Te = 0.5 → 10 eV, ne = 1e7 → 1e10 cm⁻³

**Goal**: Find best (Te, ne) combination that matches TALIF targets

**Pros**:
- ✓ Finds optimal effective Te empirically
- ✓ Validates model sensitivity
- ✓ Accounts for uncertainties in Te and ne

**Cons**:
- ⏱️ Computationally expensive (many runs)
- ⚠️ May find multiple local minima

**Implementation**: See `parameter_sweep_cg.py`

---

### 3. Two-Temperature Model (Advanced)

**Method**: Separate bulk and tail populations

**Model**:
```
f(E) = f_bulk * f_Maxwellian(T_bulk) + f_tail * f_tail(T_tail)
```

**Parameters**:
- T_bulk = 0.5 eV (bulk temperature)
- T_tail = 5 eV (tail temperature)
- f_tail = 0.1-0.3 (fraction in tail)

**Rate calculation**:
```python
k_eff = f_bulk * k(T_bulk) + f_tail * k(T_tail)
```

**Pros**:
- ✓ More physically realistic
- ✓ Captures threshold effects better
- ✓ Can be validated against EEDF measurements

**Cons**:
- ❌ Requires modifying define_rates()
- ❌ More parameters to tune (T_bulk, T_tail, f_tail)
- ⏱️ More complex

---

### 4. BOLSIG+ with Non-Thermal EEDF (Most Rigorous)

**Method**:
1. Use BOLSIG+ to solve Boltzmann equation
2. Input: E/N ratio for cathode region
3. Output: Full EEDF and rate coefficients

**Workflow**:
```
BOLSIG+ (E/N=100-200 Td) → EEDF(E) → <σv> for all reactions
```

**Pros**:
- ✓ Fully self-consistent
- ✓ No effective Te approximation needed
- ✓ Gold standard for accuracy

**Cons**:
- ❌ Requires external software (BOLSIG+)
- ❌ Must have cross-section database
- ❌ E-field may not be constant (need to solve Poisson)
- ⏱️ Significantly more complex

---

## Recommended Strategy

### Phase 1: Effective Te with Parameter Sweep ✓ (Current)

1. Run `parameter_sweep_cg.py`:
   - Test Te: 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0 eV
   - Test ne: 1e7, 5e7, 1e8, 5e8, 1e9, 5e9, 1e10 cm⁻³
   - Find best match to TALIF data

2. Analyze results:
   - If good match (within factor 2-5): **Done!**
   - If not: Proceed to Phase 2

### Phase 2: Two-Temperature Model (If Needed)

Modify `define_rates.py` to include:
```python
def define_rates_two_temp(params):
    T_bulk = 0.5  # eV
    T_tail = 5.0  # eV
    f_tail = 0.2  # 20% in tail

    # Weight electron-impact rates
    k['e_CH4_CH3_H'] = (1-f_tail)*k_bulk + f_tail*k_tail
    ...
```

### Phase 3: BOLSIG+ Integration (Future)

If very high accuracy needed:
1. Generate EEDF for E/N = 100-200 Td
2. Export <σv>(E/N) lookup table
3. Interpolate in simulation

---

## Practical Considerations for CG

### E-field Estimation

Cathode fall voltage: ~100-300 V
Cathode fall distance: ~0.1-0.3 cm
→ E-field: **500-1000 V/cm** (very high!)

But in glow discharge geometry:
- Negative glow has lower E-field: **50-150 V/cm**
- We use E = 80 V/cm as compromise

### ne Estimation

CG has **low ne** because:
- Diffusion losses to walls
- Low ionization rate (low E-field in glow)

Expected range: **1e7 - 1e9 cm⁻³**
(Much lower than bulk plasma ~1e10 cm⁻³)

### Diffusion Length

For H atoms at 400 K, P = 0.4 Torr:
```
D_H ~ 300 cm²/s
k_loss ~ 500 s⁻¹
L_diff = √(D/k) ~ 0.77 cm
```

Since L_discharge = 0.45 cm < L_diff:
- H escapes to walls before recombining
- Wall losses dominate (as in our model)

---

## Expected Te from Sweep

### Hypothesis

Based on physics, we expect:

**If sweep finds Te ~ 0.5-1.5 eV**:
- Bulk-dominated chemistry
- High-E tail not critical
- Suggests lower E-field than assumed

**If sweep finds Te ~ 2-5 eV**:
- High-E tail is important
- Ionization/dissociation enhanced
- Confirms non-thermal EEDF

**If sweep finds Te > 5 eV**:
- Very non-thermal
- May need two-temperature model
- Or: Cross-sections may be off

---

## How to Use Parameter Sweep Results

### 1. Run the sweep
```bash
python3 parameter_sweep_cg.py
```

### 2. Check output
- Best (Te, ne) combination
- Error for H, CH, C₂ individually
- Heatmaps showing parameter sensitivity

### 3. Interpret results

**If error is small (< 0.5 in log-scale)**:
→ Model is validated! ✓

**If H is off but CH/C₂ ratios are good**:
→ Check wall loss rates for H
→ Maybe increase Te slightly

**If all species are off by same factor**:
→ Total density may be wrong
→ Check pressure/temperature

**If CH << target but C₂ OK**:
→ CH consumption too fast
→ Check CH + CH4, CH + radicals rates

**If C₂ << target**:
→ C₂ production too slow
→ Check e + C2H2 → C2H rate (our new addition!)

---

## Next Steps After Sweep

1. **Document best parameters** in main.py
2. **Run full simulation** with best Te, ne
3. **Generate time evolution plots**
4. **If needed**: Implement two-temperature model
5. **Compare with SE region** (different Te, ne)

---

## References

- **Godyak & Piejak** (1990): EEDF in RF discharges
- **Lieberman & Lichtenberg** (2005): Plasma cathode sheaths
- **Hagelaar & Pitchford** (2005): BOLSIG+ solver
- **Luque (1997)**: Ar/CH₄ glow discharge (your reference paper!)

---

**Status**: Phase 1 (parameter sweep) ready to run
**Expected runtime**: ~5-15 minutes for 56 parameter combinations
