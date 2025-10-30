# Cathode Glow (CG) Physics - Experimental Constraints

## Experimental Data (Your Measurements)

### Cathode Fall
- **Voltage**: 253 V
- **Distance**: 4 mm = 0.4 cm
- **Average E-field**: 253 V / 0.4 cm = **632.5 V/cm**

### Spatial Considerations
- **CG location**: Close to cathode (within cathode fall region)
- **E-field in CG**: Likely **higher than average** (700-1000 V/cm)
- **Diffusion length**: 0-1 mm = **0-0.1 cm** (very short!)

---

## Implications

### 1. E-Field is VERY High

Previous assumption: E = 80 V/cm ❌
**Correct value**: E = **700-1000 V/cm** ✓

This affects:
- Electron acceleration (high-E tail enhanced)
- Ion drift losses (much faster)
- EEDF highly non-thermal

**Ion drift rate**:
```
k_drift = μ * E / L
```

For ArPlus (μ = 3057 cm²/V/s):
```
k_drift = 3057 * 800 / 0.45 = 5.4e6 s⁻¹ (!)
```

This is **huge** - ions are swept out in microseconds!

---

### 2. Diffusion Length is VERY Short

**Measured**: L_diff = 0-1 mm = 0-0.1 cm

This dramatically changes wall loss rates!

**Wall loss formula**:
```
k_loss = D / L_eff² = D / (L_diff)²
```

For different species at P = 0.4 Torr, T = 400 K:

| Species | D (cm²/s) | L_diff (cm) | k_loss (s⁻¹) |
|---------|-----------|-------------|--------------|
| **H**   | ~300      | 0.1         | **30,000**   |
| **CH**  | ~150      | 0.1         | **15,000**   |
| **C₂**  | ~100      | 0.1         | **10,000**   |
| **CH₃** | ~100      | 0.1         | **10,000**   |
| **CH₄** | ~80       | 0.1         | **8,000**    |

**These are 10-100× higher than our current values!**

Current (from audit): 100-800 s⁻¹ ❌
**Correct (your data)**: 8,000-30,000 s⁻¹ ✓

---

### 3. Te Must Be High (Non-Thermal EEDF)

With E/N = 800 V/cm / (9.66e15 cm⁻³) ≈ **830 Td**

This is **very high** E/N ratio → strongly non-thermal EEDF

Expected:
- Bulk electrons: ~0.5 eV (thermalized)
- High-E tail: 5-10+ eV (from acceleration)
- Effective Te: **2-5 eV** (much higher than bulk)

---

## Updated Model Parameters for CG

### Baseline Configuration

```python
params_CG = {
    'P': 0.4,              # Torr (measured)
    'Tg': 400,             # K (estimated, could be 300-500 K)
    'ne': 1e8,             # cm⁻³ (tunable: 1e7-1e9)
    'Te': 3.0,             # eV (effective, tunable: 1-7 eV)
    'E_field': 800,        # V/cm (from cathode fall: 632-1000 V/cm)
    'L_discharge': 0.45,   # cm (total gap)
    'L_diff': 0.1,         # cm (YOUR DATA: 0-1 mm)
}
```

### Tunable Ranges (for Parameter Sweep)

| Parameter | Min | Baseline | Max | Notes |
|-----------|-----|----------|-----|-------|
| **ne** | 1e7 | 5e8 | 5e9 | Low in CG, measure if possible |
| **Te** | 1.0 | 3.0 | 7.0 | Effective for non-thermal EEDF |
| **E_field** | 600 | 800 | 1000 | From cathode fall measurement |
| **Tg** | 300 | 400 | 500 | Gas heating in discharge |
| **L_diff** | 0.05 | 0.1 | 0.15 | Your measurement: 0-1 mm |

---

## Recalculated Wall Loss Rates

### Method

For species i:
```
k_loss,i = D_i / L_diff²
```

Where diffusion coefficient at P = 0.4 Torr, T = 400 K:
```
D_i = D_i,ref * (T/T_ref)^1.75 * (P_ref/P)
```

### Updated Loss Rates (L_diff = 0.1 cm)

```python
# Group 11: Loss Reactions (for L_diff ~ 0.1 cm)
k['loss_H_11_1'] = 3.0e4      # H: very fast diffusion
k['loss_CH_11_9'] = 1.5e4     # CH: fast, reactive
k['loss_CH3_11_21'] = 1.0e4   # CH3: moderate
k['loss_CH2_11_1'] = 1.2e4    # CH2: fast
k['loss_C2_11_3'] = 8.0e3     # C2: slower
k['loss_C_11_8'] = 1.2e4      # C: fast
k['loss_C2H_11_12'] = 8.0e3   # C2H: slower
k['loss_C2H2_11_19'] = 5.0e3  # C2H2: slow
k['loss_C2H4_11_20'] = 5.0e3  # C2H4: slow
k['loss_C2H6_11_5'] = 4.0e3   # C2H6: slow
k['loss_CH4_11_6'] = 4.0e3    # CH4: slow
k['loss_H2_11_2'] = 8.0e3     # H2: fast
k['loss_Ar_11_7'] = 3.0e3     # Ar: moderate
k['loss_e_11_4'] = 5.0e4      # e: very fast (ambipolar)
```

**These are 10-50× higher than previous "improved" values!**

---

## Ion Drift Losses (CRITICAL)

With E = 800 V/cm, L = 0.45 cm:

```python
k['drift_ArPlus_10_1'] = 3057 * 800 / 0.45 = 5.4e6 s⁻¹
k['drift_CH4Plus_10_2'] = 6432 * 800 / 0.45 = 1.1e7 s⁻¹
k['drift_CH3Plus_10_3'] = 4950 * 800 / 0.45 = 8.8e6 s⁻¹
k['drift_CH5Plus_10_4'] = 4762 * 800 / 0.45 = 8.5e6 s⁻¹
```

**Ions are lost in ~0.1-0.2 μs!** This is extremely fast.

---

## Physical Interpretation

### Why is CG Chemistry Possible?

Despite fast losses, chemistry happens because:

1. **Production rates are also fast**
   - e + CH4 → products: k × ne × n_CH4
   - With ne = 1e8, n_CH4 = 1.5e15: rate ~ 1e8 × 1e15 × 1e-11 = **1e12 s⁻¹**
   - This is much faster than diffusion loss!

2. **Steady-state balance**
   ```
   Production = Loss
   k_reaction × ne × n_reactant = k_loss × n_product
   ```

3. **Recycling at walls**
   - H → wall → H2 (returns to gas)
   - Some species stick, others recombine and return

### H Atom Balance

With k_loss_H = 30,000 s⁻¹, lifetime τ = 33 μs

But H is produced continuously:
```
Production: e + CH4 → CH3 + H (k = 4.2e-11 cm³/s)
           Rate = 4.2e-11 × 1e8 × 1.5e15 = 6.3e12 s⁻¹ (per cm³)

Steady-state: n_H = Production / k_loss
                  = 6.3e12 / 3e4 = 2.1e8 cm⁻³
```

**This is 40,000× LOWER than TALIF target (8.57e15 cm⁻³)!**

### Problem!

With your measured L_diff = 0.1 cm, we get:
- **k_loss too high** → densities too low
- **Unless**: Production rates are underestimated

**Possible explanations**:
1. **ne is higher** than 1e8 (need 1e10-1e11 cm⁻³?)
2. **Te is higher** (more ionization/dissociation)
3. **Wall losses partially "return"** H via recombination
4. **L_diff measured differently** than model L_eff
5. **Spatial averaging** - TALIF sees integrated density

---

## Recommended Approach: Multi-Parameter Sweep

Given uncertainties, we need to sweep:

### Primary Parameters (High Impact)
1. **Te**: 1-7 eV (non-thermal EEDF)
2. **ne**: 1e7-1e11 cm⁻³ (wide range!)
3. **L_diff**: 0.05-0.3 cm (affects all wall losses)

### Secondary Parameters (Medium Impact)
4. **E_field**: 600-1000 V/cm (from cathode fall)
5. **Tg**: 300-500 K (gas heating)

### Rate Uncertainties (to explore)
6. **Scale electron-impact rates**: 0.5-2× (cross-section uncertainties)
7. **Wall recombination coefficient**: 0-1 (H return fraction)

---

## Proposed Tunable Model Structure

```python
class TunableCGModel:
    """CG model with tunable parameters."""

    def __init__(self):
        # Fixed (measured)
        self.P = 0.4  # Torr
        self.V_cathode = 253  # V
        self.d_cathode = 0.4  # cm
        self.L_gap = 0.45  # cm

        # Tunable - Primary
        self.Te = 3.0           # eV (1-7)
        self.ne = 5e8           # cm⁻³ (1e7-1e11)
        self.L_diff = 0.1       # cm (0.05-0.3)

        # Tunable - Secondary
        self.E_field = 800      # V/cm (600-1000)
        self.Tg = 400           # K (300-500)

        # Tunable - Rate scaling
        self.scale_e_impact = 1.0      # 0.5-2.0
        self.wall_return_H = 0.0       # 0-1 (H recombination at walls)

    def update_wall_losses(self):
        """Recalculate wall losses based on L_diff."""
        # D ∝ T^1.75 / P
        D_ref = {'H': 300, 'CH': 150, 'C2': 100, ...}  # at 400 K, 0.4 Torr

        for species, D_400K in D_ref.items():
            D = D_400K * (self.Tg / 400)**1.75
            k_loss = D / self.L_diff**2
            # Store in rates dict

    def run_simulation(self):
        """Run with current parameters."""
        # Build rates with scaling
        # Run solve_ivp
        # Return (H, CH, C2) densities
```

---

## Next Steps

### 1. Update define_rates.py with Tunability

Add parameters for:
- L_diff (for wall losses)
- E_field scaling
- Rate scaling factors

### 2. Run Comprehensive Parameter Sweep

Sweep grid:
- Te: [1, 2, 3, 4, 5, 6, 7] eV (7 values)
- ne: [1e7, 5e7, 1e8, 5e8, 1e9, 5e9, 1e10, 5e10, 1e11] (9 values)
- L_diff: [0.05, 0.1, 0.15, 0.2] cm (4 values)

Total: **252 simulations** (~10-30 min on single core)

### 3. Analyze Sensitivity

Which parameter has most impact on:
- H density
- CH density
- C2 density
- H/CH ratio
- C2/CH ratio

### 4. Validate Against TALIF

Find parameter set that minimizes:
```
error = |log(H/H_target)| + |log(CH/CH_target)| + |log(C2/C2_target)|
```

---

## Key Questions for You

1. **Can ne be measured?** (e.g., via Langmuir probe, microwave interferometry?)
   - This would constrain the model significantly

2. **Is Tg known?** Or can it be estimated from rotational temperature?

3. **Wall material?** (stainless steel, glass, polymer?)
   - Affects recombination coefficient

4. **Exactly where is "CG"?** Distance from cathode?
   - Affects E-field and ne estimates

5. **How is L_diff = 0-1 mm measured?**
   - Spatial TALIF scans?
   - Intensity decay length?
   - This affects interpretation

---

## Immediate Action Items

1. ✓ Update define_rates.py with:
   - E_field = 800 V/cm (not 80!)
   - L_diff = 0.1 cm parameter
   - Recalculate all wall losses: k = D/L_diff²

2. ✓ Create tunable model class

3. ✓ Run parameter sweep with realistic ranges

4. ✓ Document sensitivity analysis

Should I proceed with implementing these updates?
