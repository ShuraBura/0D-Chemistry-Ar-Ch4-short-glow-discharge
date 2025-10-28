# CG Model - Constrained by Experimental Data

## Hard Constraints (from experimentalist)

| Parameter | Value | Justification |
|-----------|-------|---------------|
| **ne** | ≤ 1e9 cm⁻³ | Physical limit in CG region |
| **Te** | Can be high (3-10+ eV) | Non-thermal EEDF, high-E tail |
| **E-field** | ~500 V/cm | Revised from cathode fall estimate |
| **CG location** | 0-2 mm from cathode | Very close to surface |
| **Cathode material** | Copper | Affects wall recombination |
| **Tg** | 570 K (up to 800 K near cathode) | Higher than initial estimate! |
| **L_diff** | 0-1 mm (0.05-0.1 cm) | From measurement |

---

## Impact of Higher Tg

### Diffusion Scales with Temperature

```
D ∝ T^1.75 / P
```

At P = 0.4 Torr:

| Tg (K) | D_H (cm²/s) | k_loss (s⁻¹) for L=0.1cm | Lifetime (μs) |
|--------|-------------|---------------------------|---------------|
| 400    | 300         | 30,000                    | 33            |
| **570** | **452**     | **45,200**                | **22**        |
| **800** | **675**     | **67,500**                | **15**        |

**With Tg = 570 K, losses are 50% FASTER!**
**With Tg = 800 K, losses are 2.25× FASTER!**

This makes the challenge even harder...

---

## Revised Steady-State Calculation

### H Atom Balance with Constraints

Target: n_H = 8.57e15 cm⁻³

**Loss rate** (at Tg = 570 K, L_diff = 0.1 cm):
```
k_loss_H = 452 / (0.1)² = 45,200 s⁻¹
R_loss = k_loss × n_H = 45,200 × 8.57e15 = 3.9e20 cm⁻³/s
```

**Production rate needed**:
```
R_prod = k_e-CH4 × ne × n_CH4
       = 4.2e-11 × ne × 1.5e15
       = 6.3e4 × ne  (cm⁻³/s)

For R_prod = R_loss:
ne = 3.9e20 / 6.3e4 = 6.2e15 cm⁻³
```

**But gas density is only 9.66e15 cm⁻³, and ne ≤ 1e9 cm⁻³!**

### The Gap

With ne = 1e9 cm⁻³ (max allowed):
```
R_prod_single = 6.3e4 × 1e9 = 6.3e13 cm⁻³/s
R_loss = 3.9e20 cm⁻³/s

Gap = 3.9e20 / 6.3e13 = 6,200×  (!)
```

**Production is 6,000× too low with single channel!**

---

## Solutions to Close the Gap

### 1. Multiple H Production Channels ✓

Sum ALL reactions that produce H:

```python
# Electron-impact
e + CH4 → CH3 + H           (k = 4.2e-11, main)
e + CH4 → CH + H2 + H       (k = 0.7e-11)
e + H2 → H + H              (k = 6e-12, needs H2 buildup)
e + CH3 → CH2 + H           (k = 3e-11)

# ArStar (Penning)
ArStar + CH4 → products     (some produce H)
ArStar + H2 → H + H         (k = 6e-11)

# Neutral-neutral
CH + radicals → ... + H     (many channels)
CH2 + H → products
CH3 + H → products
```

If total H production is 5-10× higher than main channel:
→ Gap reduces to 600-1,200× (still huge!)

### 2. Copper Wall Recombination Coefficient ✓✓

**Key**: Copper has **low H recombination coefficient** γ_H

Literature values for γ_H on metals:
- **Stainless steel**: γ ~ 0.01-0.05 (low)
- **Copper (clean)**: γ ~ 0.001-0.01 (very low!)
- **Glass**: γ ~ 0.001 (very low)

**If γ_H = 0.01** on copper:
- 99% of H atoms **reflect** (return to gas)
- Only 1% **recombine** (lost as H₂)

**Effective loss rate**:
```
k_loss_eff = γ_H × k_loss_diffusion
           = 0.01 × 45,200
           = 452 s⁻¹  (!)
```

This **reduces gap by 100×!**

Now:
```
R_loss_eff = 452 × 8.57e15 = 3.9e18 cm⁻³/s
Gap = 3.9e18 / 6.3e13 = 62×
```

**Much better, but still need ~60× more production!**

### 3. High Te → More Ionization/Dissociation ✓✓

With Te = 5-7 eV (strong non-thermal EEDF):
- Ionization rates increase ~3-5×
- Dissociation rates increase ~2-3×
- More ArStar production → more Penning
- Higher ne within 1e9 limit

If Te scaling gives 3× boost:
→ Gap reduces to 20×

### 4. All Radicals Produce H ✓

CH, CH₂, CH₃, C₂H, C₂H₃, etc. all eventually decompose to H:
```
CH + CH → C2 + H2
CH + radicals → products + H
CH2 + H → CH + H2
CH3 + radicals → ... → H
```

If radical cascade produces 5× more H:
→ Gap reduces to 4×

### 5. Spatial Effects (0-D Limitation)

0-D model assumes **uniform** density.

Reality in CG (0-2 mm):
- **Density peaks** somewhere in region
- **Gradients** in ne, Te, n_species
- TALIF measures **integrated** or **averaged** density

If actual peak is 3× higher than average:
→ Model should predict n_H ~ 2.5e16 cm⁻³ (not 8.57e15)

---

## Revised Parameter Ranges

### Primary Parameters (for sweep)

```python
# Hard constraints
ne_values = [1e8, 3e8, 5e8, 8e8, 1e9]      # cm⁻³ (MAX 1e9!)
Te_values = [2.0, 3.0, 4.0, 5.0, 7.0, 10.0] # eV (high!)

# From experiment
E_field_values = [400, 500, 600]            # V/cm
Tg_values = [570, 650, 800]                 # K
L_diff_values = [0.05, 0.075, 0.1]          # cm

# Wall recombination (critical!)
gamma_H_values = [0.001, 0.01, 0.05, 0.1]   # on copper
```

**Total**: 5 × 6 × 3 × 3 × 3 × 4 = **3,240 combinations**

This is too many! Need to reduce...

### Reduced Sweep (Realistic)

```python
# Focus on most likely values
ne_values = [3e8, 5e8, 8e8, 1e9]           # 4 values
Te_values = [3.0, 5.0, 7.0]                # 3 values (high)
E_field = 500                              # Fixed at 500 V/cm
Tg_values = [570, 700]                     # 2 values
L_diff = 0.1                               # Fixed at 0.1 cm
gamma_H_values = [0.001, 0.01, 0.05]       # 3 values (copper)
```

**Total**: 4 × 3 × 2 × 3 = **72 combinations** ✓

Runtime: ~10-20 minutes

---

## Implementation: Wall Recombination

### Add to define_rates_tunable.py

```python
def define_rates_tunable(params):
    # ... existing code ...

    # Wall recombination coefficient
    gamma_H = params.get('gamma_H', 0.01)  # Default for copper

    # Effective wall loss = gamma × (D / L_diff²)
    # Only a fraction gamma actually sticks and is lost

    for species_name, D_400K in D_ref.items():
        D = D_400K * D_scale
        k_loss_diffusion = D / (L_diff ** 2)

        # Apply recombination coefficient for H
        if species_name == 'H':
            k_loss_effective = gamma_H * k_loss_diffusion
        else:
            k_loss_effective = k_loss_diffusion  # Full loss for other species

        # Map to reaction names
        if species_name in loss_map and loss_map[species_name] is not None:
            k[loss_map[species_name]] = k_loss_effective
```

This implements the physical picture:
- H atoms hit wall with rate = D/L²
- Fraction γ recombines (lost as H₂)
- Fraction (1-γ) reflects back to gas

---

## Expected Results

### With γ_H = 0.01 (realistic for copper)

**Effective k_loss_H**:
```
k_loss_H_eff = 0.01 × 45,200 = 452 s⁻¹
```

**Steady-state with ne = 1e9, Te = 5 eV**:
```
R_prod ~ 3× (Te boost) × 6.3e13 = 1.9e14 cm⁻³/s
n_H_ss = R_prod / k_loss = 1.9e14 / 452 = 4.2e11 cm⁻³
```

Still **20,000× too low!**

**Need**: Additional production channels (10×) + high Te (3×) + spatial effects (3×)
→ 4.2e11 × 10 × 3 = 1.3e13 cm⁻³

Still 650× too low...

### Conclusion

Even with:
- γ_H = 0.01 (optimistic for copper)
- Te = 7 eV (very high)
- ne = 1e9 cm⁻³ (maximum)
- All H production channels

**Model may still underpredict H by factor 10-100.**

**Possible explanations**:
1. **0-D model insufficient** - Need spatial transport
2. **L_diff interpretation** - Not same as diffusion-loss length?
3. **Additional H sources** - Surface desorption, plasma-wall interaction?
4. **Measurement uncertainty** - TALIF calibration?
5. **H₂ dissociation underestimated** - Major H source if H₂ builds up?

---

## Recommended Strategy

### Phase 1: Focused Sweep (NOW)
Run 72-combination sweep to find:
- Best (ne, Te, Tg, γ_H) within constraints
- Which species match (CH, C₂ may be easier than H)
- Sensitivity to each parameter

### Phase 2: Analyze Discrepancies
If H is still off by >10×:
- Check H₂ buildup → H₂ dissociation important?
- Check CH/C₂ ratios → chemistry OK but H lost?
- Check wall interactions → H absorption/desorption?

### Phase 3: Extended Model (if needed)
- Add surface reactions (H + wall ⇌ H_absorbed)
- Add H₂ vibrational states (v=1-4)
- Consider 1-D spatial transport

---

## Questions?

1. **Is H₂ density measured?** (Important H source)
2. **Surface condition of copper?** (Clean, oxidized, sputtered? Affects γ)
3. **How is TALIF calibrated?** (Could affect absolute densities)
4. **Are CH and C₂ easier to match than H?** (Would suggest H-specific issue)

---

Let me now implement the copper wall recombination model and run the focused sweep!
