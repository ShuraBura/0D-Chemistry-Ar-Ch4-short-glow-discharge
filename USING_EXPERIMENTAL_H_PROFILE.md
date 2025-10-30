# Using Experimental nH(x) Profile in 0-D Model

## Overview

You have an experimental spatial profile of H density: **nH(x)** measured over 0-1 mm from the cathode. Since your model is **0-dimensional** (predicts spatially uniform density), you need to map this profile to a single representative value.

## Quick Answer

**RECOMMENDED APPROACH**: Use the **spatial average** of your nH(x) profile:

```python
import numpy as np

# Your experimental data
x = np.array([...])    # Position (cm)
nH = np.array([...])   # H density (cm⁻³)

# Calculate spatial average
nH_avg = np.trapz(nH, x) / (x[-1] - x[0])

# Use in simulation
set_density('H', nH_avg)
```

**Why?** Your TALIF measurements are already spatially averaged over the 0-1 mm region (see SPATIAL_AVERAGING_AND_H2.md), so this is the most physically consistent comparison.

---

## Four Options for Using nH(x)

### Option 1: Spatial Average as Initial Condition ✓ RECOMMENDED

**What it does**: Use spatial average as starting point, let H evolve dynamically

**How to implement**:
```python
# In your simulation script (e.g., run_cg_optimized.py)
from load_experimental_profile import use_profile_in_simulation

# Load and process profile
nH_avg, L_diff = use_profile_in_simulation('your_profile.csv', method='average')

# Set as initial condition
def calculate_initial_densities(y0, species_list, params):
    def set_density(name, value):
        y0[species_list.index(name)] = value

    # ... other species ...
    set_density('H', nH_avg)  # Use experimental average
```

**When to use**:
- Default approach for model validation
- You want to check if model can maintain/reach experimental H density
- You want to validate H production/loss mechanisms

**Pros**:
- Tests model self-consistency
- Validates chemistry and transport
- Most physically meaningful

**Cons**:
- Model may not maintain experimental H (indicates missing physics)

---

### Option 2: Fixed H Density (Non-Evolving)

**What it does**: Prescribe H density, keep it constant throughout simulation

**How to implement**:
```python
# Modify odefun.py or odefun_optimized.py
def dydt_func(self, t, y):
    # ... calculate reactions ...

    # Fix H density
    H_idx = self.species.index('H')
    self.dydt[H_idx] = 0  # Prevent H from changing

    # Also fix other background species
    self.dydt[self.e_idx] = 0
    self.dydt[self.Ar_idx] = 0
    self.dydt[self.CH4_idx] = 0

    return self.dydt
```

Or use the provided script:
```bash
python3 run_with_fixed_H.py
```

**When to use**:
- You want to treat H as a known background species
- Focus on validating CH and C₂ chemistry only
- H density is well-measured and trusted

**Pros**:
- H matches experiment by design
- Focuses validation on other species
- Faster convergence (one less differential equation)

**Cons**:
- **Doesn't validate H chemistry** (decoupled)
- May hide issues with H production/loss
- Less physically self-consistent

---

### Option 3: Tune Drift Gain to Match Profile

**What it does**: Adjust phenomenological H source term to match experimental density

**How to implement**:
```python
# Calculate required drift gain
nH_target = nH_avg
k_loss_H = 30000  # s⁻¹ (from L_diff = 0.1 cm)

# At steady state: Production + Drift = Loss
R_loss = k_loss_H * nH_target
R_prod_chem = 1e17  # Estimate from test run (cm⁻³/s)

H_drift_gain = R_loss - R_prod_chem

# Update in odefun.py line 18
self.H_drift_gain = H_drift_gain  # Instead of 3.2e17
```

**When to use**:
- You accept that model needs a phenomenological source term
- You want to match H while still allowing dynamics
- Intermediate between fully dynamic and fully fixed

**Pros**:
- Can match experimental H density
- Allows H to respond to chemistry changes
- Pragmatic for engineering predictions

**Cons**:
- Phenomenological (not first-principles)
- Drift gain value may not be physically meaningful
- Hides gaps in fundamental chemistry model

---

### Option 4: Estimate Diffusion Length from Profile

**What it does**: Use profile shape to extract diffusion length parameter

**How to implement**:
```python
from load_experimental_profile import estimate_diffusion_length

# Load profile
x, nH = load_nH_profile('your_profile.csv')

# Fit exponential decay: nH(x) ~ exp(-x/L_diff)
L_diff_measured = estimate_diffusion_length(x, nH)

print(f"Measured L_diff: {L_diff_measured:.4f} cm")

# Use in define_rates_tunable.py
params['L_diff'] = L_diff_measured
```

**When to use**:
- Profile shows clear exponential decay
- You want to extract transport parameters from data
- Complement to density value extraction

**Pros**:
- Extracts additional physics (diffusion) from profile
- Improves transport model fidelity
- Can validate/calibrate L_diff parameter

**Cons**:
- Requires clear exponential decay (may not fit all profiles)
- Sensitive to noise in tail of profile
- 0-D model doesn't explicitly model spatial diffusion

---

## Comparison Table

| Option | H Evolves? | Validates H? | Complexity | Best For |
|--------|------------|--------------|------------|----------|
| 1. Avg as IC | Yes | **Yes** | Low | **Model validation** |
| 2. Fixed H | No | No | Low | Validating other species |
| 3. Tune drift | Yes | Partial | Medium | Pragmatic predictions |
| 4. Extract L_diff | Yes | Yes | Medium | Transport validation |

---

## Recommended Workflow

### Step 1: Prepare your experimental data

Create a CSV file with your nH(x) profile:
```csv
x_cm,nH_cm3
0.000,5.0e15
0.001,6.2e15
0.002,7.8e15
...
0.010,4.5e15
```

### Step 2: Process profile

```bash
python3 load_experimental_profile.py
```

This will calculate:
- Spatial average: `<nH>` (use in model)
- Peak value: `nH_peak` (for reference)
- Estimated L_diff from decay

### Step 3: Run simulation (Option 1)

```python
# In your simulation script
nH_avg = 8.58e15  # From spatial average
set_density('H', nH_avg)

# Run simulation
# ...

# Check if H stays near experimental value
if abs(final_H - nH_avg) / nH_avg < 0.3:
    print("✓ Model maintains H within 30%")
else:
    print("✗ Model drifts from experimental H")
```

### Step 4: If needed, try Option 2 (Fixed H)

```bash
python3 run_with_fixed_H.py
```

This focuses validation on CH and C₂ chemistry.

### Step 5: Interpret results

**If Option 1 works** (H stays near experimental):
→ ✓ Model is self-consistent!
→ H chemistry and transport are validated

**If Option 1 fails** (H drifts significantly):
→ Try Option 3 (tune drift gain)
→ Or accept that 0-D model has limitations
→ Consider adding missing physics (wall recycling, etc.)

**If Option 2 shows good CH/C₂ match**:
→ Hydrocarbon chemistry is validated
→ H-specific mechanisms need refinement

---

## Physical Interpretation

### What does spatial averaging mean?

Your 0-D model predicts:
```
n(x) = n₀  (uniform, everywhere)
```

Reality:
```
n(x) = n(x)  (varies with position)
```

TALIF measures:
```
<n> = (1/L) ∫₀ᴸ n(x) dx
```

**Therefore**: Compare model `n₀` to measured `<n>`, NOT to `n_peak`

### Example:

If your profile is:
- Peak: nH = 1.0e16 cm⁻³ at x = 0.3 mm
- Spatial average: <nH> = 8.58e15 cm⁻³
- Model predicts: nH = 8.5e15 cm⁻³

**Interpretation**: ✓ Excellent match! (within 1%)

Don't expect model to match the peak (1e16) because:
1. 0-D model is spatially uniform
2. TALIF measures the average, not the peak

---

## Code Examples

### Example 1: Load and use profile

```python
from load_experimental_profile import use_profile_in_simulation

# Process your profile
nH_fixed, L_diff = use_profile_in_simulation(
    'my_nH_profile.csv',
    method='average'  # Options: 'average', 'peak', 'rms'
)

# Use in simulation
params = {
    'H_fixed': nH_fixed,
    'L_diff': L_diff,
    # ... other parameters ...
}
```

### Example 2: Create profile from discrete points

```python
import numpy as np

# Your measurement points
x_mm = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]  # mm
nH_measured = [5e15, 7e15, 9e15, 8e15, 7e15, 6e15]  # cm⁻³

# Convert to cm
x_cm = np.array(x_mm) / 10.0

# Calculate average
nH_avg = np.trapz(nH_measured, x_cm) / (x_cm[-1] - x_cm[0])

print(f"Spatial average H density: {nH_avg:.2e} cm⁻³")
```

### Example 3: Modify existing simulation for fixed H

```python
# In odefun_optimized.py, add after line 109:

    def __call__(self, t, y):
        # ... existing code ...

        # NEW: Fix H density if requested
        if 'H_fixed' in self.params:
            H_idx = self.species.index('H')
            y[H_idx] = self.params['H_fixed']  # Enforce value
            # ... calculate reactions ...
            self.dydt[H_idx] = 0  # Prevent evolution

        return self.dydt
```

---

## FAQ

### Q1: Should I use peak or average from my profile?

**A: Use average.** Your TALIF measurement is spatially averaged over 0-1 mm. The 0-D model represents this averaged region, not the peak.

### Q2: What if my model H differs from experimental by factor 5-10?

**A:** This is common! Options:
1. Run parameter sweep to find better Te, ne, E-field
2. Check if H₂ recycling is properly included
3. Consider adding wall return physics (gamma_H tuning)
4. Accept limitation and use fixed H (Option 2)

### Q3: Can I use a 1-D profile directly in the 0-D model?

**A:** Not directly. The 0-D model fundamentally predicts uniform density. You must:
- Reduce profile to single value (average, peak, etc.)
- Or upgrade to 1-D model (major code change)

### Q4: Does fixing H invalidate the model?

**A:** Partially. Fixing H:
- ✓ Still validates CH and C₂ chemistry
- ✗ Doesn't validate H production/loss
- Use when: H is well-known, focus is on hydrocarbons

### Q5: My profile doesn't decay exponentially. Can I still use it?

**A:** Yes! You can still:
- Calculate spatial average (always valid)
- Use peak value (conservative)
- Use RMS average (weights high-density regions)

The exponential fit (Option 4) is optional for extracting L_diff.

---

## Summary

**For model validation** → Option 1 (spatial average as IC)
**For practical predictions** → Option 2 (fixed H) or Option 3 (tuned drift)
**For transport studies** → Option 4 (extract L_diff)

The key insight: **0-D model = spatially averaged reality**

Compare:
- Model output `n_species` ↔ Measured `<n_species>`

NOT:
- ~~Model output ↔ Peak density~~

---

## Tools Provided

1. **`load_experimental_profile.py`**
   - Load nH(x) from CSV
   - Calculate spatial average
   - Estimate L_diff from decay
   - Create example profiles

2. **`run_with_fixed_H.py`**
   - Run simulation with fixed H
   - Compare CH and C₂ to targets
   - Interpret results

3. **This guide (USING_EXPERIMENTAL_H_PROFILE.md)**
   - Conceptual overview
   - Implementation recipes
   - Decision flowchart

---

## References

- `SPATIAL_AVERAGING_AND_H2.md` - Why spatial averaging matters
- `EXPERIMENTAL_TARGETS.md` - Target values for CG and CSB regions
- `NEXT_STEPS_CG.md` - Overall validation workflow
- `define_rates_tunable.py` - Where L_diff parameter is used
- `odefun_optimized.py` - Where to implement fixed species
