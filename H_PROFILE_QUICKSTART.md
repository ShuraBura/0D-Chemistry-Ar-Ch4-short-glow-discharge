# Quick Start: Using Your nH(x) Profile

## TL;DR

You have measured nH(x) over 0-1 mm. Here's what to do:

### 1️⃣ Calculate the spatial average:
```python
import numpy as np
x = np.array([...])     # Your positions (cm)
nH = np.array([...])    # Your H densities (cm⁻³)
nH_avg = np.trapz(nH, x) / (x[-1] - x[0])
```

### 2️⃣ Choose your approach:

**Option A: Let H evolve (recommended for validation)**
```python
# In your simulation script
set_density('H', nH_avg)  # Initial condition
# Run simulation normally
# Check if final H ≈ nH_avg
```

**Option B: Fix H density (recommended for focusing on CH/C2)**
```python
# In odefun.py, add after line 104:
H_idx = self.species.index('H')
self.dydt[H_idx] = 0  # Keep H constant

# Or just run:
python3 run_with_fixed_H.py
```

### 3️⃣ Compare results:
- **If Option A works**: Your model is self-consistent! ✓
- **If Option A fails**: Use Option B to validate other species

---

## Decision Tree

```
Do you want to validate H production/loss chemistry?
│
├─ YES → Use Option A (dynamic H)
│         - Test if model maintains experimental H
│         - Full chemistry validation
│         - May need parameter tuning
│
└─ NO  → Use Option B (fixed H)
          - Focus on CH and C₂ validation
          - H matches by design
          - Faster, simpler
```

---

## Example Workflow

### Step 1: Prepare your data

Save your profile as CSV:
```
x_cm,nH_cm3
0.000,5.0e15
0.001,6.2e15
...
```

### Step 2: Calculate average

```bash
python3 load_experimental_profile.py
```

Output:
```
Spatial average: 8.58e15 cm⁻³
Use this value in your model!
```

### Step 3A: Dynamic H (validation mode)

```python
# In run_cg_optimized.py or parameter_sweep_cg.py
set_density('H', 8.58e15)  # From experiment

# Run simulation
result = run_simulation(params)

# Check final H
final_H = result.y[H_idx, -1]
if abs(final_H - 8.58e15) / 8.58e15 < 0.5:
    print("✓ Model validates!")
```

### Step 3B: Fixed H (pragmatic mode)

```bash
python3 run_with_fixed_H.py
```

Output will show:
- CH: Target vs Predicted
- C₂: Target vs Predicted
- H: Matches exactly (by design)

---

## Why Spatial Average?

Your TALIF measures:
```
<nH> = (1/L) ∫ nH(x) dx  (average over 0-1 mm)
```

Your 0-D model predicts:
```
nH = constant  (spatially uniform)
```

**Comparison**: Model's `nH` ↔ Experiment's `<nH>`

---

## Current Experimental Value

From EXPERIMENTAL_TARGETS.md:

**CG region (0-1 mm)**:
- H = **8.58×10¹⁵ cm⁻³** (spatially averaged)

This is your target!

---

## Files Created for You

1. **`load_experimental_profile.py`** - Load and process your nH(x) data
2. **`run_with_fixed_H.py`** - Run simulation with fixed H
3. **`USING_EXPERIMENTAL_H_PROFILE.md`** - Detailed guide
4. **`H_PROFILE_QUICKSTART.md`** - This quick reference

---

## Still Confused?

**Q: My model gives H = 1e14, but experiment is 8.58e15. What now?**

A: Try these in order:
1. Run parameter sweep: `python3 parameter_sweep_cg.py`
2. Check if H₂ dissociation is included (should be!)
3. Reduce L_diff (slower wall loss)
4. Increase Te (more dissociation)
5. Use fixed H (Option B) to focus on other species

**Q: Should I match the peak or average?**

A: **Average!** The 0-D model represents the spatially-averaged region.

**Q: Can I use both peak and average?**

A: Use average for CG region (0-1 mm). If you have measurements at different locations (CSB, bulk), treat them as separate regions with different parameters.

---

## Next Steps

1. ✓ Calculate nH_avg from your profile
2. ✓ Choose Option A or B based on your goal
3. ✓ Run simulation
4. ✓ Compare to targets (H, CH, C₂)
5. If needed: Parameter sweep or add missing physics

Good luck! 🚀
