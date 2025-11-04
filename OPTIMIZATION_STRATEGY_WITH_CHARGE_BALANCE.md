# Optimization Strategy with Charge-Consistent Ne

**Date**: 2025-11-03
**Issue**: How to handle Ne when re-optimizing with Te-dependent rates

---

## The Problem You Identified

### Current Situation:
- **Ne imposed**: 1.95×10⁹ cm⁻³ (from sheath width observations)
- **Ion density produced**: 1.54×10⁸ cm⁻³ (from chemistry with Te=1 eV)
- **Charge imbalance**: 1253%

### Your Concern:
> "If we let Ne run freely, it will just equal ion densities (~1.5×10⁸), which is much lower than the observed value. How does that help?"

**Answer**: You're right! Simply letting Ne = Σn_i⁺ would give Ne ≈ 1.5×10⁸, which contradicts sheath width observations.

---

## The Physical Connection: Ne ↔ Sheath Width

### Why Ne Matters for Sheath Width:

The **cathode sheath** separates the cathode from the plasma. Its width is:

```
d_sheath ∝ (V_cathode)^(3/4) / (n_i)^(1/2)
```

From Child-Langmuir law (space-charge limited current).

**If you measured** d_sheath ≈ 0.1-0.3 mm:
- Combined with cathode fall voltage V_c ≈ 250 V
- This implies **n_i ≈ 1-2×10⁹ cm⁻³** in the CG region

**This is your constraint!** Ne ≈ 1.95×10⁹ comes from physical observations, not arbitrary choice.

---

## The Solution: Te-Dependent Rates Can Produce Higher n_i!

### Current Model (Te = 1 eV):
```
Ionization: k_ion(Ar, Te=1eV) = 8e-12 cm³/s
Production: n_i × k_ion × n_e × n_Ar = 1.5×10⁸ × 8e-12 × ...
Result: n_i = 1.5×10⁸ cm⁻³ ❌ Too low!
```

### With Te-Dependent Rates (Te = 5 eV):
```
Ionization: k_ion(Ar, Te=5eV) = 2.5e-6 cm³/s  (300,000× higher!)
Production: Much higher!
Result: n_i can reach 1-2×10⁹ cm⁻³ ✓ Matches sheath width!
```

**The key insight**: At higher Te, ionization is MUCH stronger, so the chemistry can naturally produce the high ion densities needed to match observations.

---

## Your Second Insight: E-Field Affects Te

You said:
> "Ion densities can increase if the E field is lowered"

This is interesting! Let's think through the physics:

### Effect of E-Field on Discharge:

**Lower E-field**:
- E/N decreases
- Te typically decreases (less electron heating)
- But ionization k_ion ∝ exp(-E_ion/Te) is VERY sensitive to Te
- Lower Te → much lower ionization → lower n_i ❌

**Higher E-field**:
- E/N increases
- Te increases (more electron heating)
- Ionization increases exponentially
- Higher n_i ✓

So actually, **higher E** (not lower) gives higher n_i.

### Current Model:
- **E = 300 V/cm** (very high!)
- **E/N = 2329 Td** (extreme!)
- This SHOULD produce high n_i
- But doesn't because k_ion assumes Te = 1 eV ❌

### With Te-Dependent Rates:
- **E = 300 V/cm** → **Te = 6-8 eV**
- **k_ion increases 100,000-1,000,000×**
- **n_i reaches 1-2×10⁹** ✓ Matches observations!

---

## Recommended Optimization Strategy

### **Option A: Constrained Ne (Recommended)**

**Approach**: Keep Ne near observed value, let Te adjust

```python
# Optimization parameters:
- Te: FREE (range 1-8 eV)
- Ne: CONSTRAINED near 1.95×10⁹ ± 50%
- E-field: FREE (range 100-400 V/cm)
- 40 rates: FREE (within literature bounds)

# Objective function:
error = w_H × (H - H_target)²
      + w_CH × (CH - CH_target)²
      + w_C2 × (C2 - C2_target)²
      + w_charge × (Σn_i⁺ - n_e)²  ← NEW: Charge balance penalty!
```

**Benefits**:
- ✓ Respects sheath width observations
- ✓ Enforces charge balance
- ✓ Finds Te that produces n_i ≈ Ne
- ✓ Physically self-consistent

**Expected result**:
- Te ≈ 5-7 eV (from E/N)
- n_i ≈ 1.5-2×10⁹ (matches Ne)
- Charge balance error < 10%

---

### **Option B: Free Ne with Charge Balance Constraint**

**Approach**: Let Ne be free, but enforce quasi-neutrality

```python
# Optimization parameters:
- Te: FREE (1-8 eV)
- Ne: FREE (1e8-1e10 cm⁻³)
- E-field: FREE (100-400 V/cm)
- 40 rates: FREE

# Constraint:
ENFORCE: |Σn_i⁺ - n_e| / n_e < 0.05  (5% charge balance)

# Objective function:
error = (H - H_target)² + (CH - CH_target)² + (C2 - C2_target)²
```

**Benefits**:
- ✓ Perfect charge balance (by construction)
- ✓ No assumptions about Ne
- ✓ Finds Ne from first principles

**Risk**:
- ⚠ May find Ne ≠ 1.95×10⁹
- ⚠ May not match sheath width
- Need to validate against discharge physics

---

### **Option C: Two-Stage Optimization (Most Rigorous)**

**Stage 1**: Find (Te, E) that gives n_i ≈ 1.95×10⁹

```python
# Fix Ne = 1.95×10⁹
# Vary only Te and E-field
# Minimize: |Σn_i⁺ - Ne|
# Result: Best (Te, E) for charge balance
```

**Stage 2**: Optimize rates to match targets

```python
# Fix Te and E from Stage 1
# Vary 40 reaction rates
# Minimize: (H-target)² + (CH-target)² + (C2-target)²
```

**Benefits**:
- ✓ Separates physics (Te, E) from chemistry (rates)
- ✓ Guarantees charge balance
- ✓ Respects sheath width

**Cost**:
- ⏱️ Two optimization runs
- More complex workflow

---

## Recommended: Option A with Charge Balance Penalty

### Implementation:

```python
def objective_function_with_charge_balance(x, param_names, params_base):
    # Extract parameters
    Te = x[param_names.index('Te')]
    E_field = x[param_names.index('E_field')]
    Ne = x[param_names.index('ne')]

    # Run simulation with Te-dependent rates
    params = params_base.copy()
    params['Te'] = Te
    params['E_field'] = E_field
    params['ne'] = Ne

    k = define_rates_tunable(params)  # ← Use Te-dependent version!
    # ... run simulation ...

    # Calculate ion densities
    n_i_total = sum(n_i_plus for all positive ions)

    # Species targets
    target_error = (
        1.0 * ((H - 5.18e13) / 5.18e13)**2 +
        20.0 * ((CH - 1.0e9) / 1.0e9)**2 +
        3.0 * ((C2 - 1.3e11) / 1.3e11)**2
    )

    # Charge balance penalty
    charge_error = ((n_i_total - Ne) / Ne)**2

    # Combined objective
    total_error = target_error + 10.0 * charge_error  # ← Adjust weight

    return total_error
```

### Parameter Bounds:

```python
bounds = [
    # ... 40 rate bounds from database ...
    (100, 400),      # E-field (V/cm)
    (1.0, 8.0),      # Te (eV) ← NEW!
    (5e8, 5e9),      # Ne (cm⁻³), allow ±2× variation from 1.95e9
]
```

### Expected Outcome:

The optimizer will find:
- **Te ≈ 5-7 eV** (consistent with E/N = 2329 Td)
- **E ≈ 250-350 V/cm** (adjusts to balance ionization)
- **Ne ≈ 1.5-2.5×10⁹** (near observed value)
- **n_i ≈ Ne** (charge balanced!)
- **H, CH, C2 match targets** (chemistry optimized)

---

## Why This Works

### The Feedback Loop:

1. **High E-field** → High E/N → High Te
2. **High Te** → Strong ionization (k_ion ↑ 100,000×)
3. **Strong ionization** → High n_i
4. **High n_i** → Quasi-neutrality: n_i ≈ Ne
5. **Matches sheath width observations** ✓

### Current Model Fails Because:
- Assumes Te = 1 eV (too low)
- k_ion is weak
- Can't produce enough n_i
- Charge imbalance

### New Model Succeeds Because:
- Te is free parameter
- Finds Te ≈ 6 eV from E/N
- k_ion is strong enough
- Produces n_i ≈ 2×10⁹
- Charge balanced!

---

## Your Question About E-Field

You asked:
> "Ion densities can increase if the E field is lowered"

Actually, the opposite is true for **ionization-driven** discharges:

### Increasing E:
```
E ↑ → Te ↑ → k_ion ↑ exponentially → n_i ↑
```

### Decreasing E:
```
E ↓ → Te ↓ → k_ion ↓ exponentially → n_i ↓
```

**BUT**: There's a subtlety for **diffusion-dominated** regimes:

If E is SO high that:
- Ions are driven to walls very fast (drift velocity ↑)
- Faster than they can be replaced by ionization
- Then n_i might decrease

Typical crossover: E/N ≈ 1000-2000 Td (you're at 2329 Td, right at the edge!)

So your intuition might be correct if E = 300 V/cm is **too high** and causing excessive ion loss.

**The optimization will find the right balance!**

---

## Next Steps

1. **Modify `optimize_c2_ch_targeted.py`**:
   - Use `define_rates_tunable.py`
   - Add Te as parameter (1-8 eV)
   - Add charge balance penalty

2. **Run optimization** (~30 min with 15 iterations)

3. **Compare results**:
   - What Te does it find?
   - Is n_i ≈ Ne?
   - Do targets still match?

4. **Validate**:
   - Is Te consistent with E/N?
   - Does Ne match sheath width?
   - Is charge balance good?

---

## Summary

**Your concern is valid**: Simply letting Ne = n_i would give Ne ≈ 1.5×10⁸, contradicting sheath width observations.

**The solution**: Use Te-dependent rates! At higher Te (6-8 eV), ionization is 100,000× stronger, so chemistry CAN produce n_i ≈ 2×10⁹.

**Strategy**:
- Keep Ne as optimization parameter (constrained near 1.95×10⁹)
- Add Te as parameter (1-8 eV)
- Add charge balance penalty to objective
- Let optimizer find self-consistent (Te, E, Ne)

**Expected result**: Te ≈ 6 eV, n_i ≈ 2×10⁹, charge balanced, targets matched!

Ready to implement this?
