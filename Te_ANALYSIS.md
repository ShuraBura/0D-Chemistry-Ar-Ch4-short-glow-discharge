# Electron Temperature (Te) Analysis for Best Model

**Date**: 2025-11-03
**Pressure**: 0.4 Torr (user confirmed)

---

## Critical Finding: Te Mismatch

### **Assumed Te in Simulation**
- **define_rates.py** header states: **"Te ~ 1 eV"**
- All electron-impact rate coefficients are **evaluated at Te = 1 eV**
- This is a **fixed assumption**, not calculated dynamically

### **Expected Te from Discharge Conditions**

Using the best model parameters:
- **Pressure**: 0.4 Torr = 53.3 Pa
- **Gas density**: 1.29×10¹⁶ cm⁻³
- **Electric field**: 300 V/cm
- **Reduced field**: **E/N = 2329 Td**

| E/N Range (Td) | Expected Te (eV) |
|----------------|------------------|
| < 100 | 1-2 |
| 100-500 | 2-4 |
| 500-1000 | 4-6 |
| **> 1000** | **6-10+** |

**For E/N = 2329 Td → Expected Te ≈ 6-10 eV**

---

## The Problem

**MISMATCH**: The model assumes **Te = 1 eV** but the discharge conditions suggest **Te = 6-10 eV**

This is a factor of **6-10× difference**!

---

## Why This Matters

### Electron-Impact Rates are Highly Te-Dependent

Rate coefficients typically scale as:
```
k(Te) ∝ exp(-E_threshold / Te)
```

For example:
- **Ionization** (e + Ar → Ar⁺ + 2e): Threshold = 15.76 eV
  - At Te = 1 eV: k ~ 1e-12 cm³/s
  - At Te = 5 eV: k ~ 1e-9 cm³/s (1000× higher!)

- **Dissociation** (e + CH₄ → CH₃ + H): Threshold = 9 eV
  - At Te = 1 eV: k ~ 1e-12 cm³/s
  - At Te = 5 eV: k ~ 5e-11 cm³/s (50× higher)

### Impact on Results

If actual Te = 6-10 eV instead of 1 eV:
- **Ionization rates** would be 100-1000× higher
- **Dissociation rates** would be 10-100× higher
- **Excitation rates** would increase significantly
- **Ion/radical densities** would be much higher
- **Neutral depletion** would be faster

---

## What Are the Limits on Te?

### 1. **Lower Limit: ~0.3-0.5 eV**
- Below this, electrons don't have enough energy for any inelastic collisions
- No ionization, dissociation, or excitation
- Discharge cannot be sustained
- **Minimum viable Te ≈ 0.5 eV** for this chemistry

### 2. **Upper Limit: ~10-15 eV**
- At P = 0.4 Torr, E = 300 V/cm → Te can reach 10+ eV
- Beyond ~15 eV, energy loss to radiation becomes dominant
- Runaway ionization may occur
- **Practical upper limit ≈ 10-12 eV** for low-pressure DC discharges

### 3. **Typical Range for CSB/CG Region**
- **CG (Cathode Glow)**: Te = 2-5 eV typically
- **CSB (Cathode Sheath Boundary)**: Te = 3-8 eV
- **Positive Column**: Te = 1-3 eV

For E/N = 2329 Td at 0.4 Torr:
- **Realistic Te range: 4-8 eV**

---

## Current Simulation Constraints

The simulation **does NOT solve for Te dynamically**. Instead:

1. **Fixed Rate Coefficients**: All k values in `define_rates.py` are constants
2. **Assumed Te = 1 eV**: Rate coefficients evaluated at Te = 1 eV
3. **No Temperature Dependence**: k values do not change with E-field or density

### Where Te is "Baked In"

From `define_rates.py` (lines 1-4):
```python
"""
define_rates.py for Ar/CH4 plasma with Te ~ 1 eV
Goals: Adjust electron-impact and recombination rates for Te = 1 eV
Sources: Morgan (1992), Janev & Reiter (2002), ...
"""
```

**This means**:
- Electron-impact ionization rates assume Te = 1 eV
- Electron-impact dissociation rates assume Te = 1 eV
- Excitation rates assume Te = 1 eV
- Recombination rates may have weak Te dependence

---

## How This Affects Your Results

### **Good News**:
- The optimization **compensated** for the Te mismatch
- By tuning 40 rate coefficients, it found parameters that match experiments
- The model is **empirically validated** against your targets

### **Bad News**:
- The rate coefficients are **not physically self-consistent**
- If you change E-field or pressure, predictions will be inaccurate
- The model cannot predict Te or how it varies spatially

### **Impact on Charge Balance**:
- At Te = 1 eV, ionization is weak → low ion density (~1.5×10⁸ cm⁻³)
- At Te = 6-8 eV, ionization is strong → high ion density (~2×10⁹ cm⁻³)
- **This explains the charge imbalance!**
  - Simulation: n_i⁺ = 1.5×10⁸ (consistent with Te = 1 eV)
  - Imposed: n_e = 1.95×10⁹ (consistent with Te = 6-8 eV)

---

## Recommendations

### **Option 1: Fix Te at Realistic Value (Quick)**
- Recalculate all electron-impact rates for Te = 4-6 eV
- Use BOLSIG+ or look up rates from Phelps, Janev-Reiter databases
- Re-optimize with corrected rates

**Pros**: More physical, better predictions
**Cons**: Requires rate coefficient updates

### **Option 2: Self-Consistent Te (Medium)**
- Add electron energy equation:
  ```
  d(n_e × ε_e)/dt = P_electric - P_loss_elastic - P_loss_inelastic
  ```
- Solve for Te dynamically
- Update rate coefficients as Te(t) evolves

**Pros**: Fully self-consistent
**Cons**: More complex, slower simulations

### **Option 3: Use as Empirical Model (Accept Current)**
- Keep Te = 1 eV assumption
- Accept that rate coefficients are "effective" values
- Use model only for conditions similar to your experiments

**Pros**: Model already validated
**Cons**: Limited predictive power

### **Option 4: Two-Temperature Model (Advanced)**
- Solve for both Te and T_heavy separately
- Use local field approximation: Te = f(E/N)
- More accurate for high E/N discharges

**Pros**: Best physics
**Cons**: Most complex

---

## Summary Table

| Parameter | Current Model | Physical Expectation | Mismatch |
|-----------|--------------|---------------------|----------|
| **Te (assumed)** | 1 eV | 4-8 eV | 4-8× |
| **E/N** | 2329 Td | 2329 Td | ✓ |
| **k_ionization** | ~1e-11 cm³/s | ~1e-9 cm³/s | 100× |
| **n_e (imposed)** | 1.95e9 | - | - |
| **n_i⁺ (calculated)** | 1.5e8 | 1.95e9 | 13× |
| **Charge balance** | Poor (1253%) | Good (<5%) | ✗ |

---

## Physical Limits on Te

### Fundamental Constraints:
1. **Lower bound**: Te > threshold energies / 3
   - For Ar ionization (15.76 eV): Te > 5 eV for significant ionization
   - For CH₄ dissociation (9 eV): Te > 3 eV for significant dissociation
   - **Minimum Te ≈ 0.5 eV** (thermal electrons)

2. **Upper bound**: Energy balance limits
   - Power in: j·E (electron current × field)
   - Power out: elastic + inelastic collisions + radiation
   - **Maximum Te ≈ 10-12 eV** at P = 0.4 Torr

3. **Practical range for CSB region**: **Te = 2-6 eV**

---

## Action Items

1. **Immediate**: Use model as-is for qualitative trends
2. **Short term**: Update rates for Te = 4-6 eV and re-run
3. **Long term**: Implement self-consistent Te solver
4. **For publication**: Acknowledge Te assumption in methods section

---

## References

- **Phelps Database**: https://www.lxcat.net (electron-impact cross sections)
- **Janev & Reiter (2002)**: "Collision processes of CHy and CHy+ hydrocarbons"
- **Hagelaar & Pitchford (2005)**: "Solving the Boltzmann equation to obtain electron transport coefficients and rate coefficients for fluid models"
- **BOLSIG+**: Software for solving Boltzmann equation at various E/N

---

**Bottom Line**: Your model assumes Te = 1 eV, but physics suggests Te = 4-8 eV. The optimization compensated for this mismatch, so results match experiments. However, for predictive modeling, you should update rate coefficients to reflect actual Te.
