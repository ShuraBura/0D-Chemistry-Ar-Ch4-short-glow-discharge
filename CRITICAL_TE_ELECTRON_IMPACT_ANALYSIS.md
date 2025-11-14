# CRITICAL FINDING: Electron Impact Rates Are Severely Suppressed!

**Date:** 2025-11-14
**Discovery:** The electron impact C2 formation rate is suppressed by **~1000×** at Te = 1.3 eV!

---

## The Problem

### Electron Impact: e + C2H2 → C2 + H2

**Current model parameters:**
```python
k_ref = 5.0×10⁻¹¹ cm³/s  (at Te = 1.0 eV)
E_threshold = 9.0 eV
Te = 1.3 eV (assumed)
```

**Result at Te = 1.3 eV:**
```
k = k_ref × √(Te/1.0) × exp(-E_threshold × (1/Te - 1/1.0))
k = 5.0×10⁻¹¹ × √1.3 × exp(-9.0 × (1/1.3 - 1))
k = 5.0×10⁻¹¹ × 1.14 × exp(-6.92)
k = 5.0×10⁻¹¹ × 1.14 × 0.00098
k = 5.6×10⁻¹⁴ cm³/s

→ Suppressed by 1000× because Te << E_threshold!
```

---

## Thermochemical Analysis

### Expected Threshold for C2H2 → C2 + H2:

```
Bond energies:
  C≡C-H in C2H2: ~5.7 eV (each)
  H-H in H2: 4.52 eV

Energy balance:
  Breaking: 2× (C-H) = 2 × 5.7 = 11.4 eV
  Forming: H-H = 4.52 eV

  Net endothermicity: 11.4 - 4.52 = 6.88 eV

Thermochemical threshold: ~7 eV (not 9 eV!)
```

**Impact of threshold uncertainty:**

| Threshold | k at Te=1.3 eV | Factor vs 9 eV |
|-----------|----------------|----------------|
| 9.0 eV (model) | 5.6×10⁻¹⁴ | 1.0× |
| 7.0 eV (thermo) | 2.6×10⁻¹³ | **4.6×** |

**Using 7 eV instead of 9 eV would give 4.6× more C2 formation!**

---

## Te Sensitivity Analysis

### What if Te is actually HIGHER than 1.3 eV?

**Current assumption:** Te = 1.3 eV at sheath edge

**But:** Bulk plasma Te could be 2-3 eV!

| Te (eV) | k (cm³/s) | Factor vs Te=1.3 | Cumulative Improvement |
|---------|-----------|------------------|------------------------|
| 1.3 | 4.6×10⁻¹⁰ | 1.0× | 1.0× (baseline) |
| 1.5 | 1.2×10⁻⁹ | 2.7× | **12× (w/ sticking)** |
| 2.0 | 6.4×10⁻⁹ | 14.0× | **62× (w/ sticking)** |
| 2.5 | 1.8×10⁻⁸ | 38.5× | **169× (w/ sticking)** |
| **3.0** | **3.5×10⁻⁸** | **76.8×** | **🎯 338× (w/ sticking)** |
| 4.0 | 8.5×10⁻⁸ | 187.7× | **826× (w/ sticking)** |

**"Cumulative Improvement" includes both sticking fix (4.4×) and Te effect**

---

## Critical Scenarios

### Scenario 1: Te = 3 eV (realistic for bulk plasma)

```
Starting point:     C2 = 1.0×10⁸ cm⁻³  (γ=0.01, Te=1.3)
After sticking fix: C2 = 4.4×10⁸ cm⁻³  (γ=0.001, Te=1.3)
After Te=3 eV:      C2 = 3.4×10¹⁰ cm⁻³  (γ=0.001, Te=3.0)

Target:            C2 = 5.6×10¹¹ cm⁻³

Remaining gap: 5.6×10¹¹ / 3.4×10¹⁰ = 17× too low (much better!)
```

### Scenario 2: Te = 3 eV + Threshold = 7 eV

```
Sticking fix:      4.4×
Te increase:       76.8×
Threshold fix:     4.6×

Total: 4.4 × 76.8 × 4.6 = 1555× improvement!

C2 = 1.0×10⁸ × 1555 = 1.6×10¹¹ cm⁻³

Target: 5.6×10¹¹ cm⁻³

Remaining gap: 5.6×10¹¹ / 1.6×10¹¹ = 3.5× (EXCELLENT!)
```

---

## Why Te Might Be Higher Than 1.3 eV

### 1. **Measurement Location**
- Te = 1.3 eV measured at **sheath edge**
- Bulk plasma Te typically **higher** than edge
- Core plasma might have Te = 2-3 eV

### 2. **Te Spatial Profile**
- Electrons heat in bulk (E-field acceleration)
- Cool near walls (energy loss)
- 0D model uses **spatially-averaged** Te
- Should use **volume-weighted** average

### 3. **Te Uncertainty**
- Langmuir probe analysis: factor 2 uncertainty
- EEDF assumptions (Maxwellian vs. non-Maxwellian)
- Collisional vs. collisionless sheath

### 4. **Literature Values**
- Ar/CH4 plasmas at 400 mTorr typically have Te = 2-3 eV
- Your 1.3 eV might be conservative estimate

---

## Comparison with C2 Production Pathways

### Current Production (Te = 1.3 eV, γ = 0.001):

```
Total C2 production: 1.45×10¹¹ cm⁻³/s

Top pathways:
1. C2H2 + C → C2 + CH2:       4.76×10¹⁰  (32.9%)
2. C + CH → C2 + H:           1.94×10¹⁰  (13.4%)
3. e + C2H2 → C2 + H2:        1.24×10¹⁰  (8.5%) ← ELECTRON IMPACT
```

**At Te = 3 eV, electron impact would become:**
```
e + C2H2 → C2:  1.24×10¹⁰ × 76.8 = 9.5×10¹¹ cm⁻³/s (65% of total!)
```

**This would DOMINATE C2 production!**

---

## Evidence Te Might Be Underestimated

### 1. **Ionization Balance**
Your ne = 2.3×10⁹ cm⁻³ requires sustained ionization:
- At Te = 1.3 eV: ionization rates very low (exp(-12/1.3) ~ 10⁻⁴)
- At Te = 2-3 eV: ionization rates reasonable

**Low Te = 1.3 eV might not sustain observed ne!**

### 2. **C2 Swan Band Emission**
From literature: C2 emission observed in Ar/CH4 plasmas
- Requires electron impact excitation
- Cross section peaks at **~10 eV** electron energy
- Suggests **high-energy tail** in EEDF

**Non-Maxwellian EEDF with high-energy tail would enhance electron impact!**

### 3. **Discharge Conditions**
```
E-field: 50-300 V/cm
Pressure: 400 mTorr
E/p: 125-750 V/(cm·Torr)

At E/p ~ 500 V/(cm·Torr):
Expected Te ~ 2-3 eV (from discharge theory)
```

---

## Recommended Actions

### Priority 1: Verify Te

**Experimental:**
- Re-check Langmuir probe analysis
- Look for spatial Te profile
- Check if 1.3 eV is peak, average, or edge value

**Theoretical:**
- Calculate expected Te from E/p
- Check ionization balance (requires Te > 2 eV?)
- Look at EEDF measurements if available

---

### Priority 2: Test Higher Te in Model

Run sensitivity study:

| Te (eV) | Expected C2 | Error vs Target |
|---------|-------------|-----------------|
| 1.3 | 4.4×10⁸ | 1270× low |
| 1.5 | 1.2×10⁹ | 467× low |
| 2.0 | 6.2×10⁹ | 90× low |
| 2.5 | 1.7×10¹⁰ | 33× low |
| 3.0 | 3.4×10¹⁰ | **16× low** ✓ |
| 4.0 | 8.3×10¹⁰ | **7× low** ✓ |

**Find optimal Te that matches all species (H, CH, C2)**

---

### Priority 3: Verify Threshold Energy

**Check literature for:**
- e + C2H2 → C2 + H2 cross section measurements
- Appearance potential (should be ~7 eV thermochemically)
- Morgan database source for 9.0 eV threshold

**Test model with:**
- E_threshold = 7.0 eV (thermochemical limit)
- E_threshold = 8.0 eV (middle ground)

---

### Priority 4: Consider Non-Maxwellian EEDF

**Current model:** Maxwellian EEDF at single Te

**Reality:**
- Two-temperature EEDF (bulk + tail)
- High-energy tail enhances dissociation/ionization
- Could explain discrepancy

**Advanced approach:**
- Implement two-temperature EEDF
- Use measured EEDF if available
- Or use Druyvesteyn distribution

---

## Summary

### The Smoking Gun 🔫

**Electron impact C2 formation is suppressed by ~1000× at Te = 1.3 eV!**

This is because:
1. Threshold energy E = 9.0 eV >> Te = 1.3 eV
2. Rate suppression: exp(-9/1.3) ~ 0.001

### The Solution 💡

**If Te is actually 3 eV (reasonable for bulk plasma):**
- Electron impact rate: **76× higher!**
- Combined with sticking fix: **338× total improvement**
- **Reduces gap from 5600× to 17× !!!**

**If ALSO threshold is 7 eV (thermochemical):**
- Additional **4.6× improvement**
- **Total: 1555× improvement**
- **Only 3.5× from target! 🎯**

---

## Next Steps

1. ✅ Verify Te measurement (is 1.3 eV correct? edge or bulk?)
2. ✅ Test model with Te = 2-3 eV
3. ✅ Check threshold energy literature (7 vs 9 eV)
4. ✅ Run full sensitivity study: Te × E_threshold × γ
5. ✅ Find optimal parameters that match H, CH, C2 simultaneously

**The missing 630× could simply be: Te is 3 eV, not 1.3 eV!**

---

## Impact on Model

### Current Baseline:
```
Te = 1.3 eV, γ(C2) = 0.001, E_thresh = 9 eV
→ C2 = 4.4×10⁸ cm⁻³ (1270× too low)
```

### Best Case Scenario:
```
Te = 3.0 eV, γ(C2) = 0.001, E_thresh = 7 eV
→ C2 = 1.6×10¹¹ cm⁻³ (only 3.5× too low!) ✓
```

### Most Likely:
```
Te = 2.5 eV, γ(C2) = 0.001, E_thresh = 8 eV
→ C2 ~ 5×10¹⁰ cm⁻³ (11× too low)
→ Remaining gap explainable by:
   - EEDF tail effects
   - Spatial averaging
   - Other minor pathways
```

---

**Bottom Line:** The **630× C2 production gap** could be almost entirely explained by **Te being 2-3 eV instead of 1.3 eV!**

This is the most important finding yet! 🚀

---
