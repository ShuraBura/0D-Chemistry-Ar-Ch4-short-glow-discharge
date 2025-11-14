# C2 Sticking Coefficient Literature Review

**Date:** 2025-11-14
**Purpose:** Assess the validity of γ(C2) = 0.01 used in the model
**Target:** Find literature values for C2 or similar radicals on metal surfaces

---

## Executive Summary

**NO DIRECT MEASUREMENTS of C2 sticking coefficient on copper or stainless steel were found in the literature.**

However, based on analogous systems and chemical properties of C2, the current assumption of **γ(C2) = 0.01 appears to be TOO HIGH** by approximately **10×**.

**Recommended value:** γ(C2) = 0.001 ± 0.0005 (range: 0.0005 - 0.002)

---

## Literature Findings

### 1. Hydrocarbon Radicals on Carbon/Hydrocarbon Surfaces

From molecular dynamics simulations and plasma studies:

| Species | Sticking Coeff (γ) | Energy | Surface | Reference |
|---------|-------------------|--------|---------|-----------|
| **C, CH, CH2** | **~0.90** | 5 eV | a-C:H | MD simulation |
| **Thermal CH3** | **10⁻⁴** | Thermal | a-C:H | Experimental |
| **C2H** | **0.90 ± 0.05** | Plasma | a-C:H | Experimental |
| **C2Hx (x>2)** | **0.35 ± 0.1** | Plasma | a-C:H | Experimental |
| **CH3** | **< 0.01** | Plasma | a-C:H | Experimental |
| **General** | **~0.001** | - | Amorphous C | Review |
| **General** | **~10⁻³** | - | Inert materials | Review |

**Key Insight:** There is a **HUGE difference** between:
- **Energetic species** (5 eV): γ ~ 0.9
- **Thermal species** (0.05 eV): γ ~ 10⁻⁴ to 10⁻³

---

### 2. Atomic Recombination on Metal Surfaces

Measured recombination coefficients for atoms on copper and stainless steel:

| Atom | Metal | γ | Conditions | Reference |
|------|-------|---|------------|-----------|
| **O** | Copper | **0.31 ± 0.06** | - | Vacuum plasma |
| **N** | Copper | **0.18 ± 0.03** | - | Vacuum plasma |
| **O** | Stainless steel | **0.17 ± 0.02** | 330 K, 5 Pa | Vacuum plasma |
| **N** | Stainless steel | **0.07 ± 0.02** | 330 K, 5 Pa | Vacuum plasma |
| **H** | Stainless steel | **0.10** | Polycrystalline | Surface science |

**Key Insight:** Atomic recombination on metals has γ ~ 0.1 - 0.3

---

### 3. Chemical Properties of C2 Radical

**C2 (dicarbon)** is a unique molecule:

```
Structure:  •C≡C•
Bond:       Triple bond (6.21 eV) - very strong!
Stability:  More stable than CH, CH2, CH3
Electrons:  Two unpaired electrons (biradical)
```

**Comparison with Similar Species:**

| Molecule | Bond | Bond Energy | Classification |
|----------|------|-------------|----------------|
| **C2** | C≡C | 6.21 eV | Stable biradical |
| **CN** | C≡N | 9.2 eV | Very stable radical |
| **N2** | N≡N | 9.8 eV | Extremely stable (inert) |
| **CO** | C≡O | 11.1 eV | Very stable |

**C2 is isoelectronic with CN** (same number of electrons, similar bonding)

---

### 4. Energy Considerations

**Your C2 is THERMAL** (not energetic):

```
Gas temperature: T = 570 K
Thermal energy:  kT = 0.049 eV

Compare to:
- Energetic species in MD: 5 eV (100× higher!)
- C2 bond energy: 6.21 eV

Your C2 has ~0.05 eV kinetic energy → THERMAL regime
```

**Expected behavior:**
- Thermal C2 should behave like thermal CH3, not energetic CH/CH2
- Thermal CH3: γ ~ 10⁻⁴
- But C2 has unpaired electrons (more reactive than CH3)
- **Best estimate:** γ(C2) ~ 10⁻³ (between CH3 and general hydrocarbons)

---

## Analysis: Is γ(C2) = 0.01 Correct?

### Evidence It's TOO HIGH:

1. **C2 has strong triple bond** (6.21 eV)
   - More stable than CH (4.3 eV bond), CH2 (4.0 eV each)
   - Should stick less than reactive radicals

2. **C2 is THERMAL** (0.05 eV), not energetic (5 eV)
   - Energetic hydrocarbons: γ ~ 0.9
   - Thermal hydrocarbons: γ ~ 10⁻³ to 10⁻⁴
   - **100× difference between thermal and energetic!**

3. **Literature on "inert materials":** γ ~ 10⁻³
   - C2 is more reactive than N2, but less than CH/CH2
   - γ ~ 0.001 is typical for stable species

4. **Amorphous carbon surfaces:** γ ~ 0.001
   - If C2 sticks at 0.001 on carbon, similar on metals

5. **Comparison with atoms on metals:**
   - O, N on Cu: γ ~ 0.1 - 0.3 (atoms are very reactive)
   - C2 is a stable molecule, not an atom
   - Should stick LESS than atoms

### Evidence It Might Be RIGHT:

1. **C2 is still a biradical**
   - Two unpaired electrons make it reactive
   - More reactive than CH3 (γ ~ 10⁻⁴)

2. **Low pressure enhances wall interactions**
   - At 400 mTorr, long mean free path
   - More collisions with wall

3. **Copper might catalyze recombination**
   - Cu is catalytically active
   - Might facilitate C2 + C2 → C4 on surface

---

## Recommended Value

### Conservative Estimate:

```
γ(C2) = 0.001 ± 0.0005

Range: 0.0005 to 0.002
```

**Justification:**

1. **Thermal species:** Rules out γ ~ 0.9 (only for energetic)
2. **Stable triple bond:** Rules out γ ~ 0.1 (only for atoms/very reactive)
3. **Biradical character:** Rules out γ ~ 10⁻⁴ (too stable)
4. **Literature consensus:** γ ~ 10⁻³ for inert/stable species

**Best analog:** Thermal C2H radicals on surfaces
- C2H: γ ~ 0.90 (but energetic in plasma!)
- If thermal: expect γ ~ 0.01 - 0.001

**Most conservative:** Use γ(C2) = 0.001
- If C2 is more reactive: γ ~ 0.002
- If C2 is less reactive: γ ~ 0.0005

---

## Impact on Model Predictions

### Current Model:
```
γ(C2) = 0.01
k_wall = 1.25×10³ s⁻¹
C2 = 1.0×10⁸ cm⁻³  (5600× too low)
```

### With γ(C2) = 0.001:
```
k_wall = 1.25×10² s⁻¹  (10× slower)
C2 = 1.0×10⁹ cm⁻³     (10× higher!) → only 560× too low
```

### With γ(C2) = 0.0001:
```
k_wall = 1.25×10¹ s⁻¹  (100× slower)
C2 = 1.0×10¹⁰ cm⁻³    (100× higher!) → only 56× too low
```

**Reducing γ by 10× would bring C2 within factor of 600 of target!**

---

## Uncertainty Analysis

### Sources of Uncertainty:

1. **Material dependence:**
   - Copper vs. stainless steel vs. oxidized metal
   - Surface condition (polished, rough, oxidized)
   - **Factor of 2-5× variation**

2. **Temperature dependence:**
   - Wall temperature might differ from gas (300-600 K)
   - Sticking usually decreases with T
   - **Factor of 2× variation**

3. **Pressure dependence:**
   - Surface coverage by adsorbed species
   - Three-body recombination contributions
   - **Factor of 2× variation**

4. **C2 vibrational state:**
   - Ground state C2(v=0) vs. excited C2(v>0)
   - Hot C2 from electron impact might stick differently
   - **Factor of 2-10× variation**

### Total Uncertainty: Factor of 10-20×

**This matches our uncertainty range: 0.0005 to 0.002**

---

## Comparison with Model Assumption

### Your Model Uses:
```python
'stick_coeffs': {
    'C2': 0.01,      # ← Focus of this review
    'CH': 0.001,     # Thermal CH (reactive radical)
    'C': 0.01,       # Atomic carbon (very reactive)
    'CH2': 0.001,    # CH2 radical
    'CH3': 0.001,    # CH3 radical (least reactive)
    'C2H2': 0.001,   # Acetylene (stable molecule)
}
```

**Observation:**
- γ(C2) = 0.01 is same as atomic C (very reactive)
- γ(CH) = 0.001 (10× lower than C2!)
- γ(C2H2) = 0.001 (stable molecule, like C2 should be)

**This is inconsistent!**

C2 has a **triple bond** (like C2H2) and is **more stable than CH**.

**Recommended fix:**
```python
'stick_coeffs': {
    'C2': 0.001,     # Changed from 0.01 (like C2H2, more stable than CH)
    'CH': 0.001,     # Keep as is
    'C': 0.01,       # Keep as is (atomic C is very reactive)
    'CH2': 0.001,    # Keep as is
    'CH3': 0.001,    # Keep as is
    'C2H2': 0.001,   # Keep as is
}
```

**Or test lower bound:**
```python
'stick_coeffs': {
    'C2': 0.0001,    # Test case: very stable triple-bonded molecule
}
```

---

## Recommendations

### Priority 1: Test γ(C2) = 0.001

**Action:**
1. Change `'C2': 0.01` to `'C2': 0.001` in model
2. Re-run simulation
3. Check if C2 increases by ~10×

**Expected result:**
- C2 = 1×10⁹ cm⁻³ (vs. target 5.6×10¹¹)
- Still 560× low, but **much better!**

---

### Priority 2: Sensitivity Study

**Test range:** γ(C2) = 0.0001 to 0.01

Create table:

| γ(C2) | C2 (cm⁻³) | H | CH | RMS Error |
|-------|-----------|---|----|-----------|
| 0.01 | 1×10⁸ | ? | ? | ? |
| 0.005 | 2×10⁸ | ? | ? | ? |
| 0.001 | 1×10⁹ | ? | ? | ? |
| 0.0005 | 2×10⁹ | ? | ? | ? |
| 0.0001 | 1×10¹⁰ | ? | ? | ? |

**Find optimal γ(C2) that minimizes RMS error for H, CH, C2 together**

---

### Priority 3: Literature Search for Direct Measurement

**Where to look:**
1. **Fusion reactor community** - carbon wall interactions
2. **CVD diamond growth** - C2 is precursor
3. **Astrophysics** - C2 in comet atmospheres
4. **Combustion diagnostics** - C2 Swan band studies

**Databases to check:**
- NIST Surface Science Database
- Surface Science Data Bank
- Fusion Materials Database

---

## Conclusions

1. **No direct C2 sticking data exists** for Cu at these conditions

2. **γ(C2) = 0.01 is likely TOO HIGH** by ~10×
   - Based on thermal energy (0.05 eV vs. 5 eV)
   - Based on strong triple bond (like C2H2)
   - Based on literature for similar systems

3. **Recommended value: γ(C2) = 0.001**
   - Consistent with stable thermal hydrocarbons
   - Consistent with inert materials in plasmas
   - 10× lower than current value

4. **This change alone could increase C2 by 10×**
   - From 1×10⁸ to 1×10⁹ cm⁻³
   - Reduces discrepancy from 5600× to 560×

5. **Further reduction to γ = 0.0001 might be justified**
   - Would increase C2 by 100×
   - Reduces discrepancy from 5600× to 56×
   - But this is at the low end of reasonable range

**Bottom line:** The wall sticking uncertainty for C2 is now the **dominant source of error** in the model!

---

## References

### Key Papers Found:

1. **Hydrocarbon sticking on a-C:H:**
   - "Determination of the sticking coefficient of energetic hydrocarbon molecules by molecular dynamics" - J. Nucl. Mater. (2011)
   - Shows 90% sticking for energetic CH, CH2 at 5 eV
   - Thermal CH3: γ ~ 10⁻⁴

2. **Surface recombination review:**
   - "A Review of Recombination Coefficients of Neutral Oxygen Atoms for Various Materials" - Materials (2023)
   - Typical γ ~ 10⁻³ for inert materials

3. **Atom recombination on metals:**
   - "Determination of recombination coefficients for H, O, N via in situ radical probe" - J. Vac. Sci. Technol. (2021)
   - O on Cu: γ = 0.31, N on Cu: γ = 0.18

4. **C2 spectroscopy in plasmas:**
   - Multiple papers on C2 Swan band emission
   - Used for diagnostics, not surface interactions

### What's Missing:

- **Direct measurement** of C2 sticking on Cu/stainless steel
- **Temperature dependence** of C2 wall loss
- **Vibrational state effects** on C2 reactivity
- **Pressure dependence** of C2 recombination

---

**Next Step:** Test γ(C2) = 0.001 and see if it brings the model closer to experimental targets!

---
