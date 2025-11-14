# Comprehensive Rate and Literature Review
## Ar + CH4 Glow Discharge Plasma Chemistry

**Date:** 2025-11-14
**Purpose:** Audit all reaction rates and compare with literature for similar plasma conditions

---

## 1. CRITICAL RATE CONSTANT ISSUES IDENTIFIED

### 1.1 **H + C2H2 → C2 + H2 + H** (Line 315)
**Current model:**
```python
k['C2H2_H_C2_H2_H_cm3_7_50'] = 1e-11  # CONSTANT
```

**Correct (Baulch et al. 2005):**
```python
k = 1.67e-14 * T^1.64 * exp(-15250/T)  # cm³/(molecule·s)
```

**At T=570K:**
- Correct: k = 1.33×10⁻²¹ cm³/s
- Model: k = 1.00×10⁻¹¹ cm³/s
- **Error: 7.5 billion times too high!**

**Status:** ✅ IDENTIFIED AND CORRECTED

---

### 1.2 **C2 + H → CH + C** (Line 255)
**Current model:**
```python
k['C2_H_CH_C_cm3_7_6'] = 9.6e-11  # CONSTANT
```

**Issues:**
- 93% of C2 destruction
- Uses constant rate (no T-dependence)
- **Should this rate depend on C2 vibrational state?**

**Literature check needed:**
- Baulch et al. (2005) value for ground state C2(X¹Σ_g⁺)
- What is rate for vibrationally excited C2(v>0)?
- User's C2 is produced vibrationally hot

**Status:** ⚠️ NEEDS LITERATURE VERIFICATION

---

### 1.3 **CH + CH2 → C2 + H2 + H** (Line 288)
**Current model:**
```python
k['CH2_CH_C2_H2_H_cm3_7_26'] = 1.2e-10  # CONSTANT
```

**Importance:**
- **21% of C2 production** (dominant with corrected rates)
- Radical-radical recombination
- Should be weakly T-dependent or constant

**Literature check needed:**
- Verify this rate constant
- Check stoichiometry (produces H + H2 or just H2?)

**Status:** ⚠️ NEEDS VERIFICATION

---

### 1.4 **CH2 + CH2 → C2 + 2H2** (Line 314)
**Current model:**
```python
k['CH2_CH2_C2_H2_H2_cm3_7_58'] = 1e-11  # CONSTANT
```

**Importance:**
- **16% of C2 production**
- Radical-radical recombination
- Barrierless, should be ~constant

**Literature check needed:**
- Verify rate constant
- Common reaction in combustion chemistry

**Status:** ⚠️ NEEDS VERIFICATION

---

### 1.5 **C2H2 + C → C2 + CH2** (Line 268)
**Current model:**
```python
k['C2H2_C_C2_CH2_cm3_7_19'] = 1e-10  # CONSTANT
```

**Importance:**
- **18% of C2 production**
- C-atom insertion into C2H2
- May have small barrier

**Literature check needed:**
- Verify rate constant
- Check if T-dependence needed

**Status:** ⚠️ NEEDS VERIFICATION

---

### 1.6 **CH + H → C + H2** (Line 252)
**Current model:**
```python
k['CH_H_C_H2_cm3_7_3'] = 1.2e-10  # CONSTANT
```

**Importance:**
- CH destruction pathway
- Should be compared with CH + H measurement

**Status:** ⚠️ NEEDS VERIFICATION

---

## 2. TEMPERATURE-DEPENDENT REACTIONS (Correctly Implemented)

### ✅ **CH + CH4 → C2H4 + H** (Lines 270-274)
```python
k = 6.7e-11 × (T/293)^(-0.4)  # Thiesemann et al. 1997
```
- Valid range: 290-700 K
- Barrierless addition
- **Correctly implemented!**

---

### ✅ **H + CH4 → CH3 + H2** (Lines 283-287)
```python
k = 6e-12 × exp(-0.5 eV / kB·T)  # Ea = 0.5 eV
```
- Activation barrier correctly included
- **Correctly implemented!**

---

### ✅ **CH3 + H → CH2 + H2** (Lines 298-301)
```python
k = 6e-12 × exp(-0.65 eV / kB·T)  # Ea = 0.65 eV
```
- Activation barrier correctly included
- **Correctly implemented!**

---

### ✅ **H + C2H4 → C2H3 + H2** (Lines 331-334)
```python
k = 1.0e-11 × exp(-0.6 eV / kB·T)  # Ea = 0.6 eV
```
- Activation barrier correctly included
- **Correctly implemented!**

---

## 3. LITERATURE REVIEW: Ar + CH4 PLASMA CHEMISTRY

### 3.1 Similar Studies Found

**Plasma-Based CH4 Conversion (Wang et al., J. Phys. Chem. C 2020)**
- Modeling of CH4 conversion in different plasma sources
- Identifies dominant pathways for C2 hydrocarbon formation
- **Key finding:** Radicals recombine to form C2H6, C3H8, C4H10

**Chemical Reaction Studies in CH4/Ar DBD (Cernogora et al., J. Phys. Chem. A 2005)**
- CH4/Ar dielectric barrier discharge
- Optical emission spectroscopy: C, CH, C2 radicals detected
- **Key finding:** CH4 → CH3, CH2, CH radicals via electron impact

**Low Pressure Glow Discharge (N2-CH4)**
- Mass spectrometry and OES at 1-2 Torr
- Methane decomposition pathways studied
- **Relevant for low-pressure conditions**

---

### 3.2 Dominant Reaction Pathways from Literature

**CH4 Dissociation:**
1. e + CH4 → CH3 + H + e (primary)
2. e + CH4 → CH2 + H2 + e
3. e + CH4 → CH + H2 + H + e

**C2 Formation:**
1. **Radical recombination:** CH + CH2, CH2 + CH2, CH + CH
2. **C-atom reactions:** C + C2H2, C + CH
3. **Electron impact:** e + C2H2 → C2 + H2 + e (minor)

**C2 Destruction:**
1. **C2 + H → products** (dominant in H-rich plasma)
2. Wall loss
3. Diffusion

**CH Formation:**
1. **From C2:** C2 + H → CH + C (if C2 present)
2. **Direct:** e + CH4 → CH + products
3. **From CH2:** CH2 + H → CH + H2

---

## 4. COMPARISON: Model vs. Literature

### ✅ **Correctly Captured:**
- Electron impact dissociation of CH4
- Radical recombination pathways
- Temperature-dependent abstraction reactions
- Ion-molecule reactions

### ⚠️ **Potential Issues:**
1. **H + C2H2 → C2** rate was wrong (now corrected)
2. **C2 + H → CH + C** may be wrong for vibrationally excited C2
3. **Radical recombination rates** (CH+CH2, CH2+CH2) need verification
4. **Missing vibrational state-specific rates** for hot species

### ❓ **Not Clear from Literature:**
- Effect of suprathermal H atoms (17.5-110 eV) on chemistry
- C2 vibrational distribution and its effect on reactivity
- Spatial gradients (sheath vs. bulk chemistry)

---

## 5. SPECIFIC RATE CONSTANT RECOMMENDATIONS

### **Priority 1: Verify Against Baulch et al. (2005)**

Check these reactions in Baulch database:

1. ✅ **H + C2H2 → C2 + H2 + H** - CORRECTED
2. ⚠️ **C2 + H → CH + C** - Verify value and T-range
3. ⚠️ **CH + CH2 → C2 + H2 + H** - Verify rate
4. ⚠️ **CH2 + CH2 → C2 + 2H2** - Verify rate
5. ⚠️ **C2H2 + C → C2 + CH2** - Verify rate
6. ⚠️ **CH + H → C + H2** - Verify rate

### **Priority 2: Search for Vibrationally Excited Species**

Look for:
- C2(v) + H → CH + C rates
- CH(v) + H rates
- Effect of vibrational excitation on reactivity

**Literature sources to check:**
- Baulch et al. (2005) - Hydrocarbon combustion kinetics
- NIST Chemical Kinetics Database
- Tsang & Hampson (1986) - CH radical reactions
- Koshi et al. (1990s) - C2 chemistry
- Maas & Warnatz (1988) - High-temperature hydrocarbons

---

## 6. EXPERIMENTAL CONDITIONS: User's System

**Plasma Parameters:**
- P = 400 mTorr (±40/-20)
- Tgas = 570 K
- ne = 2.3×10⁹ cm⁻³ (can vary ±100%)
- Te = ~1.3-1.5 eV (can vary, especially higher)
- E-field = 50-300 V/cm
- Ar/CH4 mixture (~97/3%)

**Suprathermal H Atoms:**
- 15-25% of H atoms are fast
- ~90% of fast H are MG (17.5 eV)
- ~2-5% are BG (110 eV)
- **User notes:** Too fast to participate in chemistry

**Spatial Profiles:**
- C2 peaks at ~4 mm (broad peak)
- CH peaks at 2-4.5 mm (broad plateau)
- Both produced in bulk plasma (not at cathode)

**Target Densities (at sheath edge):**
- H = 2.3×10¹⁴ cm⁻³
- CH = 1.34×10⁹ cm⁻³
- C2 = 5.60×10¹¹ cm⁻³

---

## 7. MODEL PREDICTIONS WITH CORRECTED RATES

**With H + C2H2 correction:**

**C2 Production (Total: 5.41×10¹² cm⁻³/s):**
1. CH + CH2 → C2 + H2 + H - 21.36%
2. C2H2 + C → C2 + CH2 - 18.46%
3. CH2 + CH2 → C2 + 2H2 - 15.57%
4. H + C2H → C2 + H2 - 10.54%
5. e + C2H2 → C2 + H2 - 8.14%

**C2 Destruction (Total: 5.41×10¹² cm⁻³/s):**
1. C2 + H → CH + C - 92.59%
2. Wall sticking - 6.36%
3. Diffusion loss - 1.02%

**Result:**
- C2 = 2.75×10⁸ cm⁻³
- Target = 5.60×10¹¹ cm⁻³
- **2000× too low!**

---

## 8. CONCLUSIONS

### **Major Issues Found:**
1. ✅ **H + C2H2 → C2** rate was 7.5 billion times too high (FIXED)
2. ⚠️ **C2 + H → CH + C** may be too fast for vibrationally excited C2
3. ⚠️ **Radical recombination rates** need literature verification

### **Next Steps:**

**Immediate:**
1. Search Baulch et al. (2005) for C2 + H → CH + C rate
   - Check temperature range validity
   - Look for vibrational state dependence

2. Verify top 5 C2 production reaction rates:
   - CH + CH2 → C2 + H2 + H
   - C2H2 + C → C2 + CH2
   - CH2 + CH2 → C2 + 2H2
   - H + C2H → C2 + H2
   - e + C2H2 → C2 + H2

3. Search for C2(v) + H reaction studies
   - Theoretical calculations (RRKM, MD simulations)
   - Experimental measurements at elevated T

**Long-term:**
- Implement vibrational state-specific chemistry
- Consider two-temperature model (thermal vs. suprathermal H)
- Add spatial resolution (1D radial model)

---

## 9. KEY REFERENCES TO CONSULT

1. **Baulch et al. (2005)** - "Evaluated Kinetic Data for Combustion Modeling: Supplement II"
   - Primary source for C2, CH, C chemistry
   - Temperature ranges: 300-3000 K

2. **NIST Chemical Kinetics Database** - https://kinetics.nist.gov/
   - Comprehensive rate constant compilation
   - Search by reaction formula

3. **Tsang & Hampson (1986)** - "Chemical Kinetic Data Base for Combustion Chemistry. Part I. Methane and Related Compounds"
   - CH radical reactions
   - Reliable T-dependent rates

4. **Koshi et al. (1992)** - "Reactions of C2 with H2, H, and CH3"
   - Specific C2 chemistry
   - Shock tube measurements

5. **Maas & Warnatz (1988)** - "Ignition processes in hydrogen-oxygen mixtures"
   - C1-C2 hydrocarbon chemistry
   - Validated kinetic models

6. **Wang et al. (2020)** - "Plasma-Based CH4 Conversion into Higher Hydrocarbons and H2"
   - Recent plasma chemistry modeling
   - Comparison of different plasma sources

---

**END OF REVIEW**
