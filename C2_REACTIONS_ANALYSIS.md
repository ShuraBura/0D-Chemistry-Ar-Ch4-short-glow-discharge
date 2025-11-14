# Complete C2 Reaction Analysis

**Date:** 2025-11-14
**Purpose:** Catalog ALL C2 reactions in model and identify which are valid

---

## C2 FORMATION PATHWAYS (Products contain C2)

| Line | Reaction | Rate Constant | Status | Notes |
|------|----------|---------------|--------|-------|
| 72 | **e + C2H2 → C2 + H2** | scale_electron_impact | ✓ VALID | Electron dissociation |
| 141 | **CH + CH → C2 + H2** | 2.16e-10 | ❓ CHECK | Baulch lists CH+CH→C2H2, not C2! |
| 176 | **C2H+ + e → C2 + H** | scale_recombination | ✓ VALID | Dissociative recombination |
| 204 | **C + CH → C2 + H** | 1.2e-10 | ✓ VALID | Baulch confirmed |
| 209 | **CH + C → C2 + H** | 1.2e-10 | ✓ DUPLICATE | Same as line 204 |
| 219 | **C2H2 + C → C2 + CH2** | 1.0e-10 | ⚠️ FIX RATE | Should be 2.0e-10 (Baulch) |
| 224 | **CH + C → C2 + H2** | 1.0e-10 | ✓ VALID | Different products than 204 |
| 226 | **CH2 + CH → C2 + H2 + H** | 1.2e-10 | ❌ REMOVE | NOT in Baulch ("garbage") |
| 243 | **CH + C → C2 + H** | 1.0e-10 | ✓ DUPLICATE | Same as 204/209 |
| 244 | **CH + CH → C2 + H2** | 1.0e-10 | ❓ CHECK | Baulch: CH+CH→C2H2 (not C2!) |
| 247 | **C2H + H → C2 + H2** | 1.0e-10 | ✓ VALID | Radical + H abstraction |
| 250 | **C2H2 + H → C2 + H2 + H** | 1.0e-11 | ❌ FIX RATE | 7.5e9× too high! |
| 254 | **C + C2H3 → C2 + CH3** | 1.0e-10 | ✓ VALID | C-atom insertion |
| 258 | **CH2 + CH2 → C2 + H2 + H2** | 1.0e-11 | ❌ FIX PRODUCTS | Should be C2H2 + H2 |
| 269 | **C + C + M → C2 + M** | three-body | ✓ VALID | Minor at low P |

---

## C2 DESTRUCTION PATHWAYS (Reactants contain C2)

| Line | Reaction | Rate Constant | Status | Notes |
|------|----------|---------------|--------|-------|
| 100 | **Ar* + C2 → Ar + C2** | 1.0e-10 | ✓ VALID | Quenching (elastic) |
| 206 | **C2 + H → CH + C** | 9.6e-11 | ❌ REMOVE | Endothermic! Not at 570K |
| 252 | **C2 + CH → C3 + H** | 1.0e-10 | ✓ VALID | C3 formation |
| 306 | **C2 → (wall)** | stick coeff | ✓ VALID | Wall loss |
| 347 | **C2 → (loss)** | diffusion | ✓ VALID | Diffusion loss |

---

## CRITICAL ISSUES IDENTIFIED

### 1. **CH + CH → C2 + H2** (Lines 141, 244)

**Problem:** Baulch (2005) lists **CH + CH → C2H2** (not C2!)

**Baulch Data:**
```
CH + CH → C2H2        k = 1.0e-10 cm³/s (ground state ³CH)
```

**Model has:**
```python
Line 141: CH + CH → C2 + H2    k = 2.16e-10  (ion reaction?)
Line 244: CH + CH → C2 + H2    k = 1.0e-10   (neutral reaction?)
```

**Need to verify:**
- Is line 141 for ions (CH+ + CH)?
- Line 244 should probably be CH + CH → C2H2 (not C2 + H2)

**Action:** Check if CH + CH can produce C2 directly, or only C2H2

---

### 2. **C2 + H → CH + C** (Line 206) - MAIN DESTROYER

**Problem:** This is **93% of C2 destruction**, but it's ENDOTHERMIC at 570K!

**Thermochemistry:**
```
C2 (6.21 eV bond) + H → CH + C (exothermic)
But: Reverse barrier is HIGH at low T
At 570K: kT = 0.05 eV << activation barrier
```

**Baulch Status:** NOT listed for T < 1000K (endothermic region)

**Action:** **REMOVE or set k=0 for T < 1000K**

---

### 3. **CH2 + CH → C2 + H2 + H** (Line 226)

**Problem:** NOT in Baulch, likely unphysical

**Action:** **REMOVE**

---

### 4. **CH2 + CH2 → C2 + H2 + H2** (Line 258)

**Problem:** Wrong products! Should make C2H2, not C2

**Baulch:**
```
³CH2 + ³CH2 → C2H2 + H2    k ~ 5.3e-11 cm³/s at 298K
```

**Action:** **FIX products: CH2 + CH2 → C2H2 + H2**

---

### 5. **H + C2H2 → C2 + H2 + H** (Line 250)

**Problem:** Rate 7.5 billion times too high!

**Correct (Baulch):**
```python
k = 1.67e-14 × T^1.64 × exp(-15250/T)  # cm³/s
At T=570K: k = 1.33e-21 cm³/s
```

**Model has:** k = 1.0e-11 cm³/s (constant)

**Action:** **Implement Arrhenius form**

---

### 6. **C2H2 + C → C2 + CH2** (Line 219)

**Problem:** Rate 2× too low

**Baulch:** k = 2.0e-10 cm³/s
**Model:** k = 1.0e-10 cm³/s

**Action:** **UPDATE to 2.0e-10**

---

## REAL C2 FORMATION AT LOW T (570K)

Based on Baulch (2005) and NIST database, VALID pathways at 570K:

### ✓ **Electron Impact**
```python
e + C2H2 → C2 + H2 + e         # Line 72 ✓
```

### ✓ **Ion Chemistry**
```python
C2H+ + e → C2 + H              # Line 176 ✓
```

### ✓ **C-Atom Reactions**
```python
C + CH → C2 + H                # Lines 204, 209, 243 ✓
C + CH → C2 + H2               # Line 224 ✓
C2H2 + C → C2 + CH2            # Line 219 (needs rate fix)
C + C2H3 → C2 + CH3            # Line 254 ✓
```

### ✓ **Radical Abstraction**
```python
C2H + H → C2 + H2              # Line 247 ✓
```

### ✓ **Three-Body (minor at low P)**
```python
C + C + M → C2 + M             # Line 269 ✓
```

### ❌ **NOT VALID at 570K**
```python
CH + CH → C2 + H2              # Lines 141, 244 (produces C2H2, not C2)
CH2 + CH → C2 + H2 + H         # Line 226 (doesn't exist)
CH2 + CH2 → C2 + H2 + H2       # Line 258 (produces C2H2, not C2)
H + C2H2 → C2 + H2 + H         # Line 250 (rate 7.5e9× wrong)
```

---

## CORRECTIONS TO IMPLEMENT

### Phase 1: Remove Invalid Reactions
```python
# define_rates.py
# Line 255: Comment out or set to zero
# k['C2_H_CH_C_cm3_7_6'] = 9.6e-11  # REMOVED: endothermic at 570K
k['C2_H_CH_C_cm3_7_6'] = 0.0  # Disabled for T < 1000K

# Line 288: Comment out
# k['CH2_CH_C2_H2_H_cm3_7_26'] = 1.2e-10  # REMOVED: not in literature
k['CH2_CH_C2_H2_H_cm3_7_26'] = 0.0  # Disabled
```

### Phase 2: Fix Rate Constants
```python
# define_rates.py
# Line 268: Fix C2H2 + C rate
k['C2H2_C_C2_CH2_cm3_7_19'] = 2.0e-10  # Was 1.0e-10 (Baulch 2005)

# Line 252: Fix CH + H rate
k['CH_H_C_H2_cm3_7_3'] = 2.0e-10  # Was 1.2e-10 (Baulch 2005)

# Line 315: Fix H + C2H2 rate (implement Arrhenius)
# k['C2H2_H_C2_H2_H_cm3_7_50'] = 1e-11  # OLD CONSTANT
k['C2H2_H_C2_H2_H_cm3_7_50'] = 1.67e-14 * Tgas**1.64 * np.exp(-15250/Tgas)
```

### Phase 3: Fix Stoichiometry
```python
# build_reactions.py
# Line 258: Fix products
# OLD: push(sto('CH2', 2), sto('C2', 1, 'H2', 2), ...)
# NEW: push(sto('CH2', 2), sto('C2H2', 1, 'H2', 1), ...)
```

### Phase 4: Verify CH + CH Reaction
- Check Baulch for CH + CH pathways
- Determine if line 141 is ionic or neutral
- Correct products if needed (C2H2 vs C2 + H2)

---

## EXPECTED IMPACT

### Before Corrections:
**C2 Production:**
- CH + CH2 → C2 - 21% ❌
- C2H2 + C → C2 - 18% (rate too low)
- CH2 + CH2 → C2 - 16% ❌
- H + C2H2 → C2 - 10% (rate 7.5e9× too high)

**C2 Destruction:**
- C2 + H → CH + C - 93% ❌

### After Corrections:
**C2 Production:**
- Lose 21% + 16% = 37% from fake pathways
- Gain 2× from C2H2 + C (18% → 36%)
- Lose ~99.9% from H + C2H2 (10% → ~0%)
- **Net production change: -11% (slight decrease)**

**C2 Destruction:**
- Lose 93% from C2 + H removal
- Remaining: wall (6%) + diffusion (1%) + C3 formation (minor)
- **Destruction decreases by ~15×!**

### Net Result:
**C2 should INCREASE by ~10-15× after corrections!**

Current: C2 = 2.75×10⁸ cm⁻³
Predicted: C2 ~ 3-4×10⁹ cm⁻³
Target: C2 = 5.60×10¹¹ cm⁻³

**Still ~150× too low, but MUCH better!**

---

## Next Steps

1. ✅ Document all C2 reactions
2. ⬜ Implement Phase 1 corrections (remove invalid reactions)
3. ⬜ Implement Phase 2 corrections (fix rates)
4. ⬜ Implement Phase 3 corrections (fix stoichiometry)
5. ⬜ Verify CH + CH reaction products
6. ⬜ Recalculate C2 pathways
7. ⬜ Analyze remaining gap (if any)

---
