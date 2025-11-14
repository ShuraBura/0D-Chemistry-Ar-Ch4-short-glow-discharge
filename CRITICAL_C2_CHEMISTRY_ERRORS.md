# CRITICAL C2 CHEMISTRY ERRORS FOUND

**Date:** 2025-11-14
**Status:** Chemistry model has fundamental errors in C2 formation/destruction

---

## Summary of Baulch (2005) Verification Results

We cross-checked the main C2 reactions against **Baulch et al. (2005)** - the gold standard for combustion kinetics. The results are shocking:

### ❌ **Top 3 C2 Production Pathways Are WRONG**

| Reaction | Model Rate | Baulch Status | % of C2 Prod. |
|----------|-----------|---------------|---------------|
| **CH + CH2 → C2 + H2 + H** | 1.2e-10 | ❌ NOT IN BAULCH ("garbage") | 21% |
| **C2H2 + C → C2 + CH2** | 1.0e-10 | ✓ Should be 2.0e-10 | 18% |
| **CH2 + CH2 → C2 + H2 + H2** | 1.0e-11 | ❌ WRONG PRODUCTS! | 16% |

### ❌ **Main C2 Destruction Is WRONG**

| Reaction | Model Rate | Baulch Status | % of C2 Dest. |
|----------|-----------|---------------|---------------|
| **C2 + H → CH + C** | 9.6e-11 | ❌ NOT IN BAULCH (endothermic) | 93% |

---

## Detailed Findings

### 1. **C2 + H → CH + C** - LINE 255 in define_rates.py

**Current Model:**
```python
k['C2_H_CH_C_cm3_7_6'] = 9.6e-11  # cm³/s
```

**Baulch Status:** ❌ **NOT FOUND**

**Why it's wrong:**
- This reaction is **endothermic** at low temperatures
- C2 + H → CH + C requires breaking the strong C≡C bond (6.21 eV)
- At T=570K, thermal energy kT = 0.05 eV << bond energy
- Baulch (2005) does NOT list this reaction for T < 1000K
- **This reaction should NOT occur at 570K!**

**Impact:**
- **93% of all C2 destruction** in current model
- Artificially suppressing C2 density by ~10×

**Action Required:**
- **REMOVE or set to zero** for T < 1000K
- Search for REAL C2 destruction pathways at low T

---

### 2. **CH + CH2 → C2 + H2 + H** - LINE 288 in define_rates.py

**Current Model:**
```python
k['CH2_CH_C2_H2_H_cm3_7_26'] = 1.2e-10  # cm³/s
```

**Baulch Status:** ❌ **NOT FOUND** ("garbage reaction")

**Why it's wrong:**
- Not listed in Baulch (2005)
- Not a known reaction in combustion chemistry
- Implausible products (3-body breakup without stabilization)

**Impact:**
- **21% of C2 production** in current model
- Overestimating C2 from this non-existent pathway

**Action Required:**
- **REMOVE this reaction**
- Replace with validated CH + CH2 pathways from literature

---

### 3. **CH2 + CH2 → C2 + H2 + H2** - LINE 323 in define_rates.py

**Current Model:**
```python
k['CH2_CH2_C2_H2_H2_cm3_7_58'] = 1e-11  # cm³/s (produces C2!)
```

**Baulch Status:** ❌ **WRONG PRODUCTS!**

**Correct Reaction (Baulch 2005):**
```python
³CH2 + ³CH2 → C2H2 + H2  # Produces C2H2, NOT C2!
```

**Baulch Rate:**
- k = 3×10⁻⁹ × exp(-6000/T) cm³/s for 1000-3000 K
- Or k = 5.3×10⁻¹¹ cm³/s at 298 K (other sources)

**Impact:**
- **16% of C2 production** assigned to wrong product
- Should produce C2H2, not C2

**Action Required:**
- **FIX products:** CH2 + CH2 → C2H2 + H2
- Update rate constant to literature value

---

### 4. **C2H2 + C → C2 + CH2** - LINE 268 in define_rates.py

**Current Model:**
```python
k['C2H2_C_C2_CH2_cm3_7_19'] = 1e-10  # cm³/s
```

**Baulch Status:** ✓ **IN BAULCH, but rate is wrong**

**Correct Rate (Baulch 2005):**
```python
k = 2.0e-10 cm³/s
```

**Impact:**
- **18% of C2 production**
- Current rate is **2× too low**

**Action Required:**
- **UPDATE:** k = 2.0e-10 cm³/s

---

### 5. **CH + H → C + H2** - LINE 252 in define_rates.py

**Current Model:**
```python
k['CH_H_C_H2_cm3_7_3'] = 1.2e-10  # cm³/s
```

**Baulch Status:** ✓ **IN BAULCH, but rate is wrong**

**Correct Rate (Baulch 2005):**
```python
k = 2.0e-10 cm³/s
```

**Impact:**
- CH destruction pathway
- Affects C atom production for C2H2 + C → C2 reaction

**Action Required:**
- **UPDATE:** k = 2.0e-10 cm³/s

---

### 6. **H + C2H2 → C2 + H2 + H** - LINE 315 in define_rates.py

**Current Model:**
```python
k['C2H2_H_C2_H2_H_cm3_7_50'] = 1e-11  # CONSTANT (WRONG!)
```

**Baulch Status:** ✓ **IN BAULCH, but rate is 7.5 billion times too high!**

**Correct Rate (Baulch 2005):**
```python
k = 1.67e-14 × T^1.64 × exp(-15250/T)  # cm³/s
```

**At T=570K:**
- Correct: k = 1.33×10⁻²¹ cm³/s
- Model: k = 1.00×10⁻¹¹ cm³/s
- **Error: 7.5 billion times too high!**

**Action Required:**
- **ALREADY IDENTIFIED** - needs Arrhenius form implementation

---

## Critical Implications

### Before Corrections:
The model predicted C2 pathways were:

**C2 Production:**
1. CH + CH2 → C2 + H2 + H - **21%** (doesn't exist!)
2. C2H2 + C → C2 + CH2 - **18%** (rate 2× too low)
3. CH2 + CH2 → C2 + H2 + H2 - **16%** (wrong products!)
4. H + C2H2 → C2 + H2 + H - **10%** (rate 7.5e9× too high!)

**C2 Destruction:**
1. C2 + H → CH + C - **93%** (doesn't exist at low T!)

### After Corrections:
- Top 3 C2 production pathways are eliminated or corrected
- Main C2 destruction pathway is eliminated
- **Need to find REAL C2 chemistry at 570K!**

---

## What REALLY Produces C2 at Low Temperatures?

### Known C2 Formation Pathways at T < 1000K:

Based on literature search (NIST, Baulch, combustion databases):

#### ✓ **Electron Impact Dissociation**
```python
e + C2H2 → C2 + H2 + e      # Already in model (line 110)
e + C2H4 → C2 + 2H2 + e     # Check if present
```

#### ✓ **Ion-Molecule Reactions**
```python
C2H+ + e → C2 + H           # Dissociative recombination (line 228)
```

#### ✓ **Radical Recombination**
```python
CH + CH → C2H2              # Baulch: k = 1.0e-10 (line 309, but products wrong!)
C + CH → C2 + H             # Baulch: k = 1.2e-10 (line 253) ✓
```

#### ✓ **C-atom Reactions**
```python
C2H2 + C → C2 + CH2         # Baulch: k = 2.0e-10 (line 268, needs fix) ✓
C2H + C → C2 + CH           # Check if present
```

#### ❓ **Possible Missing Pathways:**
```python
CH + C → C2 + H             # Barrierless radical recombination?
C2H + H → C2 + H2           # Already in model (line 312)
C + C + M → C2 + M          # Three-body recombination (unlikely at low P)
```

---

## Immediate Action Items

### Priority 1: Remove/Fix Incorrect Reactions
1. ❌ **REMOVE:** C2 + H → CH + C (line 255) - endothermic, doesn't occur at 570K
2. ❌ **REMOVE:** CH + CH2 → C2 + H2 + H (line 288) - not in literature
3. ✓ **FIX PRODUCTS:** CH2 + CH2 → C2H2 + H2 (line 323) - produces C2H2, not C2
4. ✓ **FIX RATE:** C2H2 + C → C2 + CH2 - change to 2.0e-10 (line 268)
5. ✓ **FIX RATE:** CH + H → C + H2 - change to 2.0e-10 (line 252)
6. ✓ **FIX RATE:** H + C2H2 → C2 - implement Arrhenius form (line 315)

### Priority 2: Search for Real C2 Pathways
1. Check electron-impact reactions for C2 formation
2. Verify ion chemistry produces C2
3. Search for missing neutral-neutral C2 formation
4. Check if CH + CH → C2H2 or CH + CH → C2 + H2

### Priority 3: Recalculate
1. Run model with corrections
2. Analyze new C2 pathways
3. Compare with experimental targets

---

## Expected Impact

**With these corrections:**
- C2 production will decrease (lose 21% + 16% fake pathways, gain 2× from C2H2+C)
- C2 destruction will decrease dramatically (lose 93% fake destruction)
- **Net effect:** C2 should INCREASE significantly

**Prediction:**
- Current C2 = 2.75×10⁸ cm⁻³ (2000× too low)
- After removing C2 + H destruction: C2 should increase by ~15×
- Need to verify REAL C2 formation pathways to match target

---

## Next Steps

**Choose one:**

**Option A:** Correct all errors NOW and recalculate
- Fast implementation (~30 min)
- See immediate impact
- May reveal more missing chemistry

**Option B:** Literature search FIRST for real C2 pathways
- Systematic approach
- Understand what SHOULD be there
- Then implement all at once

**Option C:** BOTH - correct known errors, then search for gaps
- Pragmatic approach
- Quick wins + thorough understanding

**RECOMMENDATION: Option C**

1. Fix the 6 identified errors (15 minutes)
2. Run model and analyze new C2 pathways (10 minutes)
3. Identify remaining gaps and search literature (30 minutes)
4. Implement missing chemistry if needed (30 minutes)

---

**Files to Modify:**
- `define_rates.py` - rate constants (lines 252, 255, 268, 288, 315, 323)
- `build_reactions.py` - reaction stoichiometry (CH2 + CH2 products)
- `species_list.py` - verify all species present

**Total Time Estimate:** ~90 minutes

---
