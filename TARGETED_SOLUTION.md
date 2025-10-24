# Targeted Solution to CH Problem

## Key Reactions to Target (Based on Chemistry Diagnostic)

### 1. **C2 + H → CH + C** (59% of CH production!)
```
Current rate: 9.6×10⁻¹¹ cm³/s
Literature range: [8×10⁻¹¹, 1.2×10⁻¹⁰]
Current position: Mid-range (can reduce to 8×10⁻¹¹)
Impact: 59% of CH production (2.42×10¹⁵ cm⁻³ s⁻¹)
```

**Action:** Reduce from 9.6×10⁻¹¹ → 8×10⁻¹¹ (20% reduction)
**Expected:** CH production reduces by ~12% (20% of 59%)

---

### 2. **C2H2 + H → C2 + H2 + H** (91% of C2 production!)
```
Current rate: 1×10⁻¹¹ cm³/s
Literature range: [8×10⁻¹², 1.2×10⁻¹¹]
Current position: Mid-range (can reduce to 8×10⁻¹²)
Impact: 91% of C2 production (3.57×10¹⁵ cm⁻³ s⁻¹)
```

**Action:** Reduce from 1×10⁻¹¹ → 8×10⁻¹² (20% reduction)
**Expected:** C2 production reduces by ~18% (20% of 91%)
- This cascades: less C2 → less CH production from C2+H reaction

---

### 3. **CH + CH → C2 + H2** (Could help consume CH)
```
Current rate: 2.16×10⁻¹⁰ cm³/s (ABOVE MAX!)
Literature range: [1.2×10⁻¹⁰, 1.8×10⁻¹⁰]
Current position: 20% above maximum
Impact: Only 0.01% of C2 production (insignificant)
```

**Action:** Reduce from 2.16×10⁻¹⁰ → 1.8×10⁻¹⁰ (bring to max)
**Why it doesn't help:**
- CH density is high (4.59×10¹⁰)
- But this reaction rate is tiny (3.79×10¹¹ cm⁻³ s⁻¹)
- Only 0.09% of CH loss (negligible)
- The problem: reaction is too slow to matter even with high CH

---

## Recommended Strategy

### Phase 1: Break the C2H2 → C2 → CH Loop

**Modify these 2 key rates to their minimum literature values:**

1. **C2H2_H_C2_H2_H_cm3_7_50**: 1.0×10⁻¹¹ → 8.0×10⁻¹² (20% reduction)
   - Reduces C2 production by ~18%

2. **C2_H_CH_C_cm3_7_6**: 9.6×10⁻¹¹ → 8.0×10⁻¹¹ (17% reduction)
   - Reduces CH production by ~10%

**Combined effect:**
- Less C2H2 converts to C2 (18% less C2)
- Less C2 converts to CH (10% direct + cascade from C2 reduction)
- **Expected total CH reduction: ~20-30%**

**From:** CH = 4.59×10¹⁰ (46x too high)
**To:** CH ≈ 3.2-3.7×10¹⁰ (32-37x too high)

Still not enough, but **major progress!**

---

### Phase 2: Increase CH Wall Losses

Current CH wall loss: 11% of total CH destruction (4.59×10¹⁴ cm⁻³ s⁻¹)

**loss_CH_11_9:**
```
Current: 1.0×10⁴ s⁻¹ (at MAXIMUM)
Range: [1×10³, 1×10⁴]
```

**Problem:** Already at maximum!

**Alternative - increase C2 wall losses:**

**loss_C2_11_3:**
```
Current: 5.85×10¹⁴ cm⁻³ s⁻¹ (15% of C2 destruction)
Range: [1×10⁻⁴, 2×10³] s⁻¹
```

If we prevent C2 from accumulating, it can't make CH!

**Action:** Increase C2 wall loss from current to 2×10³ (max)
**Impact:** More C2 lost to walls → less available to make CH

---

### Phase 3: Multi-Parameter Optimization Targeting Loop

Run optimization targeting these specific reactions:

**Tunable rates (focus on loop):**
1. C2H2_H_C2_H2_H_cm3_7_50 [8e-12, 1.2e-11] - C2 formation
2. C2_H_CH_C_cm3_7_6 [8e-11, 1.2e-10] - CH formation
3. loss_C2_11_3 [1e-4, 2e3] - C2 removal
4. loss_CH_11_9 [1e3, 1e4] - CH removal
5. CH production reactions from CH4:
   - e_CH4_CH_H2_H_vib_cm3_1_3 [2e-11, 1e-10]
6. C2H2 production reactions:
   - CH_CH3_C2H2_H2_cm3_7_23 [8e-11, 1.2e-10]
   - CH3_CH3_C2H2_H2_H2_cm3_7_49

**Plus Ne and E as before**

---

## Expected Outcome from Phase 1 (Quick Test)

Let me test just reducing those 2 key rates:

| Parameter | Current | New | Change |
|-----------|---------|-----|--------|
| C2H2 + H → C2 | 1.0×10⁻¹¹ | 8.0×10⁻¹² | -20% |
| C2 + H → CH | 9.6×10⁻¹¹ | 8.0×10⁻¹¹ | -17% |

**Run simulation with just these changes and see impact on:**
- C2H2 density
- C2 density (should drop ~18%)
- CH density (should drop ~20-30%)
- H density (might improve since less consumed by C2 reactions)

---

## Implementation Plan

### Immediate (5 minutes):
1. Modify define_rates.py to use minimum values for these 2 reactions
2. Run baseline simulation
3. Check if CH improves

### Short-term (30 minutes):
4. If Phase 1 helps, add C2 wall loss increase
5. Run again
6. Check progress

### If still needed (1-2 hours):
7. Run focused optimization on just the 6-8 key reactions in the loop
8. Keep Ne and E as additional parameters
9. See if targets achievable

---

## Quick Test Script

```python
# Modify just these 2 rates, run simulation, report CH

from define_rates import define_rates
# ... modify k['C2H2_H_C2_H2_H_cm3_7_50'] = 8e-12
# ... modify k['C2_H_CH_C_cm3_7_6'] = 8e-11
# Run simulation, check CH density
```

Would you like me to:
1. **Run Phase 1 quick test** - just reduce those 2 key rates to minimum?
2. **Create focused optimization** - target just the loop reactions?
3. **Something else?**
