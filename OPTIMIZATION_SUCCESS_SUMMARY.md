# Optimization Success Summary

## Major Achievement

**Successfully reduced CH error from 46x to 6.84x** - an 85% improvement!

```
Progress Timeline:
  Baseline (Ne=3.3e9, uncorrected rates):     CH = 46x too high,   f(x) = 40,396
  2-param optimization (Ne+E only):           CH = 44x too high,   f(x) = 37,225
  32-param warm-start (FINAL):                CH = 6.84x too high, f(x) = 684  ✓
```

## Final Optimized Parameters

### Physical Parameters
- **Ne**: 1.284×10⁹ cm⁻³ (down 61% from experimental 3.3×10⁹)
- **E_field**: 300 V/cm (at maximum boundary, up from 50 V/cm)

### Target Species Results
| Species | Current | Target | Ratio | Status |
|---------|---------|--------|-------|--------|
| **H**   | 2.52×10¹³ | 5.18×10¹³ | 0.49x | Need +2x |
| **CH**  | 6.84×10⁹  | 1.00×10⁹  | **6.84x** | **Down from 46x!** |
| **C2**  | 5.11×10¹⁰ | 1.30×10¹¹ | 0.39x | Need +2.5x |

## Chemistry Analysis - How CH Was Reduced

### Baseline (CH = 46x too high):
**C2 + H → CH + C**: Contributed 59% of CH production
- C2 was 5.6x too high, feeding excessive CH production
- C2H2 + H → C2: 91% of C2 production (feedback loop)

### Optimized (CH = 6.84x too high):
**CH Production breakdown:**
1. CH2 + H → CH + H2: 48.7% (neutral chemistry)
2. **C2 + H → CH + C**: 24.0% (reduced from 59%!)
3. C + H → CH: 12.9%
4. Electron-impact from CH4: 14.4%

**Key: C2 was reduced by maximizing C2 removal:**
- C2 wall sticking (stick_C2_9_9): 58.6% of C2 loss (set to MAX 6250)
- C2 volume loss (loss_C2_11_3): 18.7% of C2 loss (set to MAX 2000)
- C2 + H → CH: 22.6% of C2 loss (still feeding CH, but less C2 available)

**Strategy:** Instead of trying to reduce C2 production (C2H2 + H → C2 is 96.5% and hard to change),
the optimizer maximized C2 removal rates to starve the C2 + H → CH pathway.

## Key Optimized Rates

### Maximized Loss Rates (reduce intermediates):
```
loss_C2_11_3:        2000.0    (MAX - volume loss of C2)
loss_CH_11_9:        10000.0   (MAX - volume loss of CH)
stick_C2_9_9:        6250.0    (MAX - wall sticking of C2)
stick_CH_9_3:        6250.0    (MAX - wall sticking of CH)
stick_H_9_1:         3890.0    (HIGH - wall sticking of H)
stick_C_9_10:        6250.0    (MAX - wall sticking of C)
stick_CH3_9_2:       5820.0    (HIGH - wall sticking of CH3)
stick_C2H2_9_11:     2000.0    (MAX - wall sticking of C2H2)
```

### Minimized Production Rates:
```
e_CH4_CH_H2_H_vib_cm3_1_3:      2e-11  (MIN - electron-impact CH production)
ArStar_CH4_CH3_H_cm3_3_1:       5e-11  (MIN - metastable dissociation)
e_Ar_ArStar_cm3_1_7:            4e-11  (LOW - Ar metastable production)
stick_ArStar_9_5:               71.4   (MIN - keep metastables low)
```

## Optimization Details

- **Algorithm**: Differential evolution with warm start
- **Parameters optimized**: 32 (30 reaction rates + E_field + Ne)
- **Iterations**: 15 (completed successfully)
- **Evaluations**: 558
- **Time**: 3.7 minutes
- **Objective function**: f(x) = 684.19

### Warm Start Strategy (based on f(x)=1097 from step 4):
- Ne biased toward 1.3-2.0×10⁹ (lower end)
- E_field biased toward 200-280 V/cm (higher end)
- C2H2→C2 and C2→CH rates biased toward minimums

## What Worked

1. **High E field (300 V/cm)**: Enhanced electron-impact dissociation, created more radicals
2. **Low Ne (1.28×10⁹)**: Reduced overall ionization and ion-driven chemistry
3. **Maximized removal rates**: Wall sticking and volume loss for C2, CH, C
4. **Broke C2H2→C2→CH feedback**: By removing C2 faster than it could make CH

## Remaining Challenges

1. **H is too low** (0.49x target): Need to increase H production or reduce H loss
   - Currently: H drift gain = 3.2×10¹⁷ cm⁻³ s⁻¹
   - H wall sticking at 3890 (near max) - consuming H

2. **CH still 6.84x too high**: Further reduction needed
   - CH2 + H → CH (48.7%) and C2 + H → CH (24%) still dominant
   - May need to reduce CH2 and C2 further

3. **C2 is too low** (0.39x target): Overshot the reduction
   - Wall sticking and volume loss at maximums
   - May need to balance C2 removal vs. C2 target

## Next Steps to Consider

1. **Multi-objective optimization**: Balance H, CH, C2 simultaneously
   - Current weighting: H=1x, CH=20x, C2=3x
   - May need to adjust weights or use Pareto optimization

2. **Tune H drift gain**: Currently fixed at 3.2×10¹⁷
   - Could make this a tunable parameter
   - Or adjust H wall sticking rate

3. **Fine-tune Ne and E**:
   - E_field hit maximum boundary (300 V/cm)
   - Could expand range or lock at 300 and reoptimize

4. **Investigate CH2 chemistry**: Now the dominant CH source (48.7%)
   - CH2 + H → CH + H2
   - Reducing CH2 might help CH without affecting C2

## Files Generated

- `optimization_results/FINAL_RESULT.json`: Complete parameter set
- `optimization_results/best_iteration_0000_f684.2.json`: Detailed chemistry analysis
- `optimization_run.log`: Full optimization log

## Conclusion

The warm-start optimization was **highly successful**, achieving:
- **85% reduction in CH error** (46x → 6.84x)
- **Completion in 3.7 minutes** (558 evaluations)
- **Comprehensive logging** of all chemistry for further analysis

The optimizer discovered a smart strategy: **maximize C2 removal** to break the
C2H2 → C2 → CH feedback loop, rather than trying to prevent C2 production.

This demonstrates that the optimization framework is working correctly and finding
physically meaningful solutions within literature-constrained bounds.
