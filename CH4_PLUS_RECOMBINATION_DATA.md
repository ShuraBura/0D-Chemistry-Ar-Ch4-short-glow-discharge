# CH4+ Dissociative Recombination - Literature Data

## Source
Thomas et al., "Dissociative Recombination of CH4+", J. Phys. Chem. A 2013, 117, 9999-10005
DOI: 10.1021/jp400353x

## Rate Coefficient

**Temperature-dependent formula**:
```
k(Te) = 1.71 × 10⁻⁶ (Te/300)^(-0.66) cm³/s
```

Valid for electron temperatures: 10 ≤ Te ≤ 1000 K

**At 300 K (room temperature)**:
```
k = 1.71 × 10⁻⁶ cm³/s
```

**Uncertainty**: ±0.02 × 10⁻⁶ cm³/s (coefficient) and ±0.02 (exponent)

## Product Channels (at low collision energies)

| Products | Branching Fraction |
|----------|-------------------|
| CH4 | 0.00 ± 0.00 |
| CH3 + H | 0.18 ± 0.03 |
| CH2 + 2H | 0.51 ± 0.03 |
| CH2 + H2 | 0.06 ± 0.01 |
| CH + H2 + H | 0.23 ± 0.01 |
| CH + 2H2 | 0.02 ± 0.01 |

**Key finding**: Two or more C–H bonds are broken in ~80% of all collisions!

## Dominant Channel for Our Model

For simplified chemistry modeling, use the dominant channel:
```
CH4+ + e → CH2 + 2H    (51% branching ratio)
```

Or combine similar channels:
```
CH4+ + e → CH3 + H     (18%)  
CH4+ + e → CH2 + H2    (57% total: CH2 + 2H + CH2 + H2)
CH4+ + e → CH + H2 + H (23%)
```

## Implications for Our Model

- **Current state**: CH4+ is 38.4% of all ions with NO recombination!
- **With this reaction**: Will create major electron sink
- **Expected impact**: Electron imbalance will reduce from 97.8% → ~40-50%
- **Ion/e ratio**: Should improve from 0.52x → 2-3x

## Implementation Strategy

For V6 optimization, we'll add:
1. Primary reaction: `CH4Plus_e_CH2_2H_cm3_6_XX` with k ~ 1.7 × 10⁻⁶ cm³/s
2. Make it tunable within literature bounds: [1.0e-6, 3.0e-6] cm³/s (±50% for uncertainty)
3. Alternative: Add multiple channels with branching ratios
