# H2+ Species Setup Instructions

## Summary

H3+ chemistry has been added to the model, but **H2Plus must be added to species lists** in all configuration files.

## What Was Added

### Rate Coefficients (COMPLETE ✓)
1. **Group 2**: `e + H2 → H2+ + 2e` (ionization, Te-dependent)
2. **Group 5**:
   - `H2+ + H2 → H3+ + H` (fast formation, k=2e-9)
   - `H3+ + CH4 → CH5+ + H2` (proton transfer, k=1.5e-9)
   - `H3+ + H2 → H2+ + H2` (reverse, k=6.4e-10)
3. **Group 6**:
   - `H2+ + e → H + H` (recombination, Te-dependent)
   - `H3+ + e → H2 + H` (dominant recombination, Te-dependent)
   - `H3+ + e → H + H + H` (3-body recombination, Te-dependent)
4. **Group 9**: `H2+` wall sticking
5. **Group 10**: `H2+` drift loss

### Files Updated (COMPLETE ✓)
- `define_rates.py`: All H2+/H3+ rates added
- `define_rates_tunable.py`: All H2+/H3+ rates added with Te-dependence
- `build_reactions.py`: All H2+/H3+ reactions added to reaction builder

## What Needs To Be Done

### Add H2Plus to Species Lists

You need to add `'H2Plus'` to the species list and ion_species list in **ALL** configuration files:

#### Files to Update:
1. `main.py` (lines 79-85)
2. `run_optimized.py`
3. `run_cg_optimized.py`
4. `test_quick.py`
5. `parameter_sweep_cg.py`
6. `sweep_cg_reduced.py`
7. `sweep_cg_constrained.py`
8. `benchmark.py`

#### Where to Add H2Plus:

**Species List** (add after 'H3Plus'):
```python
'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
            'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
            'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'H2Plus', 'C2H3', 'C3H2', 'CHPlus', ...
```

**Ion Species List** (add after 'H3Plus'):
```python
'ion_species': ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH3Minus',
                'H3Plus', 'H2Plus', 'CHPlus', 'CH2Plus', 'C2H5Plus', ...
```

**Mobilities Dict** (add H2Plus mobility):
```python
'mobilities': {
    'ArPlus': 3057.28,
    'CH4Plus': 6432,
    'CH3Plus': 4949.6,
    'CH5Plus': 4761.6,
    'ArHPlus': 2969.6,
    'CH2Plus': 4949.6,
    'C2H5Plus': 4949.6,
    'C2H4Plus': 4949.6,
    'C2H3Plus': 4949.6,
    'C2HPlus': 5000,
    'H3Plus': 5000,
    'H2Plus': 6000,  # <-- ADD THIS (similar to H3+, molecular ion)
    'CHPlus': 5000,
    'CH3Minus': 2000,
    'HMinus': 5000,
}
```

**Mobility Rationale**: H2+ is a light molecular ion, mobility ~6000 cm²/(V·s) at 0.4 Torr
- Similar to H3+ (5000)
- Higher than heavier ions (CH4+, CH3+ ~5000-6400)

## Why H2+ Is Critical

Without H2+, there is **NO pathway to form H3+**:
- H3+ exists in your species list
- H3+ has drift and stick reactions
- But H3+ has ZERO formation!

The formation chain is:
```
e + H2 → H2+ + 2e         (ionization)
H2+ + H2 → H3+ + H        (very fast, k=2e-9)
```

Without H2+ ionization, H3+ cannot form, making it an orphaned species.

## Testing After Adding H2Plus

After updating all files, test with:
```bash
python3 main.py
```

Check that:
1. H2+ density is non-zero (should be ~1e7-1e9 cm⁻³)
2. H3+ density increases (was zero before, should now be ~1e6-1e8 cm⁻³)
3. No errors about missing species

## Impact on Results

**Expected changes**:
- H2+ will appear at low density (~1% of ArPlus)
- H3+ will now form properly (was impossible before)
- H3+ + CH4 → CH5+ + H2 provides additional CH5+ formation pathway
- Overall ion chemistry becomes more realistic

## Note on Tg=570K

Your Tg=570K is significantly above the 400K calibration point. Consider:
1. Group 7 (neutral-neutral) rates may be faster at 570K
2. Recommended: Implement Tg-dependence for Group 7 in future
3. For now: H3+ chemistry addition is the priority

---

**NEXT STEP**: Add H2Plus to species lists in all configuration files listed above.
