# CRITICAL ISSUE: Charge Balance Violation

## Problem

**ALL optimization results violate charge neutrality requirement:**

| Result | Ne (cm⁻³) | Positive Ions (cm⁻³) | Ion/e Ratio | Imbalance |
|--------|-----------|---------------------|-------------|-----------|
| f=684  | 1.28×10⁹  | 9.91×10⁷            | 0.08x       | 96.5%     |
| f=814  | 1.95×10⁹  | 1.50×10⁸            | 0.08x       | 96.5%     |
| All others | ~10⁹   | ~10⁸                | 0.08-0.11x  | 93-96%    |

**Expected:** Ion density should be 1-10x electron density in sheath-negative glow transition region.

**Actual:** Ion density is 10-13x LOWER than electron density!

## Why This Happens

The simulation fixes `Ne` as an input parameter (1.28×10⁹ - 1.95×10⁹ cm⁻³), but:

1. **Ionization is too weak:** Electron-impact ionization rates don't produce enough ions to balance
2. **0D limitation:** No spatial transport to maintain quasineutrality
3. **Ion loss is high:** Wall recombination removes ions faster than they're created

### Ion Production vs. Loss

**Ion Sources:**
- e + Ar → Ar⁺ + 2e (main source)
- e + CH₄ → CH₄⁺ + 2e
- Secondary ionization from metastables

**Ion Sinks:**
- Wall recombination (very fast)
- Ion-electron recombination
- Charge exchange with neutrals

**Result:** Ion density reaches ~10⁸ cm⁻³ steady state, but we're forcing Ne = 10⁹ cm⁻³.

## Physical Consequences

A plasma with such severe charge imbalance **cannot physically exist:**

1. **Electric fields:** Would be enormous (~kV/cm) to separate charges
2. **Plasma oscillations:** Would quickly restore neutrality
3. **Sheath formation:** Requires ion density ≈ electron density in bulk

**Our simulation is unphysical in current state.**

## Why Optimization Still "Worked"

The optimizer found parameters that minimized the objective function for H, CH, C2 densities, BUT:
- Did not check charge balance
- Used unphysical electron density
- Neutral chemistry (H, CH, C2) doesn't directly depend on charge balance

**The neutral densities might be approximately correct even with wrong charge balance!**

## Solutions

### Option 1: Calculate Ne from Charge Neutrality (RECOMMENDED)
```python
# Instead of fixing Ne, calculate it:
Ne = sum(positive_ion_densities) - sum(negative_ion_densities)

# Where negative ions include HMinus, CH3Minus, etc.
```

**Pros:**
- Automatically enforces quasineutrality
- Physically consistent
- Matches how real plasmas behave

**Cons:**
- Ne becomes a result, not an input
- Need to verify against experimental Ne = 3.3×10⁹ ± 50%

### Option 2: Increase Ionization Rates
Adjust electron-impact ionization cross-sections until:
```
sum(positive_ions) ≈ Ne (within 10%)
```

**Pros:**
- Keeps Ne as input (matches experiment)
- Simple to implement

**Cons:**
- May exceed literature bounds for ionization rates
- Doesn't address fundamental 0D limitation

### Option 3: Add Ionization Constraint to Optimization
Add penalty to objective function:
```python
charge_penalty = abs(sum(positive_ions) - Ne) / Ne
objective = species_error + 1000 * charge_penalty
```

**Pros:**
- Forces optimizer to find charge-balanced solutions
- Can still control Ne range

**Cons:**
- May not find solution if chemistry can't produce enough ions
- Might require unrealistic rate constants

### Option 4: Accept 0D Model Limitation
Acknowledge that sheath-region physics cannot be accurately represented in 0D:
- Use simulation only for neutral species (H, CH, C2, hydrocarbons)
- Don't compare electron/ion densities to experiment
- Note that real plasma would have different Ne

**Pros:**
- Neutral densities may still be valid
- Avoids forcing unphysical balance

**Cons:**
- Cannot validate full model against experiment
- Missing important physics

## Recommendation

**Immediate:** Implement Option 1 (calculate Ne from charge balance)

**Steps:**
1. Modify ODE solver to compute Ne from ion balance at each timestep
2. Compare resulting Ne to experimental 3.3×10⁹ cm⁻³
3. If too low, increase ionization rates within literature bounds
4. If still can't reach experimental Ne, use Option 4 (accept limitation)

**Long-term:** Consider 1D fluid model that includes:
- Spatial transport (ambipolar diffusion)
- Sheath physics
- Electric field self-consistency

This would properly capture the physics of the sheath-negative glow transition region.

## Impact on Current Results

**Good news:** Neutral species densities (H, CH, C2) are probably still approximately correct because:
- Dominated by neutral-neutral reactions
- Electron-impact reactions contribute <15% to H/CH/C2 production
- Main chemistry (C2H2 + H → C2, CH2 + H → CH, etc.) doesn't involve charges

**Bad news:** Cannot validate electron density or fully validate model without fixing charge balance.

## Next Steps

1. ✓ Document issue (this file)
2. [ ] Implement Ne calculation from charge balance
3. [ ] Re-run best optimization with corrected Ne
4. [ ] Compare corrected Ne to experiment (should be 3.3×10⁹ ± 50%)
5. [ ] If needed, adjust ionization rates to match experimental Ne
6. [ ] Create final optimized parameters with proper charge balance
