# CSB (Cathode Sheath Boundary) Region Model

## Decision: Start with CSB Region

**Date**: 2025-10-30

### Rationale for CSB Selection

**Advantages of CSB:**
1. **Narrower parameter ranges**: Te, Ne, and E are better constrained
2. **Clear charge balance target**: ni ~ 3-7×ne (measurable validation criterion)
3. **Better-defined plasma conditions**: More thermal, easier to characterize

**Trade-offs:**
- H is mostly diffused from CG (not locally produced)
- Need to handle H as an input/boundary condition

**Advantages of CG (not chosen for now):**
- H is produced locally (chemistry-dominated)
- Better suited for pure 0-D chemistry model

**Why we chose CSB despite H diffusion:**
- Tighter constraints make validation more tractable
- Charge balance provides clear physics check
- Can iterate faster with well-defined parameters
- Once CSB chemistry is validated, can extend to CG

---

## CSB Parameter Ranges

**Note**: CSB = SE (Sheath Edge) in experimental targets

### From Experimental Targets (SE/CSB)
```python
params_csb = {
    'P': 0.4,              # Torr (fixed)
    'Tg': 400,             # K (assumed)
    'ne': 2e10,            # cm⁻³ (from EXPERIMENTAL_TARGETS.md)
    'Te': 0.8,             # eV (from EXPERIMENTAL_TARGETS.md)
    'E_field': 50,         # V/cm (from EXPERIMENTAL_TARGETS.md)
    'L_discharge': 0.45,   # cm
}
```

**Target densities (experimental):**
- H:  6.35e14 cm⁻³
- CH: 9.27e8 cm⁻³
- C₂: 5.56e11 cm⁻³

**Key ratios:**
- H/CH = 6.85e5
- C₂/CH = 600
- H/C₂ = 1.14e3

---

## Validation Criteria

### 1. Charge Balance (Primary)

**Target**: ni ~ 3-7×ne

**Implementation**:
```python
# Calculate total positive ion density
ni_total = sum([n[species] for species in positive_ions])

# Calculate charge balance ratio
charge_ratio = ni_total / ne

# Validation check
is_valid = 3.0 <= charge_ratio <= 7.0
```

**Why this matters:**
- CSB is quasineutral but with ambipolar effects
- Ion density exceeds electron density due to diffusion
- Clear experimental constraint from measurements

### 2. Species Density Ratios (Secondary)

Focus on **ratios** rather than absolute values:
- C₂/CH ~ 600 (C₂ much more stable than CH)
- H/CH ~ 6.85e5 (H dominates)

### 3. Physics Consistency

- Electron heating: E-field + collisions
- Ion losses: ambipolar diffusion
- Neutral chemistry: H, C, C₂ balance

---

## Implementation Plan

### Phase 1: Baseline CSB Model
1. Set CSB parameters (Te, ne, E from user input)
2. Run steady-state simulation
3. Check charge balance: is 3 ≤ ni/ne ≤ 7?
4. Record species ratios

### Phase 2: Parameter Sweep
If charge balance fails:
- Sweep Te (±20% around nominal)
- Sweep ne (±30% around nominal)
- Sweep E-field (±20% around nominal)
- Find parameter sets satisfying ni ~ 3-7×ne

### Phase 3: Chemistry Validation
Once charge balance is satisfied:
- Compare C₂/CH, H/CH ratios to experiment
- Check if H density is reasonable (may be set by diffusion)
- Validate electron heating and power balance

---

## Open Questions

1. **What are the specific Te, ne, E ranges for CSB?**
   - Need experimental bounds or estimates

2. **How should we handle H diffusion from CG?**
   - Option A: Fix H density at target (6.35e14)
   - Option B: Let H evolve, compare ratio to target
   - Option C: Add H flux term (requires CG coupling)

3. **What is the negative ion density?**
   - CH3⁻, H⁻ contribute to charge balance
   - Need to include in ni calculation?

---

## Next Steps

**Immediate:**
1. Define CSB parameter ranges (user input)
2. Clarify charge balance definition: ni = sum(all positive ions)?
3. Decide how to handle negative ions in charge balance

**Once parameters defined:**
1. Create `run_csb_baseline.py` script
2. Implement charge balance validation
3. Run parameter sweep
4. Document results

---

## Notes

- CSB selection prioritizes **validation over perfect physics**
- Tighter constraints → faster iteration
- Can extend to CG once CSB is validated
- Charge balance is the key differentiator vs. CG approach
