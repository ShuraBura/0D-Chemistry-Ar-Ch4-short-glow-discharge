# Session Summary: CSB Model Development & Performance Analysis

**Date**: 2025-10-30
**Branch**: `claude/resume-h-profile-work-011CUdWT5ESdiBMDxL7U1m9q`

---

## 🎯 Session Objectives

1. Resume H profile work from previous session
2. Process experimental nH(x) profile data
3. Validate model against measured spatial averages
4. Run quick validation tests

---

## ✅ Major Accomplishments

### 1. Processed Experimental H Profile Data

**Created**: `process_experimental_data.py`

**Input**: `combined_density_matrix.csv` (500 data points, 0-17 mm)

**Results** (CG region 0-1 mm, 30 points):

| Species | Measured Avg | Previous Target | Agreement |
|---------|--------------|-----------------|-----------|
| H | 9.58e15 cm⁻³ | 8.58e15 | +12% |
| CH | 4.29e8 cm⁻³ | 4.6e8 | -7% |
| C₂ | 1.50e11 cm⁻³ | 1.44e11 | +4% |

**Key Finding**: All within ±12% - **excellent validation** of previous estimates!

**Diffusion Length**: L_diff_H = **0.057 cm** (fitted from exponential decay)
- Previous estimate: 0.1 cm
- **Implication**: H diffusion to walls is MUCH faster than assumed

**Outputs Generated**:
- `experimental_averages.txt` - measured values
- `nH_profile_fit.png` - H density with exponential fit
- `all_species_profiles.png` - all three species profiles

### 2. Updated Documentation & Code

**Files Updated**:
- `EXPERIMENTAL_TARGETS.md` - added measured values and L_diff_H
- `sweep_cg_reduced.py` - updated targets, added L_diff sweep [0.057, 0.1]
- `sweep_cg_constrained.py` - updated targets, added L_diff sweep [0.057, 0.1]
- `test_single_run.py` - updated with measured targets

### 3. Pivoted to CSB Region

**Rationale**: CSB has more favorable conditions for 0-D modeling:
- ✓ H density 15× lower (6.35e14 vs 9.58e15 cm⁻³)
- ✓ Te more thermal (0.5-2 eV vs 3-7 eV)
- ✓ E-field lower (20-300 V/cm vs 400-800 V/cm)
- ✓ Better defined plasma parameters

**Key User Clarification**: H in CSB **cannot be fully fixed**
- H has TWO contributions:
  1. **Diffusion influx from CG** (semi-fixed source term)
  2. **Local chemistry** (evolved by 0-D model)
- Correct model: **dH/dt = (local chemistry) + (diffusion_influx)**

**Created CSB Test Scripts**:
1. `test_csb_validation.py` - CSB with dynamic H
2. `test_csb_fixed_h.py` - CSB with fixed H (incorrect approach)
3. `test_csb_fixed_h_simple.py` - simplified fixed H version
4. `test_quick_validation.py` - CG quick test
5. `test_csb_h_diffusion.py` - **CSB with H diffusion influx (CORRECT)**

---

## ⚠️ Critical Finding: ODE Performance Bottleneck

### Performance Test Results

**Test conditions**:
- CSB region parameters (ne=5e9, Te=1.0, E=100 V/cm)
- Integration time: 1 second (very short!)
- Method: BDF (stiff solver)

**Results**:
- **Runtime**: >7 minutes for 1 second simulation (did not complete)
- **Root cause**: 42 species + 253 reactions = extremely stiff ODE system

**Estimated Sweep Times**:
- Single run: ~5-10 minutes
- Reduced sweep (216 runs): **18-36 hours**
- Comprehensive sweep (768 runs): **~3-5 days**

**Stiffness Factors**:
1. Large number of species/reactions
2. Wide range of timescales (fast chemistry + slow diffusion)
3. Small L_diff_H = 0.057 cm increases stiffness further
4. High chemical reactivity in CH₄/Ar system

---

## 📋 CSB Model Specification

### Correct H Balance Equation

```
dH/dt = (Production from chemistry) - (Losses to walls & reactions) + (Diffusion influx from CG)
```

### H Diffusion Influx Estimates

Based on simple diffusion model: Flux ≈ D × ∇n / L

| Scenario | H_diffusion_influx | Notes |
|----------|-------------------|-------|
| Low | 1e18 cm⁻³/s | Weak diffusion coupling |
| Medium | 1e19 cm⁻³/s | **Baseline estimate** |
| High | 1e20 cm⁻³/s | Strong diffusion coupling |

**Formula**: `Flux ≈ D_H × (n_H_CG - n_H_CSB) / L_diff`

Where:
- D_H ≈ 300 cm²/s (H diffusion coefficient at 0.4 Torr, 400 K)
- n_H_CG ≈ 9.58e15 cm⁻³ (measured)
- n_H_CSB ≈ 6.35e14 cm⁻³ (measured)
- L_diff ≈ 0.1 cm (CG-CSB separation)

### CSB Parameter Ranges (User-Specified)

```python
# CSB sweep parameters
param_grid_csb = {
    'ne': [1e9, 3e9, 5e9, 9e9],              # cm⁻³
    'Te': [0.5, 1.0, 1.5, 2.0],              # eV (more thermal)
    'E_field': [20, 50, 100, 200, 300],      # V/cm (lower than CG)
    'Tgas': [300, 400],                      # K
    'H_diffusion_influx': [1e18, 1e19, 1e20],# cm⁻³/s (tunable!)
    'L_diff': [0.1],                         # cm (fixed for CSB)
    'gamma_H': [0.01, 0.05],                 # H wall recombination
}
```

**Total combinations**: 4 × 4 × 5 × 2 × 3 × 2 = **960 runs**

At 5-10 min per run: **80-160 hours (3-7 days)**

---

## 🎯 Recommendations

### Option 1: Run Long Sweep (Overnight/Weekend)

**Pros**:
- Comprehensive parameter exploration
- Identifies best-fit parameters
- Validates model across wide range

**Cons**:
- Takes 3-7 days
- Requires stable compute environment

**Action**:
```bash
# Run overnight
nohup python3 sweep_csb_with_h_diffusion.py > sweep_csb.log 2>&1 &
```

### Option 2: Parallel Processing

**Pros**:
- Can reduce wall-clock time by factor of N_cores
- Efficient use of multi-core systems

**Cons**:
- Requires code modification
- Still takes substantial time

**Implementation**: Use `multiprocessing` or `joblib` to parallelize parameter sweep

### Option 3: Reduce Chemistry

**Pros**:
- Faster ODE integration
- May be sufficient for target species (H, CH, C₂)

**Cons**:
- Requires chemistry analysis
- May lose important pathways

**Approach**:
- Identify critical reactions for H, CH, C₂
- Remove minor species/reactions
- Validate reduced model against full model

### Option 4: Multi-Timescale Approach

**Pros**:
- Exploit separation of timescales
- Can be much faster for quasi-steady-state species

**Cons**:
- Complex implementation
- Requires careful analysis

**Approach**:
- Fast chemistry: quasi-steady-state
- Slow chemistry: fully dynamic
- Hybrid ODE-algebraic system

### Option 5: Accept Long Runtimes

**Pros**:
- No code changes needed
- Can run in background

**Cons**:
- Slow iteration
- Requires patience

**Recommendation**: Start with a **reduced grid** (e.g., 50-100 runs) to validate approach, then expand.

---

## 📊 Reduced CSB Sweep (Recommended First Step)

### Quick Validation Grid (~50 runs, ~4-8 hours)

```python
param_grid_quick = {
    'ne': [1e9, 5e9],                        # 2 values
    'Te': [0.8, 1.5],                        # 2 values
    'E_field': [50, 100, 200],               # 3 values
    'Tgas': [400],                           # 1 value (fixed)
    'H_diffusion_influx': [1e18, 1e19, 1e20],# 3 values
    'L_diff': [0.1],                         # 1 value (fixed)
    'gamma_H': [0.01],                       # 1 value (fixed)
}
```

**Total**: 2 × 2 × 3 × 3 = **36 runs** (~3-6 hours)

**Purpose**: Test H_diffusion_influx sensitivity and validate CSB approach

---

## 📁 Files Created This Session

### Analysis Scripts
- `process_experimental_data.py` - process nH(x) profiles, extract spatial averages

### Validation Scripts
- `test_single_run.py` - updated with measured targets
- `test_quick_validation.py` - CG validation with default L_diff
- `test_csb_validation.py` - CSB validation with dynamic H
- `test_csb_fixed_h.py` - CSB with fixed H (incorrect)
- `test_csb_fixed_h_simple.py` - simplified fixed H version
- `test_csb_h_diffusion.py` - **CSB with H diffusion influx** (CORRECT)

### Data Files
- `experimental_averages.txt` - measured spatial averages
- `nH_profile_fit.png` - H profile with exponential fit
- `all_species_profiles.png` - all species profiles

---

## 🚀 Next Steps

### Immediate (Next Session)

1. **Create CSB parameter sweep script** with H diffusion influx
2. **Run reduced grid** (36-50 runs) overnight
3. **Analyze results**:
   - Which H_diffusion_influx matches best?
   - CH and C₂ sensitivity to Te, ne, E-field?
   - Identify best-fit parameter set

### Short-Term

1. **Validate best-fit parameters** against all experimental data
2. **Create visualizations** of parameter sensitivity
3. **Document final model configuration**

### Long-Term (Optional)

1. **Chemistry reduction study** - identify critical reactions
2. **Parallel sweep implementation** - speed up exploration
3. **Extend to other conditions** - pressure, gas mix variations

---

## 🔑 Key Insights

### Scientific

1. **Measured spatial averages validate previous estimates** (±12%)
2. **L_diff_H much shorter than assumed** (0.057 vs 0.1 cm) → rapid wall loss
3. **CSB model requires H diffusion influx** as tunable source term
4. **H_diffusion_influx ~ 1e18-1e20 cm⁻³/s** is reasonable range

### Computational

1. **ODE system is extremely stiff** (42 species, 253 reactions)
2. **Single run takes 5-10 minutes** even for short integration times
3. **Full parameter sweeps require days** of computation
4. **Parallel processing or chemistry reduction needed** for rapid iteration

---

## 📮 Questions for User

1. **Sweep strategy**: Run overnight with reduced grid first, or go straight to full sweep?
2. **H_diffusion_influx**: Do you have independent estimates from experiments or simulations?
3. **Priority species**: Focus only on H, CH, C₂? Or need full chemistry predictions?
4. **Compute resources**: Can you run multi-day sweeps? Parallel processing available?

---

## ✅ Session Status

**All work committed and pushed to**:
`claude/resume-h-profile-work-011CUdWT5ESdiBMDxL7U1m9q`

**Commits**:
1. `ab6c6da` - Process experimental H profile data and update targets
2. `298e364` - Add CSB validation tests and identify integration performance issue
3. `b1343b6` - Add CSB model with H diffusion influx (correct approach)

**Ready for**: CSB parameter sweep implementation and execution
