# Python Conversion of Ar/CH4 Plasma Simulation

This is a **Python port** of the MATLAB 0D plasma chemistry simulation, with significant optimizations for performance.

## Overview

- **0D chemistry model** for Argon-Methane glow discharge plasma
- **41 chemical species** (electrons, ions, neutrals, excited states)
- **238 reactions** (electron impact, ion-neutral, recombination, etc.)
- **Stiff ODE system** solved using SciPy's BDF solver

## Files

### Core Simulation
- `define_rates.py` - All reaction rate constants (280+ rates)
- `build_reactions.py` - Reaction network builder
- `odefun.py` - Original ODE function (loop-based)
- `odefun_optimized.py` - **Optimized ODE function** (vectorized, sparse matrices)
- `run_optimized.py` - **Main script** to run the simulation
- `main.py` - Alternative main script with plotting

### Testing & Benchmarking
- `test_quick.py` - Quick test (0.1 second simulation)
- `benchmark.py` - Performance comparison (original vs optimized)

## Quick Start

### Run the optimized simulation:
```bash
python3 run_optimized.py
```

### Run benchmark:
```bash
python3 benchmark.py
```

## Performance

### Optimized Version
- **0.40 seconds** to simulate 100 seconds of plasma chemistry
- **250x faster than real-time**
- **806 function evaluations**
- **3x faster** per function call vs loop-based version

### Key Optimizations
1. **Sparse CSR matrices** for stoichiometry → memory efficient
2. **Vectorized NumPy operations** → eliminated Python loops
3. **Pre-computed structures** → no repeated calculations during solving

## Results (t = 100 s)

### Key Radicals
- H: 3.48×10¹³ cm⁻³
- CH: 5.30×10¹⁰ cm⁻³
- C2: 7.59×10¹¹ cm⁻³
- CH3: 3.96×10¹³ cm⁻³

### Major Ions
- e⁻: 1.00×10¹⁰ cm⁻³ (fixed)
- Ar⁺: 4.54×10⁸ cm⁻³
- CH5⁺: 1.39×10⁹ cm⁻³

### Stable Molecules
- H2: 8.55×10¹³ cm⁻³
- CH4: 1.45×10¹⁵ cm⁻³ (fixed)
- C2H2: 1.07×10¹³ cm⁻³
- C2H6: 1.89×10¹³ cm⁻³

## Dependencies

```bash
pip install numpy scipy matplotlib
```

## Comparison: MATLAB vs Python

| Aspect | MATLAB | Python (Optimized) |
|--------|--------|-------------------|
| Language | MATLAB | Python 3 |
| ODE Solver | ode15s | scipy BDF |
| Performance | Baseline | ~3x faster per call |
| Memory | Dense arrays | Sparse matrices |
| Dependencies | MATLAB license | Free (NumPy/SciPy) |

## Technical Details

### Plasma Conditions
- Pressure: 400 mTorr
- Gas temperature: 400 K
- Electron temperature: 1 eV
- Electric field: 50 V/cm
- Ar:CH4 ratio: 85:15

### Numerical Method
- **Solver**: BDF (Backward Differentiation Formula)
- **Tolerances**: rtol=1e-6, atol=1e-7
- **Time span**: 0 → 100 seconds
- **Stiffness**: High (requires implicit solver)

## Future Optimizations

Possible further improvements:
- **Numba JIT compilation** → 5-10x additional speedup
- **Cython** for critical loops
- **JAX** for automatic differentiation and GPU support
- **Analytical Jacobian** (pre-computed symbolic)

## Author

Converted from MATLAB to Python with optimizations.
