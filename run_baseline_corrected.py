#!/usr/bin/env python3
"""
Baseline simulation with:
1. Corrected electron density: 3.3e9 cm^-3 (experimental)
2. All rates within literature bounds

This establishes the baseline before optimization.
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


def define_rates_corrected(params):
    """
    Load rates with all values corrected to be within literature bounds.
    """
    # Start with original rates
    k = define_rates(params)

    # Load literature database
    db = get_complete_rate_database()

    # Correct any out-of-bounds rates
    corrections = 0
    for name, rate_db in db.items():
        if name in k:
            if k[name] < rate_db.min:
                k[name] = rate_db.min
                corrections += 1
            elif k[name] > rate_db.max:
                k[name] = rate_db.max
                corrections += 1

    return k, corrections


def calculate_initial_densities(params):
    """Calculate initial species densities."""
    species = params['species']
    ns = len(species)
    y0 = np.ones(ns) * 1e3

    def set_density(name, value):
        try:
            idx = species.index(name)
            y0[idx] = value
        except ValueError:
            pass

    set_density('e', params['ne'])
    set_density('Ar', 0.85 * 9.66e15)
    set_density('CH4', 0.15 * 9.66e15)
    set_density('ArPlus', 1e7)
    set_density('CH4Plus', 1e5)
    set_density('CH3Plus', 1e5)
    set_density('CH5Plus', 1e3)
    set_density('ArHPlus', 5e5)
    set_density('CH3Minus', 5e4)
    set_density('H2', 1e12)
    set_density('ArStar', 5e6)
    set_density('H', 1e11)
    set_density('C2', 5e7)
    set_density('CH', 5e4)
    set_density('C2H4', 5e7)
    set_density('C2H6', 1e6)
    set_density('CH2', 1e11)
    set_density('C2H2', 1e12)
    set_density('C2H5', 1e6)
    set_density('CH3', 5e7)
    set_density('C', 5e7)

    return y0


print("=" * 80)
print(" BASELINE SIMULATION - CORRECTED RATES")
print("=" * 80)
print("\nChanges from previous simulation:")
print("  1. Electron density: 1.0e10 → 3.3e9 cm^-3 (experimental)")
print("  2. All 29 out-of-bounds rates corrected to literature limits")
print("     - CH+CH→C2+H2: 2.16e-10 → 1.80e-10 (20% reduction)")
print("     - 25 recombination rates: reduced to lit. max")
print("     - 4 other rates: brought within bounds")

# Setup parameters
params = {
    'E_field': 50,
    'L_discharge': 0.45,
    'ne': 3.3e9,  # EXPERIMENTAL VALUE
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'C2H2Star'],
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'CH3Minus': 3000, 'HMinus': 3000
    }
}

print("\nBuilding reaction network with corrected rates...")
params['k'], num_corrections = define_rates_corrected(params)
print(f"  ✓ Corrected {num_corrections} rates to be within literature bounds")

params['R'], params['tags'] = build_reactions(params)
print(f"  Species: {len(params['species'])}")
print(f"  Reactions: {len(params['R'])}")

y0 = calculate_initial_densities(params)

# Create optimized ODE function
print("\nInitializing optimized ODE solver...")
ode_func = PlasmaODE_Optimized(params)
print("  ✓ Sparse matrices built")
print("  ✓ Vectorization ready")

# Run simulation
print("\n" + "=" * 80)
print(" Running Baseline Simulation: 0 → 100 seconds")
print("=" * 80)

t_span = (0, 100)
t_eval = [0, 0.1, 1, 10, 100]

print(f"\nSolver: BDF (stiff)")
print(f"Tolerances: rtol=1e-6, atol=1e-7")
print("\nStarting...")

start_time = time.time()

sol = solve_ivp(
    ode_func,
    t_span,
    y0,
    method='BDF',
    t_eval=t_eval,
    rtol=1e-6,
    atol=1e-7
)

elapsed = time.time() - start_time

print("\n" + "=" * 80)
print(f" ✓ COMPLETED in {elapsed:.2f} seconds")
print("=" * 80)
print(f"  Status: {sol.message}")
print(f"  Function evaluations: {sol.nfev}")

# Display key results
print("\n" + "=" * 80)
print(" BASELINE RESULTS (t = 100 s)")
print("=" * 80)

# Get final densities
y_final = sol.y[:, -1]
species = params['species']

def get_density(name):
    try:
        idx = species.index(name)
        return y_final[idx]
    except ValueError:
        return 0.0

# Targets
targets = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11
}

print("\nTarget Species:")
print("-" * 80)
for sp in ['H', 'CH', 'C2']:
    current = get_density(sp)
    target = targets[sp]
    ratio = current / target
    status = "✓" if 0.8 <= ratio <= 1.2 else "✗"

    print(f"{sp:5s}: {current:.2e} cm^-3  (target: {target:.2e})")
    print(f"       Ratio: {ratio:6.2f}x  {status}")

print("\nOther Key Species:")
print("-" * 80)
for sp in ['e', 'CH3', 'C2H2', 'C2H4', 'H2']:
    dens = get_density(sp)
    print(f"{sp:5s}: {dens:.2e} cm^-3")

print("\n" + "=" * 80)
print(" COMPARISON TO PREVIOUS SIMULATION")
print("=" * 80)

# Previous results (with Ne=1e10, uncorrected rates)
prev = {
    'e': 1.0e10,
    'H': 3.4772e13,
    'CH': 5.3023e10,
    'C2': 7.5862e11
}

print("\n                Previous      Current       Change")
print("-" * 60)
for sp in ['e', 'H', 'CH', 'C2']:
    current = get_density(sp)
    previous = prev[sp]
    change = (current - previous) / previous * 100
    arrow = "↑" if change > 0 else "↓"
    print(f"{sp:5s}:      {previous:.2e}    {current:.2e}    {change:+6.1f}% {arrow}")

print("\n" + "=" * 80)
print(" ANALYSIS")
print("=" * 80)

ch_current = get_density('CH')
ch_target = targets['CH']
ch_ratio = ch_current / ch_target

print(f"\nCH Status:")
print(f"  Current: {ch_current:.2e} cm^-3")
print(f"  Target:  {ch_target:.2e} cm^-3")
print(f"  Ratio:   {ch_ratio:.1f}x too high")

if ch_ratio > 10:
    print(f"\n  ⚠ CH still severely elevated ({ch_ratio:.0f}x)")
    print("  → Electron density correction and rate fixes not sufficient")
    print("  → Need multi-parameter optimization")
elif ch_ratio > 2:
    print(f"\n  ⚠ CH moderately elevated ({ch_ratio:.1f}x)")
    print("  → Significant improvement, but optimization still needed")
else:
    print(f"\n  ✓ CH within reasonable range ({ch_ratio:.1f}x)")

print("\n" + "=" * 80)
print(" NEXT STEPS")
print("=" * 80)
print("\n1. This baseline establishes starting point with:")
print("   - Experimental electron density (3.3e9)")
print("   - All rates within literature bounds")
print("\n2. Next: Run constrained optimization")
print("   - Tune ~20-30 key rates within literature ranges")
print("   - Tune E field [10, 200] V/cm")
print("   - Minimize weighted error for H, CH, C2")
print("\n✓ Baseline complete, ready for optimization")
