#!/usr/bin/env python3
"""
Performance benchmark: Original vs Optimized ODE function
"""

import numpy as np
import time

from define_rates import define_rates
from build_reactions import build_reactions
from odefun import PlasmaODE
from odefun_optimized import PlasmaODE_Optimized

print("=" * 70)
print(" Performance Benchmark: Original vs Optimized")
print("=" * 70)

# Setup
params = {
    'E_field': 50, 'L_discharge': 0.45, 'ne': 1e10,
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

params['k'] = define_rates(params)
params['R'], params['tags'] = build_reactions(params)

# Create test state
y_test = np.ones(len(params['species'])) * 1e10

print(f"\nSystem size:")
print(f"  Species: {len(params['species'])}")
print(f"  Reactions: {len(params['R'])}")
print(f"\nBenchmark: 1000 ODE function calls")

# Original version
print("\n--- Original (loop-based) ---")
ode_orig = PlasmaODE(params)
start = time.time()
for _ in range(1000):
    result = ode_orig(0, y_test)
time_orig = time.time() - start
print(f"  Time: {time_orig:.3f} s")
print(f"  Per call: {time_orig/1000*1000:.3f} ms")

# Optimized version
print("\n--- Optimized (vectorized) ---")
ode_opt = PlasmaODE_Optimized(params)
start = time.time()
for _ in range(1000):
    result = ode_opt(0, y_test)
time_opt = time.time() - start
print(f"  Time: {time_opt:.3f} s")
print(f"  Per call: {time_opt/1000*1000:.3f} ms")

# Comparison
print("\n" + "=" * 70)
print(" Results")
print("=" * 70)
print(f"  Speedup: {time_orig/time_opt:.1f}x faster")
print(f"  Time saved per call: {(time_orig-time_opt)/1000*1000:.3f} ms")
print("=" * 70)
