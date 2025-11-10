#!/usr/bin/env python3
"""
REPRODUCE the optimization result EXACTLY
Use identical code path as optimization
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

# Load result
with open('optimization_results_charge_balanced/best_f27.0.json', 'r') as f:
    result = json.load(f)

print("Reproducing optimization result...")
print(f"Target H from file: {result['all_densities']['H']:.2e} cm⁻³")
print()

# Import exactly as optimization does
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from build_reactions import build_reactions

# Define pressure helper
def pressure_to_density(P_mTorr, T=400):
    k_B = 1.380649e-23  # J/K
    P_Pa = P_mTorr * 133.322 / 1000
    return P_Pa / (k_B * T) / 1e6  # cm⁻³

n_total = pressure_to_density(500.0)

# Species list EXACTLY as in optimization
species = ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
           'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
           'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
           'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
           'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star']

# Build params EXACTLY as optimization
params_base = {
    'E_field': 50,
    'L_discharge': 0.45,
    'ne': 2.3e9,
    'Te': 1.0,
    'species': species,
    'mobilities': {
        'ArPlus': 1.54e3, 'CH4Plus': 1.54e3, 'CH3Plus': 1.54e3,
        'CH5Plus': 1.54e3, 'ArHPlus': 1.54e3, 'H3Plus': 1.54e3,
        'CH2Plus': 1.54e3, 'C2H5Plus': 1.54e3, 'C2H4Plus': 1.54e3,
        'C2H3Plus': 1.54e3, 'C2HPlus': 1.54e3, 'H2Plus': 1.54e3,
        'CHPlus': 1.54e3, 'CH3Minus': 1.54e3, 'HMinus': 1.54e3
    },
    'pressure': 500.0,
    'T': 400.0,
    'Tgas': 400.0
}

# Override with result params
params = params_base.copy()
params['E_field'] = result['E_field']
params['ne'] = result['Ne']
params['Te'] = result['Te']

# Build rates
k = define_rates(params)

# Override with optimized rates
for name, val in result['rate_values'].items():
    if name in k:
        k[name] = val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Create ODE
ode_func = PlasmaODE_Optimized(params)

print(f"ODE H_drift_gain: {ode_func.H_drift_gain:.2e} cm⁻³/s")
print()

# Initial conditions EXACTLY as optimization
ns = len(species)
y0 = np.ones(ns) * 1e3

def set_density(name, value):
    try:
        idx = species.index(name)
        y0[idx] = value
    except ValueError:
        pass

n_Ar = 0.85 * n_total
n_CH4 = 0.15 * n_total

set_density('e', params['ne'])
set_density('Ar', n_Ar)
set_density('CH4', n_CH4)
set_density('ArPlus', 1e7)
set_density('CH4Plus', 1e5)
set_density('CH3Plus', 1e5)
set_density('CH5Plus', 1e3)
set_density('ArHPlus', 5e5)
set_density('CH3Minus', 5e4)
set_density('H2', 1e12)
set_density('ArStar', 5e6)
set_density('H', 1e11)  # Initial H
set_density('C2', 5e7)
set_density('CH', 5e4)
set_density('C2H4', 5e7)
set_density('C2H6', 1e6)
set_density('CH2', 1e11)
set_density('C2H2', 1e12)
set_density('C2H5', 1e6)
set_density('CH3', 5e7)
set_density('C', 5e7)

H_idx = species.index('H')
print(f"Initial H: {y0[H_idx]:.2e} cm⁻³")
print()

# Integrate EXACTLY as optimization
print("Integrating (0, 100) seconds...")
sol = solve_ivp(
    ode_func,
    (0, 100),
    y0,
    method='BDF',
    rtol=1e-5,
    atol=1e-6,
    max_step=10.0
)

if not sol.success:
    print(f"✗ Integration failed: {sol.message}")
    exit(1)

y_final = sol.y[:, -1]
H_final = y_final[H_idx]

print(f"✓ Integration successful")
print(f"Final time: {sol.t[-1]:.2f}s")
print(f"Number of timesteps: {len(sol.t)}")
print()

print("=" * 80)
print("RESULT COMPARISON")
print("=" * 80)
print()

H_from_file = result['all_densities']['H']
print(f"H from result file: {H_from_file:.2e} cm⁻³")
print(f"H from integration: {H_final:.2e} cm⁻³")
print(f"Ratio: {H_final / H_from_file:.1f}×")
print()

if abs(H_final - H_from_file) / H_from_file < 0.01:
    print("✓✓✓ MATCH! Results are consistent.")
else:
    print("✗✗✗ MISMATCH!")
    print()
    print("The saved result does NOT match what the integration produces!")
    print("This means:")
    print("  1. The result file may be corrupted")
    print("  2. Or there's a bug in how results are saved")
    print("  3. Or the drift term wasn't active when this result was saved")
