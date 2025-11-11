#!/usr/bin/env python3
"""
Check if best result has reached true steady state
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

def pressure_to_density(pressure_mTorr, T_K=400):
    kB = 1.38064852e-23
    Torr_to_Pa = 133.322
    P_Pa = pressure_mTorr * 1e-3 * Torr_to_Pa
    n_m3 = P_Pa / (kB * T_K)
    return n_m3 * 1e-6

# Load best result
with open('optimization_results_charge_balanced_fixed_weights/best_f109.6.json', 'r') as f:
    best = json.load(f)

print("=" * 80)
print("STEADY STATE CHECK FOR BEST RESULT")
print("=" * 80)
print()

# Setup params
params = {
    'P': 500.0,
    'Te': best['Te'],
    'ne': best['Ne'],
    'E_field': best['E_field'],
    'L_discharge': 0.45,
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
    },
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar',
                'C2H4', 'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C',
                'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H',
                'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus',
                'H2Plus', 'C2H2Star'],
}

n_total = pressure_to_density(500.0)
k = define_rates(params)

# Apply tuned rates
for rate_name, rate_val in best.get('rate_values', {}).items():
    if rate_name in k:
        k[rate_name] = rate_val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Set initial conditions from best result
species = params['species']
y0 = np.array([best['all_densities'].get(s, 1e3) for s in species])

print(f"Parameters:")
print(f"  Te = {best['Te']:.2f} eV")
print(f"  Ne = {best['Ne']:.2e} cm⁻³")
print(f"  E  = {best['E_field']:.1f} V/cm")
print()

print(f"Optimizer reported (at t=100s):")
print(f"  H  = {best['target_densities']['H']:.2e} cm⁻³")
print(f"  CH = {best['target_densities']['CH']:.2e} cm⁻³")
print(f"  C2 = {best['target_densities']['C2']:.2e} cm⁻³")
print()

# Run for MUCH longer to ensure true steady state
print("Running extended integration to t=1000s...")
ode_func = PlasmaODE_Optimized(params)

sol = solve_ivp(
    ode_func,
    (0, 1000),  # Much longer integration
    y0,
    method='BDF',
    rtol=1e-7,
    atol=1e-9,
    max_step=0.5
)

print(f"Integration completed: {len(sol.t)} steps")
print()

idx_H = species.index('H')
idx_CH = species.index('CH')
idx_C2 = species.index('C2')

y_final = sol.y[:, -1]
H_final = y_final[idx_H]
CH_final = y_final[idx_CH]
C2_final = y_final[idx_C2]

print(f"True steady state (at t={sol.t[-1]:.0f}s):")
print(f"  H  = {H_final:.2e} cm⁻³ (factor: {H_final/best['target_densities']['H']:.1f}×)")
print(f"  CH = {CH_final:.2e} cm⁻³ (factor: {CH_final/best['target_densities']['CH']:.1f}×)")
print(f"  C2 = {C2_final:.2e} cm⁻³ (factor: {C2_final/best['target_densities']['C2']:.1f}×)")
print()

# Check rates
dydt_final = ode_func(sol.t[-1], y_final)
print(f"Steady state verification (dn/dt at t={sol.t[-1]:.0f}s):")
print(f"  dH/dt  = {dydt_final[idx_H]:.2e} (rel: {dydt_final[idx_H]/H_final:.2e})")
print(f"  dCH/dt = {dydt_final[idx_CH]:.2e} (rel: {dydt_final[idx_CH]/CH_final:.2e})")
print(f"  dC2/dt = {dydt_final[idx_C2]:.2e} (rel: {dydt_final[idx_C2]/C2_final:.2e})")
print()

if abs(dydt_final[idx_H]/H_final) < 1e-10:
    print("✓ H at true steady state")
else:
    print(f"⚠ H NOT at steady state (dH/dt/H = {dydt_final[idx_H]/H_final:.2e})")

if abs(dydt_final[idx_CH]/CH_final) < 1e-10:
    print("✓ CH at true steady state")
else:
    print(f"⚠ CH NOT at steady state (dCH/dt/CH = {dydt_final[idx_CH]/CH_final:.2e})")

if abs(dydt_final[idx_C2]/C2_final) < 1e-10:
    print("✓ C2 at true steady state")
else:
    print(f"⚠ C2 NOT at steady state (dC2/dt/C2 = {dydt_final[idx_C2]/C2_final:.2e})")
