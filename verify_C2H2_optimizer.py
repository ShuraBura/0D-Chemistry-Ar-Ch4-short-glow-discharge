#!/usr/bin/env python3
"""
Check what happens at true steady state for C2H2-tunable-loss optimizer result
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

# Load C2H2 optimizer result
with open('optimization_results_C2H2_tunable_loss/best_f248.8.json', 'r') as f:
    best = json.load(f)

print("=" * 80)
print("C2H2-TUNABLE-LOSS OPTIMIZER - STEADY STATE VERIFICATION")
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
print(f"  Ni/Ne = {best.get('Ni_over_Ne', 0):.2f}")
print()

print(f"Optimizer reported (at t=500s):")
print(f"  H    = {best['target_densities']['H']:.2e} cm⁻³ ({best['target_densities']['H']/2.52e14*100:.1f}% of target)")
print(f"  CH   = {best['target_densities']['CH']:.2e} cm⁻³ ({best['target_densities']['CH']/1.0e9*100:.1f}% of target)")
print(f"  C2   = {best['target_densities']['C2']:.2e} cm⁻³ ({best['target_densities']['C2']/5.6e11*100:.1f}% of target)")
print(f"  C2H2 = {best['all_densities']['C2H2']:.2e} cm⁻³ ({best['all_densities']['C2H2']/5e12*100:.1f}% of 5e+12 needed)")
print()

print(f"C2H2 loss rates:")
rv = best['rate_values']
print(f"  stick_C2H2: {rv.get('stick_C2H2_9_11', 'N/A'):.2e} s⁻¹")
print(f"  loss_C2H2:  {rv.get('loss_C2H2_11_19', 'N/A'):.2e} s⁻¹")
print()

# Run for MUCH longer to verify steady state
print("Running extended integration to t=2000s...")
ode_func = PlasmaODE_Optimized(params)

sol = solve_ivp(
    ode_func,
    (0, 2000),
    y0,
    method='BDF',
    rtol=1e-7,
    atol=1e-9,
    max_step=1.0
)

print(f"Integration completed: {len(sol.t)} steps to t={sol.t[-1]:.0f}s")
print()

idx_H = species.index('H')
idx_CH = species.index('CH')
idx_C2 = species.index('C2')
idx_C2H2 = species.index('C2H2')

y_final = sol.y[:, -1]
H_final = y_final[idx_H]
CH_final = y_final[idx_CH]
C2_final = y_final[idx_C2]
C2H2_final = y_final[idx_C2H2]

print(f"True steady state (at t={sol.t[-1]:.0f}s):")
print(f"  H    = {H_final:.2e} cm⁻³ (factor from t=500s: {H_final/best['target_densities']['H']:.1f}×)")
print(f"  CH   = {CH_final:.2e} cm⁻³ (factor from t=500s: {CH_final/best['target_densities']['CH']:.1f}×)")
print(f"  C2   = {C2_final:.2e} cm⁻³ (factor from t=500s: {C2_final/best['target_densities']['C2']:.1f}×)")
print(f"  C2H2 = {C2H2_final:.2e} cm⁻³ (factor from t=500s: {C2H2_final/best['all_densities']['C2H2']:.1f}×)")
print()

print(f"Achievement vs Target:")
print(f"  H:    {H_final/2.52e14*100:6.1f}% of target")
print(f"  CH:   {CH_final/1.0e9*100:6.1f}% of target")
print(f"  C2:   {C2_final/5.6e11*100:6.1f}% of target")
print(f"  C2H2: {C2H2_final/5e12*100:6.1f}% of 5e+12 needed")
print()

# Check rates
dydt_final = ode_func(sol.t[-1], y_final)
print(f"Steady state verification (dn/dt at t={sol.t[-1]:.0f}s):")
print(f"  dH/dt    = {dydt_final[idx_H]:.2e} (rel: {abs(dydt_final[idx_H]/H_final):.2e})")
print(f"  dCH/dt   = {dydt_final[idx_CH]:.2e} (rel: {abs(dydt_final[idx_CH]/CH_final):.2e})")
print(f"  dC2/dt   = {dydt_final[idx_C2]:.2e} (rel: {abs(dydt_final[idx_C2]/C2_final):.2e})")
print(f"  dC2H2/dt = {dydt_final[idx_C2H2]:.2e} (rel: {abs(dydt_final[idx_C2H2]/C2H2_final):.2e})")
print()
