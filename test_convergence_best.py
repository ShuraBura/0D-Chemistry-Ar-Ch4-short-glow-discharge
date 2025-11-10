#!/usr/bin/env python3
"""
Test if the best result reaches steady state when integrated
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from build_reactions import build_reactions

print("=" * 80)
print("CONVERGENCE TEST FOR BEST RESULT")
print("=" * 80)
print()

# Load result
with open('optimization_results_charge_balanced/best_f41.3.json', 'r') as f:
    result = json.load(f)

# Setup exactly as optimizer does
species_list = ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star']

mobilities = {
    'ArPlus': 1.54e3, 'CH4Plus': 1.54e3, 'CH3Plus': 1.54e3,
    'CH5Plus': 1.54e3, 'ArHPlus': 1.54e3, 'H3Plus': 1.54e3,
    'CH2Plus': 1.54e3, 'C2H5Plus': 1.54e3, 'C2H4Plus': 1.54e3,
    'C2H3Plus': 1.54e3, 'C2HPlus': 1.54e3, 'H2Plus': 1.54e3,
    'CHPlus': 1.54e3, 'CH3Minus': 1.54e3, 'HMinus': 1.54e3
}

params = {
    'E_field': result['E_field'],
    'L_discharge': 0.45,
    'ne': result['Ne'],
    'Te': result['Te'],
    'species': species_list,
    'T': 400.0,
    'Tgas': 400.0,
    'pressure': 500.0,
    'mobilities': mobilities
}

print("Building ODE with optimized rates...")
k = define_rates(params)
for name, val in result['rate_values'].items():
    if name in k:
        k[name] = val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)
ode = PlasmaODE_Optimized(params)

H_idx = species_list.index('H')
print(f"✓ ODE created (H_drift_gain = {ode.H_drift_gain:.2e})")
print()

# Build initial state from result
y0 = np.array([result['all_densities'][sp] for sp in species_list])

print(f"Starting from result state:")
print(f"  H = {y0[H_idx]:.2e} cm⁻³")
print()

# Check dH/dt at this state
dydt0 = ode(0, y0)
print(f"dH/dt at result state: {dydt0[H_idx]:.2e} cm⁻³/s")
print(f"Relative rate: {dydt0[H_idx]/y0[H_idx]:.2e} s⁻¹")
print()

if abs(dydt0[H_idx]) < 1e12:
    print("✓ Already at steady state")
else:
    print("✗ NOT at steady state!")
    print()
    print("Integrating to find TRUE steady state...")
    print()

    # Integrate with tighter tolerances
    sol = solve_ivp(
        ode,
        (0, 100),
        y0,
        method='BDF',
        t_eval=np.logspace(-6, 2, 200),
        rtol=1e-6,
        atol=1e-8,
        max_step=1.0  # Smaller max_step
    )

    if not sol.success:
        print(f"✗ Integration failed: {sol.message}")
    else:
        H_evolution = sol.y[H_idx, :]
        t = sol.t

        print(f"✓ Integration completed")
        print(f"  t_final = {sol.t[-1]:.2f} s")
        print(f"  nsteps = {len(sol.t)}")
        print()

        print("H evolution:")
        print(f"  t=0:      {H_evolution[0]:.2e} cm⁻³")
        print(f"  t=1ms:    {H_evolution[np.argmin(np.abs(t - 0.001))]:.2e} cm⁻³")
        print(f"  t=10ms:   {H_evolution[np.argmin(np.abs(t - 0.01))]:.2e} cm⁻³")
        print(f"  t=100ms:  {H_evolution[np.argmin(np.abs(t - 0.1))]:.2e} cm⁻³")
        print(f"  t=1s:     {H_evolution[np.argmin(np.abs(t - 1.0))]:.2e} cm⁻³")
        print(f"  t=10s:    {H_evolution[np.argmin(np.abs(t - 10.0))]:.2e} cm⁻³")
        print(f"  t=100s:   {H_evolution[-1]:.2e} cm⁻³")
        print()

        # Check final dH/dt
        dydt_final = ode(sol.t[-1], sol.y[:, -1])
        print(f"dH/dt at t=100s: {dydt_final[H_idx]:.2e} cm⁻³/s")

        if abs(dydt_final[H_idx]) < 1e12:
            print("✓ Reached steady state")
        else:
            print(f"✗ Still evolving (relative rate: {dydt_final[H_idx]/sol.y[H_idx,-1]:.2e} s⁻¹)")

        print()
        print("=" * 80)
        print("CONCLUSION")
        print("=" * 80)
        print()

        factor = H_evolution[-1] / H_evolution[0]
        print(f"H changed by {factor:.1f}× during integration")
        print(f"  Initial: {H_evolution[0]:.2e} cm⁻³")
        print(f"  Final:   {H_evolution[-1]:.2e} cm⁻³")
        print()

        if factor > 10:
            print("The optimizer result was FAR from steady state!")
            print("The integration in the optimizer did NOT converge.")
        else:
            print("H is relatively stable.")
