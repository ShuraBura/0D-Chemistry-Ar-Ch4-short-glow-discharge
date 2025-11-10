#!/usr/bin/env python3
"""
Integrate for full 100s and see where H ends up
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("=" * 80)
print("LONG INTEGRATION TEST")
print("=" * 80)
print()

# Load result
with open('optimization_results_charge_balanced/best_f27.0.json', 'r') as f:
    result = json.load(f)

# Setup
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from build_reactions import build_reactions

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
    'CHPlus': 1.54e3,
    'CH3Minus': 1.54e3, 'HMinus': 1.54e3
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

print("Building ODE...")
k = define_rates(params)
for name, val in result['rate_values'].items():
    if name in k:
        k[name] = val

params['k'] = k
params['R'], params['tags'] = build_reactions(params)
ode = PlasmaODE_Optimized(params)

print(f"✓ ODE created (H_drift_gain = {ode.H_drift_gain:.2e})")
print()

# Build state vector
y0 = np.array([result['all_densities'][sp] for sp in species_list])
H_idx = species_list.index('H')

print(f"Initial H: {y0[H_idx]:.2e} cm⁻³")
print()

# Integrate for 100s like optimization does
print("Integrating for 100 seconds (as optimization does)...")
t_eval = np.logspace(-6, 2, 200)  # 1 µs to 100s, log-spaced

sol = solve_ivp(
    ode,
    (0, 100),
    y0,
    method='BDF',
    t_eval=t_eval,
    rtol=1e-5,
    atol=1e-6,
    max_step=10.0
)

if not sol.success:
    print(f"✗ Integration failed: {sol.message}")
    exit(1)

print(f"✓ Integration successful")
print()

# Extract H evolution
H_evolution = sol.y[H_idx, :]
t = sol.t

print(f"H evolution:")
print(f"  t=0:      {H_evolution[0]:.2e} cm⁻³")
print(f"  t=1ms:    {H_evolution[np.argmin(np.abs(t - 0.001))]:.2e} cm⁻³")
print(f"  t=10ms:   {H_evolution[np.argmin(np.abs(t - 0.01))]:.2e} cm⁻³")
print(f"  t=100ms:  {H_evolution[np.argmin(np.abs(t - 0.1))]:.2e} cm⁻³")
print(f"  t=1s:     {H_evolution[np.argmin(np.abs(t - 1.0))]:.2e} cm⁻³")
print(f"  t=10s:    {H_evolution[np.argmin(np.abs(t - 10.0))]:.2e} cm⁻³")
print(f"  t=100s:   {H_evolution[-1]:.2e} cm⁻³")
print()

print(f"Final / Initial: {H_evolution[-1] / H_evolution[0]:.1f}×")
print()

# Check if at steady state
dydt_final = ode(sol.t[-1], sol.y[:, -1])
print(f"dH/dt at t=100s: {dydt_final[H_idx]:.2e} cm⁻³/s")

if abs(dydt_final[H_idx]) < 1e12:
    print("✓ Reached steady state")
else:
    print(f"✗ Still changing (relative rate: {dydt_final[H_idx] / sol.y[H_idx, -1]:.2e} s⁻¹)")

print()

# Plot
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.semilogx(t, H_evolution, 'b-', linewidth=2)
plt.axhline(2.52e14, color='r', linestyle='--', label='Target: 2.52e14')
plt.xlabel('Time (s)')
plt.ylabel('H density (cm⁻³)')
plt.title('H Evolution (full time)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(2, 2, 2)
plt.plot(t, H_evolution, 'b-', linewidth=2)
plt.axhline(2.52e14, color='r', linestyle='--', label='Target')
plt.xlabel('Time (s)')
plt.ylabel('H density (cm⁻³)')
plt.title('H Evolution (linear time)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(2, 2, 3)
plt.loglog(t, np.abs(H_evolution - H_evolution[-1]) + 1e8, 'b-', linewidth=2)
plt.xlabel('Time (s)')
plt.ylabel('|H - H_final| (cm⁻³)')
plt.title('Convergence to final H')
plt.grid(True, alpha=0.3)

plt.subplot(2, 2, 4)
# Plot dH/dt vs time
dH_dt = np.diff(H_evolution) / np.diff(t)
t_mid = (t[:-1] + t[1:]) / 2
plt.semilogx(t_mid, dH_dt, 'b-', linewidth=2)
plt.axhline(0, color='r', linestyle='--')
plt.xlabel('Time (s)')
plt.ylabel('dH/dt (cm⁻³/s)')
plt.title('H Growth Rate')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('H_integration_100s.png', dpi=150)
print("Saved plot: H_integration_100s.png")
print()

print("=" * 80)
print("CONCLUSION")
print("=" * 80)
print()

if H_evolution[-1] > 10 * H_evolution[0]:
    print("H is GROWING RAPIDLY throughout the 100s integration!")
    print("The optimization is NOT seeing the steady state.")
    print()
    print("Possible reasons:")
    print("  1. Integration terminates before 100s (solver issues)")
    print("  2. max_step=10s prevents fine resolution")
    print("  3. Solver declares 'success' prematurely")
elif H_evolution[-1] < 0.1 * H_evolution[0]:
    print("H is DECAYING!")
else:
    print(f"H changed by {H_evolution[-1] / H_evolution[0]:.1f}× over 100s")
