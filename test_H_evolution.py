#!/usr/bin/env python3
"""
Test H evolution during ODE integration to find bottleneck
"""

import numpy as np
from scipy.integrate import solve_ivp
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Load best result
with open('optimization_results_charge_balanced/best_f61.8.json', 'r') as f:
    result = json.load(f)

# Import necessary modules
from odefun_optimized import PlasmaODE_Optimized
from define_species import species
from build_reactions import build_reactions

print("=" * 80)
print("H EVOLUTION TEST")
print("=" * 80)
print()

# Build params from result
params = {
    'species': species,
    'k': result['rate_values']
}
params['R'], params['tags'] = build_reactions(params)

# Create ODE
ode = PlasmaODE_Optimized(params)

print(f"H_drift_gain in ODE: {ode.H_drift_gain:.2e} cm⁻³/s")
print()

# Initial conditions (from result)
y0 = np.array([result['all_densities'][sp] for sp in species])
H_idx = species.index('H')

print(f"Initial H: {y0[H_idx]:.2e} cm⁻³")
print()

# Integrate and monitor H
t_span = (0, 100)
t_eval = np.logspace(-6, np.log10(100), 200)  # Log-spaced time points

print("Integrating ODE...")
sol = solve_ivp(
    ode,
    t_span,
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

print(f"✓ Integration successful ({len(sol.t)} time points)")
print()

# Extract H evolution
H_evolution = sol.y[H_idx, :]
t = sol.t

# Find when H reaches 90% of final value
H_final = H_evolution[-1]
H_90pct = 0.9 * H_final
idx_90 = np.argmax(H_evolution >= H_90pct)
t_90 = t[idx_90] if idx_90 > 0 else t[-1]

print("H EVOLUTION SUMMARY:")
print(f"  Initial H: {H_evolution[0]:.2e} cm⁻³")
print(f"  Final H:   {H_final:.2e} cm⁻³")
print(f"  Change:    {H_final/H_evolution[0]:.2f}×")
print(f"  Time to 90% of final: {t_90:.2e} s")
print()

# Check if H is still changing significantly at end
dH_dt_final = (H_evolution[-1] - H_evolution[-10]) / (t[-1] - t[-10])
print(f"dH/dt at end: {dH_dt_final:.2e} cm⁻³/s")
print(f"Relative rate: {dH_dt_final/H_final:.2e} s⁻¹")
print()

if abs(dH_dt_final/H_final) > 1e-6:
    print("⚠  H is still changing! NOT at steady state!")
else:
    print("✓ H appears to be at steady state")

print()

# Calculate expected steady-state H
k_wall = result['rate_values'].get('stick_H_9_1', 0)
H_drift = ode.H_drift_gain
H_chem = 3.01e15  # From balance analysis

H_expected = (H_drift + H_chem) / k_wall
print(f"Expected H from balance:")
print(f"  Drift source: {H_drift:.2e} cm⁻³/s")
print(f"  Chem source:  {H_chem:.2e} cm⁻³/s")
print(f"  Wall sink (k={k_wall:.0f} s⁻¹)")
print(f"  Expected H:   {H_expected:.2e} cm⁻³")
print(f"  Actual H:     {H_final:.2e} cm⁻³")
print(f"  Ratio:        {H_final/H_expected:.3f}")
print()

# Plot H evolution
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.semilogx(t, H_evolution, 'b-', linewidth=2)
plt.axhline(H_expected, color='r', linestyle='--', label=f'Expected: {H_expected:.2e}')
plt.axhline(H_final, color='g', linestyle='--', label=f'Actual: {H_final:.2e}')
plt.xlabel('Time (s)')
plt.ylabel('H density (cm⁻³)')
plt.title('H atom evolution during ODE integration')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(2, 1, 2)
plt.loglog(t, np.abs(H_evolution - H_final), 'b-', linewidth=2)
plt.xlabel('Time (s)')
plt.ylabel('|H - H_final| (cm⁻³)')
plt.title('Convergence to final H')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('H_evolution_test.png', dpi=150)
print("Saved plot to: H_evolution_test.png")
print()

# Check drift term at final state
y_final = sol.y[:, -1]
dydt_final = ode(sol.t[-1], y_final)
print(f"At final state (t={sol.t[-1]:.1f} s):")
print(f"  dH/dt = {dydt_final[H_idx]:.2e} cm⁻³/s")
print(f"  H = {y_final[H_idx]:.2e} cm⁻³")
print()

# Manually calculate H balance at final state
wall_loss = k_wall * y_final[H_idx]
print(f"Manual H balance at final state:")
print(f"  Drift source:  {H_drift:.2e} cm⁻³/s")
print(f"  Chem source:   ~{H_chem:.2e} cm⁻³/s")
print(f"  Wall loss:     {wall_loss:.2e} cm⁻³/s")
print(f"  Net (drift+chem-wall): {H_drift + H_chem - wall_loss:.2e} cm⁻³/s")
print()

if abs(dydt_final[H_idx]) < 1e12:
    print("✓ dH/dt ≈ 0: System is at steady state")
    print()
    print("CONCLUSION: ODE integration is working correctly.")
    print("The 'low' H density is the actual steady state!")
    print()
    print("This means:")
    print("  1. Drift term IS being applied")
    print("  2. There are additional H sinks beyond wall loss")
    print("  3. OR chemical H production is lower than estimated")
else:
    print("✗ dH/dt ≠ 0: System NOT at steady state")
    print("Need to investigate why convergence failed")
