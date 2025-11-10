#!/usr/bin/env python3
"""
Check characteristic timescales to determine required integration time
"""

import json

with open('optimization_results_charge_balanced/best_f61.8.json', 'r') as f:
    result = json.load(f)

print("=" * 80)
print("CHARACTERISTIC TIMESCALES")
print("=" * 80)
print()

# H timescale
H_density = result['all_densities']['H']
k_wall_H = result['rate_values'].get('stick_H_9_1', 0)
tau_H = 1 / k_wall_H if k_wall_H > 0 else float('inf')

print(f"H atom:")
print(f"  Density: {H_density:.2e} cm⁻³")
print(f"  Loss rate: {k_wall_H:.1f} s⁻¹")
print(f"  Lifetime: {tau_H*1000:.2f} ms")
print(f"  5× lifetime: {5*tau_H:.2f} s")
print()

# CH timescale
CH_density = result['all_densities']['CH']
k_wall_CH = result['rate_values'].get('stick_CH_9_3', 0)
k_loss_CH = result['rate_values'].get('loss_CH_11_9', 0)
k_total_CH = k_wall_CH + k_loss_CH
tau_CH = 1 / k_total_CH if k_total_CH > 0 else float('inf')

print(f"CH radical:")
print(f"  Density: {CH_density:.2e} cm⁻³")
print(f"  Loss rate: {k_total_CH:.1f} s⁻¹ (wall: {k_wall_CH:.1f}, ambipolar: {k_loss_CH:.1f})")
print(f"  Lifetime: {tau_CH*1000:.2f} ms")
print(f"  5× lifetime: {5*tau_CH:.2f} s")
print()

# C2 timescale
C2_density = result['all_densities']['C2']
k_wall_C2 = result['rate_values'].get('stick_C2_9_9', 0)
k_loss_C2 = result['rate_values'].get('loss_C2_11_3', 0)
k_total_C2 = k_wall_C2 + k_loss_C2
tau_C2 = 1 / k_total_C2 if k_total_C2 > 0 else float('inf')

print(f"C2:")
print(f"  Density: {C2_density:.2e} cm⁻³")
print(f"  Loss rate: {k_total_C2:.1f} s⁻¹ (wall: {k_wall_C2:.1f}, ambipolar: {k_loss_C2:.1f})")
print(f"  Lifetime: {tau_C2*1000:.2f} ms")
print(f"  5× lifetime: {5*tau_C2:.2f} s")
print()

# C2H6 timescale (slowest hydrocarbon)
C2H6_density = result['all_densities'].get('C2H6', 0)
k_wall_C2H6 = result['rate_values'].get('stick_C2H6_9_14', 0)
tau_C2H6 = 1 / k_wall_C2H6 if k_wall_C2H6 > 0 else float('inf')

print(f"C2H6 (slowest stable):")
print(f"  Density: {C2H6_density:.2e} cm⁻³")
print(f"  Loss rate: {k_wall_C2H6:.1f} s⁻¹")
print(f"  Lifetime: {tau_C2H6*1000:.2f} ms")
print(f"  5× lifetime: {5*tau_C2H6:.2f} s")
print()

# Longest timescale
max_tau = max(tau_H, tau_CH, tau_C2, tau_C2H6)
species_max = ['H', 'CH', 'C2', 'C2H6'][[tau_H, tau_CH, tau_C2, tau_C2H6].index(max_tau)]

print("=" * 80)
print("RECOMMENDATION")
print("=" * 80)
print()
print(f"Longest timescale: {max_tau*1000:.2f} ms ({species_max})")
print(f"5× longest: {5*max_tau:.2f} s")
print(f"10× longest: {10*max_tau:.2f} s")
print()
print(f"Current integration time: 100 s")
print(f"Ratio to longest timescale: {100/max_tau:.0f}× ")
print()

if 100 < 5*max_tau:
    print(f"⚠  Integration time is INSUFFICIENT!")
    print(f"   Recommended: {5*max_tau:.0f} s minimum (5× longest timescale)")
    print(f"   Better: {10*max_tau:.0f} s (10× longest timescale)")
elif 100 < 10*max_tau:
    print(f"Integration time is marginal.")
    print(f"   Current: {100:.0f} s")
    print(f"   Recommended: {10*max_tau:.0f} s (10× longest timescale)")
else:
    print(f"✓ Integration time is sufficient for steady state")

print()
print("Update in optimize_charge_balanced_500mTorr.py:")
print(f"  solve_ivp(..., (0, {int(10*max_tau)}), ...)")
