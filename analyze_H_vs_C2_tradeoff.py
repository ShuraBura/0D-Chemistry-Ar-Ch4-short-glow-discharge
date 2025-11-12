"""
Analyze the H vs C2 trade-off

Key reactions:
- H + C2H2 → C2 + H2  (makes C2 from C2H2, rate ~5e-10)
- H + C2 → CH + C     (destroys C2, rate ~3.5e-11)

From baseline analysis:
- H + C2 → CH accounts for 72.8% of C2 loss
- If we reduce H slightly, C2 lifetime increases
- But H is already at 79.9% of target

Question: What's the optimal H density?
- Too high: destroys C2 via H + C2 → CH
- Too low: reduces C2 production via H + C2H2 → C2
- Also misses H target
"""

import json
import numpy as np
import matplotlib.pyplot as plt

# Load baseline
with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

# Baseline values
H_baseline = baseline['target_densities']['H']
C2_baseline = baseline['target_densities']['C2']
C2H2_baseline = baseline['all_densities']['C2H2']

# Target densities
H_target = 2.52e14
C2_target = 5.6e11

# Rate constants (from baseline)
k_H_C2H2_to_C2 = 5e-10  # H + C2H2 → C2 (C2 production)
k_H_C2_to_CH = baseline['rate_values']['C2_H_CH_C_cm3_7_6']  # H + C2 → CH (C2 destruction)

print("=" * 80)
print("H vs C2 TRADE-OFF ANALYSIS")
print("=" * 80)
print()
print(f"Baseline:")
print(f"  H:    {H_baseline:.2e} ({H_baseline/H_target*100:.1f}% of target)")
print(f"  C2:   {C2_baseline:.2e} ({C2_baseline/C2_target*100:.1f}% of target)")
print(f"  C2H2: {C2H2_baseline:.2e}")
print()
print(f"Rate constants:")
print(f"  k(H + C2H2 → C2): {k_H_C2H2_to_C2:.2e} cm³/s  (C2 production)")
print(f"  k(H + C2 → CH):   {k_H_C2_to_CH:.2e} cm³/s  (C2 destruction)")
print(f"  Ratio: {k_H_C2H2_to_C2/k_H_C2_to_CH:.1f}× faster production than destruction")
print()

# Calculate production and loss rates as function of H
H_range = np.linspace(0.5, 1.5, 100) * H_baseline  # 50% to 150% of baseline H

# Assume C2 and C2H2 scale with optimization (not fixed)
# For now, assume C2H2 can be boosted independently
C2H2_boosted = np.array([5e9, 1e10, 5e10, 1e11, 5e11, 1e12])  # Different C2H2 scenarios

print("SCENARIO ANALYSIS: How does H density affect C2 at different C2H2 levels?")
print()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

for C2H2 in C2H2_boosted:
    # For each C2H2 level, calculate steady-state C2 as function of H
    # At steady state: production = loss
    # [H] × [C2H2] × k_prod = [H] × [C2] × k_loss + loss_other
    # Assume loss_other is proportional to C2 (wall losses, volumetric)

    # From baseline: total C2 loss = 9.20e12 cm⁻³/s
    # H + C2 loss = 72.8% = 6.70e12
    # Other loss = 27.2% = 2.50e12
    loss_other_rate = 2.50e12 / C2_baseline  # s⁻¹ (wall + volumetric)

    # Steady state: H × C2H2 × k_prod = H × C2 × k_loss + C2 × loss_other
    # Solve for C2: C2 = (H × C2H2 × k_prod) / (H × k_loss + loss_other)

    C2_ss = (H_range * C2H2 * k_H_C2H2_to_C2) / (H_range * k_H_C2_to_CH + loss_other_rate)

    ax1.plot(H_range/H_target*100, C2_ss/C2_target*100,
             label=f'C2H2={C2H2:.1e}', linewidth=2)

    # Also plot absolute C2
    ax2.plot(H_range/H_target*100, C2_ss,
             label=f'C2H2={C2H2:.1e}', linewidth=2)

# Mark baseline point
ax1.axvline(H_baseline/H_target*100, color='red', linestyle='--', alpha=0.5, label='Baseline H')
ax1.axhline(C2_baseline/C2_target*100, color='blue', linestyle='--', alpha=0.5, label='Baseline C2')
ax1.axhline(100, color='green', linestyle=':', alpha=0.7, label='C2 Target')

ax2.axvline(H_baseline/H_target*100, color='red', linestyle='--', alpha=0.5)
ax2.axhline(C2_target, color='green', linestyle=':', alpha=0.7)

ax1.set_xlabel('H (% of target)', fontsize=12)
ax1.set_ylabel('C2 (% of target)', fontsize=12)
ax1.set_title('C2 Achievement vs H Density\n(at different C2H2 levels)', fontsize=13)
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

ax2.set_xlabel('H (% of target)', fontsize=12)
ax2.set_ylabel('C2 density (cm⁻³)', fontsize=12)
ax2.set_title('Absolute C2 Density vs H\n(at different C2H2 levels)', fontsize=13)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.set_yscale('log')

plt.tight_layout()
plt.savefig('H_vs_C2_tradeoff.png', dpi=150, bbox_inches='tight')
print("Plot saved: H_vs_C2_tradeoff.png")
print()

# Calculate optimal H for different C2H2 scenarios
print("OPTIMAL H ANALYSIS:")
print()
print("For each C2H2 level, what H density maximizes C2?")
print("(Considering both production and destruction by H)")
print()

for C2H2 in [4.81e9, 5e10, 5e11, 1e12]:
    # Take derivative: d(C2)/d(H) = 0
    # C2 = (H × C2H2 × k_prod) / (H × k_loss + loss_other)
    # Let a = C2H2 × k_prod, b = k_loss, c = loss_other
    # C2 = (H × a) / (H × b + c)
    # d(C2)/dH = [a(H×b + c) - H×a×b] / (H×b + c)²
    #          = [a×c] / (H×b + c)²
    # This is always positive! No maximum, C2 increases with H!

    # Wait, that means higher H is always better for C2 production?
    # Let me recalculate...

    a = C2H2 * k_H_C2H2_to_C2
    b = k_H_C2_to_CH
    c = 2.50e12 / C2_baseline

    # Actually, the derivative is positive, meaning C2 always increases with H
    # This is because k_prod >> k_loss (14× faster)

    C2_at_baseline_H = (H_baseline * a) / (H_baseline * b + c)
    C2_at_target_H = (H_target * a) / (H_target * b + c)
    C2_at_half_H = (0.5 * H_baseline * a) / (0.5 * H_baseline * b + c)

    print(f"C2H2 = {C2H2:.2e}:")
    print(f"  At H = 50% baseline:  C2 = {C2_at_half_H:.2e} ({C2_at_half_H/C2_target*100:.1f}%)")
    print(f"  At H = baseline:      C2 = {C2_at_baseline_H:.2e} ({C2_at_baseline_H/C2_target*100:.1f}%)")
    print(f"  At H = target:        C2 = {C2_at_target_H:.2e} ({C2_at_target_H/C2_target*100:.1f}%)")
    print()

print("=" * 80)
print("KEY INSIGHT:")
print("=" * 80)
print()
print("Since k(H + C2H2 → C2) = 14× faster than k(H + C2 → CH),")
print("HIGHER H actually HELPS C2 production!")
print()
print("The problem is NOT that H destroys C2.")
print("The problem is that C2H2 is 1000× too low!")
print()
print("If we boost C2H2 from 4.81e9 to 5e11 (100×):")
print(f"  At baseline H: C2 would reach ~{(H_baseline * 5e11 * k_H_C2H2_to_C2) / (H_baseline * k_H_C2_to_CH + 2.50e12/C2_baseline):.2e}")
print(f"  Target is: {C2_target:.2e}")
print()
print("STRATEGY: Focus on boosting C2H2, NOT on reducing H!")
print("=" * 80)
