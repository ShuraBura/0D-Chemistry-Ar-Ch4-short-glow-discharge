#!/usr/bin/env python3
"""
Physical analysis of H wall sticking constraints
"""

import numpy as np

print("=" * 80)
print("H WALL STICKING - PHYSICAL CONSTRAINTS")
print("=" * 80)
print()

# Constants
k_B = 1.380649e-16  # erg/K
T_gas = 400  # K
m_H = 1.008 * 1.66054e-24  # g

# Calculate H thermal velocity
v_thermal = np.sqrt(8 * k_B * T_gas / (np.pi * m_H))  # cm/s
print(f"H thermal velocity at {T_gas}K: {v_thermal:.3e} cm/s ({v_thermal/1e5:.2f} km/s)")
print()

# Wall loss rate formula: k_wall = γ × v̄ × (A/V) / 4
# For cylinder with L >> R: A/V ≈ 2/R
# So: k_wall ≈ γ × v̄ / (2R)
# Therefore: γ/R = 2 × k_wall / v̄

print("=" * 80)
print("STICKING COEFFICIENT FROM WALL LOSS RATE")
print("=" * 80)
print()
print("For cylindrical discharge (L >> R):")
print("k_wall ≈ γ × v̄ / (2R)")
print("Therefore: γ/R = 2 × k_wall / v̄")
print()

# Current constraint in optimization
k_wall_min = 100  # s⁻¹
k_wall_nom = 300  # s⁻¹
k_wall_max = 500  # s⁻¹

# Current value from best optimization
k_wall_current = 289  # s⁻¹ (from balance analysis)

# Required to balance drift with target H
H_drift = 3.2e17  # cm⁻³/s (from odefun.py)
H_target = 2.4e14  # cm⁻³
k_wall_needed = H_drift / H_target  # s⁻¹

print(f"Current optimization constraints:")
print(f"  k_wall range: {k_wall_min}-{k_wall_max} s⁻¹")
print(f"  k_wall nominal: {k_wall_nom} s⁻¹")
print()

print(f"Current optimization result:")
print(f"  k_wall = {k_wall_current:.1f} s⁻¹")
print(f"  γ/R = {2*k_wall_current/v_thermal:.6f} cm⁻¹")
print()

print(f"Required to balance drift (3.2e17 cm⁻³/s) at target H (2.4e14 cm⁻³):")
print(f"  k_wall needed = {k_wall_needed:.1f} s⁻¹")
print(f"  γ/R = {2*k_wall_needed/v_thermal:.6f} cm⁻¹")
print()

# Check physical plausibility for different tube radii
print("=" * 80)
print("STICKING COEFFICIENTS FOR DIFFERENT DISCHARGE GEOMETRIES")
print("=" * 80)
print()

radii = [0.5, 1.0, 1.5, 2.0, 2.5]  # cm

print("Typical H sticking coefficients on surfaces:")
print("  - Clean metals: γ ~ 0.1-1.0")
print("  - Oxides: γ ~ 0.001-0.01")
print("  - Hydrocarbons/polymers: γ ~ 0.0001-0.006")
print("  → In CH4 plasma, walls coated with hydrocarbons: γ ~ 0.001-0.006")
print()

print(f"{'R (cm)':<10} {'Current k_wall':<20} {'γ (current)':<15} {'Needed k_wall':<20} {'γ (needed)':<15} {'Status':<20}")
print("-" * 100)

for R in radii:
    gamma_current = 2 * k_wall_current / v_thermal * R
    gamma_needed = 2 * k_wall_needed / v_thermal * R

    # Check if physically realistic (hydrocarbon surface)
    current_ok = 0.0001 <= gamma_current <= 0.006
    needed_ok = 0.0001 <= gamma_needed <= 0.006

    if current_ok and needed_ok:
        status = "Both OK"
    elif current_ok:
        status = "Needed too high"
    elif needed_ok:
        status = "Current marginal"
    else:
        status = "Both unrealistic"

    print(f"{R:<10.1f} {k_wall_current:<20.1f} {gamma_current:<15.6f} {k_wall_needed:<20.1f} {gamma_needed:<15.6f} {status:<20}")

print()
print("=" * 80)
print("CONCLUSION")
print("=" * 80)
print()

# Calculate what drift would be consistent with physical sticking
gamma_max = 0.006  # Maximum realistic γ for hydrocarbons
R_typical = 1.5  # cm (typical tube radius)
k_wall_physical_max = gamma_max * v_thermal / (2 * R_typical)
H_drift_physical = k_wall_physical_max * H_target

print(f"Maximum physically realistic H sticking for hydrocarbon surface:")
print(f"  γ_max = {gamma_max:.4f}")
print(f"  R = {R_typical:.1f} cm")
print(f"  k_wall_max = {k_wall_physical_max:.1f} s⁻¹")
print()

print(f"At target H = {H_target:.2e} cm⁻³:")
print(f"  Maximum H drift that can be balanced: {H_drift_physical:.2e} cm⁻³/s")
print()

print(f"Current H drift: {H_drift:.2e} cm⁻³/s")
print(f"Ratio: {H_drift/H_drift_physical:.1f}× too large")
print()

print("RECOMMENDATIONS:")
print(f"  1. Reduce H_drift_gain from {H_drift:.2e} to ~{H_drift_physical:.2e} cm⁻³/s")
print(f"  2. OR accept lower H target: ~{H_drift/k_wall_physical_max:.2e} cm⁻³")
print(f"  3. OR constrain k_wall to physically realistic max: ~{k_wall_physical_max:.0f} s⁻¹")
