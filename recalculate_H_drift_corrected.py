#!/usr/bin/env python3
"""
Recalculate H drift for corrected target density
"""

print("=" * 80)
print("H DRIFT RECALCULATION - CORRECTED TARGET")
print("=" * 80)
print()

# Corrected experimental value
H_target_corrected = 2.52e14  # cm⁻³ (radially averaged at y=4mm)
H_target_old = 3.82e14  # cm⁻³ (was incorrect)

print(f"Experimental H density (y=4mm, averaged -6<r<6mm):")
print(f"  Corrected value: {H_target_corrected:.2e} cm⁻³")
print(f"  Previous value:  {H_target_old:.2e} cm⁻³")
print(f"  Ratio: {H_target_old/H_target_corrected:.2f}× (was overestimated)")
print()

# From best optimization result (f=61.8)
k_wall_current = 319  # s⁻¹ (wall sticking rate from balance analysis)
H_chem_production = 3.01e15  # cm⁻³/s (chemical H production, relatively constant)

print("From current best result (f=61.8):")
print(f"  H wall sticking rate: {k_wall_current:.0f} s⁻¹")
print(f"  Chemical H production: {H_chem_production:.2e} cm⁻³/s")
print()

# Calculate required drift for corrected target
wall_loss_at_target = k_wall_current * H_target_corrected
H_drift_required = wall_loss_at_target - H_chem_production

print("Required drift for corrected target:")
print(f"  Wall loss at H={H_target_corrected:.2e}: {wall_loss_at_target:.2e} cm⁻³/s")
print(f"  Chemical production: {H_chem_production:.2e} cm⁻³/s")
print(f"  Required drift: {H_drift_required:.2e} cm⁻³/s")
print()

# Compare to current drift values
H_drift_old = 3.2e17  # Original (too high)
H_drift_current = 1.06e17  # First correction
H_drift_new = H_drift_required

print("Drift value comparison:")
print(f"  Original (wrong):        {H_drift_old:.2e} cm⁻³/s")
print(f"  First correction:        {H_drift_current:.2e} cm⁻³/s ({H_drift_current/H_drift_old:.2f}× of original)")
print(f"  New (for 2.52e14 target): {H_drift_new:.2e} cm⁻³/s ({H_drift_new/H_drift_old:.2f}× of original)")
print(f"  Change from current:     {H_drift_new/H_drift_current:.2f}×")
print()

# Sanity check with current H density
H_current = 1.42e13  # cm⁻³
total_loss_current = k_wall_current * H_current
drift_plus_chem_current = H_drift_current + H_chem_production

print("Sanity check with current result:")
print(f"  Current H: {H_current:.2e} cm⁻³")
print(f"  Total production (drift+chem): {drift_plus_chem_current:.2e} cm⁻³/s")
print(f"  Total loss (wall): {total_loss_current:.2e} cm⁻³/s")
print(f"  Expected steady-state H: {drift_plus_chem_current/k_wall_current:.2e} cm⁻³")
print(f"  Actual H / Expected: {H_current/(drift_plus_chem_current/k_wall_current):.3f}")
print()

if H_current < 0.5 * (drift_plus_chem_current/k_wall_current):
    print("⚠  Current H is much lower than expected from steady-state balance!")
    print("   This suggests:")
    print("   1. ODE integration hasn't converged to steady state")
    print("   2. OR there are additional H sinks not captured")
    print("   3. OR drift term has implementation issues")
    print()

print("=" * 80)
print("RECOMMENDATION")
print("=" * 80)
print()
print(f"Update code:")
print(f"  1. H target: {H_target_corrected:.2e} cm⁻³ (was {H_target_old:.2e})")
print(f"  2. H_drift_gain: {H_drift_new:.2e} cm⁻³/s (currently {H_drift_current:.2e})")
print()
print(f"This should give H ≈ {H_target_corrected:.2e} cm⁻³ at steady state")
