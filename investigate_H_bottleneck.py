#!/usr/bin/env python3
"""
Investigate H bottleneck - why is H 42× lower than expected?
"""

import json

print("=" * 80)
print("H BOTTLENECK INVESTIGATION")
print("=" * 80)
print()

# Load best result
with open('optimization_results_charge_balanced/best_f27.0.json', 'r') as f:
    result = json.load(f)

H_actual = result['all_densities']['H']
H_target = 2.52e14
Te = result['Te']
Ne = result['Ne']
E = result['E_field']

print(f"Best result (f=27.0):")
print(f"  Te: {Te:.3f} eV")
print(f"  Ne: {Ne:.2e} cm⁻³")
print(f"  E: {E:.1f} V/cm")
print(f"  H: {H_actual:.2e} cm⁻³ ({H_actual/H_target:.2%} of target)")
print()

print("=" * 80)
print("H BALANCE FROM ANALYSIS TOOL")
print("=" * 80)
print()

# From balance analysis output
H_chem_production = 2.070e15  # cm⁻³/s
H_wall_loss = 1.908e15  # cm⁻³/s
H_chem_consumption = 2.002e15  # cm⁻³/s (total including wall)

k_wall = H_wall_loss / H_actual  # s⁻¹
k_total = H_chem_consumption / H_actual  # s⁻¹

print(f"Chemical production: {H_chem_production:.2e} cm⁻³/s")
print(f"Wall loss:           {H_wall_loss:.2e} cm⁻³/s")
print(f"Total consumption:   {H_chem_consumption:.2e} cm⁻³/s")
print()
print(f"Derived rates:")
print(f"  k_wall:  {k_wall:.1f} s⁻¹")
print(f"  k_total: {k_total:.1f} s⁻¹")
print()

print("Balance analysis says: ✓ At steady state (P/C = 1.034)")
print("This is TRUE for chemical reactions only, but...")
print()

print("=" * 80)
print("ADDING THE DRIFT TERM")
print("=" * 80)
print()

# Drift term from code
H_drift = 7.74e16  # cm⁻³/s (current value in odefun_optimized.py)

print(f"H_drift_gain (from code): {H_drift:.2e} cm⁻³/s")
print()

# Total H production including drift
H_total_production = H_drift + H_chem_production
H_total_loss = H_chem_consumption

print(f"Total H production (drift + chem): {H_total_production:.2e} cm⁻³/s")
print(f"Total H loss:                      {H_total_loss:.2e} cm⁻³/s")
print(f"Imbalance:                         {H_total_production/H_total_loss:.1f}×")
print()

# Expected steady-state H
H_expected = H_total_production / k_total
print(f"Expected H at steady state:")
print(f"  H = (drift + chem) / k_total")
print(f"  H = {H_total_production:.2e} / {k_total:.1f}")
print(f"  H = {H_expected:.2e} cm⁻³")
print()

print(f"Actual H:   {H_actual:.2e} cm⁻³")
print(f"Expected H: {H_expected:.2e} cm⁻³")
print(f"Discrepancy: {H_expected/H_actual:.1f}× (actual is {H_expected/H_actual:.1f}× too low!)")
print()

print("=" * 80)
print("POSSIBLE EXPLANATIONS")
print("=" * 80)
print()

print("1. Drift term is NOT being applied")
print("   → Need to verify H_drift_gain is actually added to dH/dt")
print()

print("2. There are MASSIVE H sinks not captured in balance analysis")
print("   → Balance tool only tracks chemical reactions explicitly")
print("   → Could be missing implicit sinks?")
print()

print("3. ODE integration hasn't reached true steady state")
print("   → But timescale analysis shows 100s >> 3ms (31862×)")
print("   → And balance shows P/C ≈ 1 for chemical reactions")
print()

print("4. Drift term should be MUCH smaller")
print("   → Current: 7.74e16 cm⁻³/s")
print("   → To match actual H: ~{:.2e} cm⁻³/s".format(H_actual * k_total - H_chem_production))
print(f"   → That would be {(H_actual * k_total - H_chem_production)/H_drift:.3f}× current value")
print()

# Calculate what drift should be for current H
H_current_loss = H_actual * k_total
drift_needed_for_current_H = H_current_loss - H_chem_production

print("=" * 80)
print("VERIFICATION TEST")
print("=" * 80)
print()

print("If H is at steady state with current density:")
print(f"  Total H loss = {H_actual:.2e} × {k_total:.1f} = {H_current_loss:.2e} cm⁻³/s")
print(f"  Chemical production = {H_chem_production:.2e} cm⁻³/s")
print(f"  Required drift = {drift_needed_for_current_H:.2e} cm⁻³/s")
print()

if drift_needed_for_current_H < 0:
    print("⚠  Required drift is NEGATIVE!")
    print("   This means chemical production EXCEEDS total loss!")
    print("   The drift term is being IGNORED or CANCELLED!")
else:
    print(f"  Ratio to actual drift: {drift_needed_for_current_H/H_drift:.4f}")
    print()
    if drift_needed_for_current_H / H_drift < 0.1:
        print("⚠  Required drift is ~100× smaller than actual drift!")
        print("   This strongly suggests the drift term is NOT being applied correctly.")

print()
print("=" * 80)
print("NEXT STEPS")
print("=" * 80)
print()
print("1. Verify drift term is in ODE: Check dH/dt includes H_drift_gain")
print("2. Test with zero drift: Run with H_drift_gain=0 and see if H changes")
print("3. Monitor H evolution: Track H vs time during integration")
print("4. Check for implicit H sinks in the ODE formulation")
