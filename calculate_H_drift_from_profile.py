#!/usr/bin/env python3
"""
Calculate H drift source term from spatial profile data
"""

import numpy as np

print("=" * 80)
print("H DRIFT SOURCE FROM SPATIAL PROFILE")
print("=" * 80)
print()

# Experimental spatial profile (user provided)
# nH = 4.66e15 * exp(-y/0.46159) + 3.81e14 (y in mm, nH in cm⁻³)
A_exp = 4.66e15  # cm⁻³
lambda_decay = 0.46159  # mm
B_exp = 3.81e14  # cm⁻³

print("Spatial H profile (experimental):")
print(f"  nH(y) = {A_exp:.2e} × exp(-y/{lambda_decay:.5f}) + {B_exp:.2e} cm⁻³")
print(f"  where y is in mm")
print(f"  Peak at y = 0.1 mm: nH ~ {A_exp * np.exp(-0.1/lambda_decay) + B_exp:.2e} cm⁻³")
print()

# Calculate H density at various positions
positions = [0.1, 0.5, 1.0, 2.0, 4.0, 10.0]  # mm
print("H density at different positions:")
print(f"{'y (mm)':<10} {'nH (cm⁻³)':<15} {'Fraction of peak':<20}")
print("-" * 50)
for y in positions:
    nH = A_exp * np.exp(-y/lambda_decay) + B_exp
    frac = nH / (A_exp * np.exp(-0.1/lambda_decay) + B_exp)
    print(f"{y:<10.1f} {nH:<15.2e} {frac:<20.3f}")
print()

# At y = 4 mm (negative glow region)
y_NG = 4.0  # mm
nH_NG = A_exp * np.exp(-y_NG/lambda_decay) + B_exp
print(f"At y = {y_NG} mm (negative glow):")
print(f"  nH = {nH_NG:.2e} cm⁻³")
print()

# Calculate H flux from gradient
# Flux = -D × ∂n/∂y
# For diffusion-dominated transport, but here we have drift too

# From the exponential profile:
# ∂n/∂y = -A_exp/λ × exp(-y/λ)
# Convert to cm: ∂n/∂y (in 1/cm⁴) = -A_exp/(λ_mm/10) × exp(-y/λ)

lambda_cm = lambda_decay / 10  # Convert mm to cm
y_NG_cm = y_NG / 10  # Convert mm to cm

gradient_at_NG = -A_exp / lambda_cm * np.exp(-y_NG/lambda_decay)  # cm⁻⁴
print(f"H density gradient at y = {y_NG} mm:")
print(f"  ∂n/∂y = {gradient_at_NG:.2e} cm⁻⁴")
print()

# Estimate diffusion coefficient for H at 500 mTorr
# D ≈ λ_mfp × v_thermal / 3
# At 500 mTorr, λ_mfp ~ 83 μm = 0.0083 cm
# v_thermal ~ 2.9e5 cm/s
lambda_mfp = 0.0083  # cm
v_thermal = 2.9e5  # cm/s
D_H = lambda_mfp * v_thermal / 3  # cm²/s
print(f"H diffusion coefficient estimate:")
print(f"  D ~ {D_H:.1f} cm²/s")
print()

# Diffusive flux at y = 4 mm
flux_diff = -D_H * gradient_at_NG  # cm⁻²/s
print(f"Diffusive H flux at y = {y_NG} mm:")
print(f"  Φ = -D × ∂n/∂y = {flux_diff:.2e} cm⁻²/s")
print()

# For 0D model: convert flux to volumetric source
# Source = Flux × A_sheath / V_NG
# Need geometry info - assume cylinder with L = 1.7 cm (discharge length)

# Actually, in a 0D model, we need to think about the NET H entering the NG volume
# from the sheath boundary. The key is: how much H is produced in the sheath per unit time
# and enters the NG?

print("=" * 80)
print("ESTIMATING 0D SOURCE TERM")
print("=" * 80)
print()

# The gradient near the sheath boundary (y ~ 0) is much steeper
y_boundary = 0.1  # mm (sheath edge)
gradient_at_boundary = -A_exp / lambda_cm * np.exp(-y_boundary/lambda_decay)  # cm⁻⁴
flux_at_boundary = -D_H * gradient_at_boundary  # cm⁻²/s

print(f"At sheath-NG boundary (y = {y_boundary} mm):")
print(f"  ∂n/∂y = {gradient_at_boundary:.2e} cm⁻⁴")
print(f"  Diffusive flux = {flux_at_boundary:.2e} cm⁻²/s")
print()

# For a cylindrical discharge
R = 1.5  # cm (typical radius)
L = 1.7  # cm (discharge length)
A_cathode = np.pi * R**2  # cm² (cathode area)
V_NG = np.pi * R**2 * L  # cm³ (NG volume, rough approximation)

H_source_volumetric = flux_at_boundary * A_cathode / V_NG  # cm⁻³/s

print(f"Cylindrical geometry (R = {R} cm, L = {L} cm):")
print(f"  Cathode area: {A_cathode:.2f} cm²")
print(f"  NG volume (approx): {V_NG:.2f} cm³")
print(f"  Volumetric source = flux × A / V = {H_source_volumetric:.2e} cm⁻³/s")
print()

# Compare to current hardcoded value
H_drift_current = 3.2e17  # cm⁻³/s
print(f"Current H_drift_gain in odefun.py: {H_drift_current:.2e} cm⁻³/s")
print(f"Ratio (current / calculated): {H_drift_current / H_source_volumetric:.1f}×")
print()

print("=" * 80)
print("ALTERNATIVE: LOSS RATE FROM STEADY-STATE PROFILE")
print("=" * 80)
print()

# At steady state in the NG, production = loss
# If nH ~ 3.81e14 cm⁻³ at y = 4 mm, and this is relatively constant in NG,
# then: Source_drift + Source_chemistry = Loss_wall + Loss_chemistry

# We can estimate the total H in the NG region and its lifetime
nH_NG_avg = nH_NG  # Average H in NG
tau_NG = 1.0e-3  # Assume ~1 ms lifetime (to be refined)

# Total H loss rate needed = nH / τ
H_loss_rate = nH_NG_avg / tau_NG  # cm⁻³/s

print(f"If average H in NG is {nH_NG_avg:.2e} cm⁻³")
print(f"And H lifetime is {tau_NG*1e3:.1f} ms:")
print(f"  Required source/sink rate: {H_loss_rate:.2e} cm⁻³/s")
print()

# From our optimization, we know chemical production is ~4e15 cm⁻³/s
# and wall loss dominates at ~1.2e16 cm⁻³/s
H_chem_production = 4.1e15  # cm⁻³/s (from balance analysis)
H_wall_loss = 1.2e16  # cm⁻³/s (from balance analysis at nH = 4.2e13)

# Scale wall loss to nH_NG
k_wall_current = 289  # s⁻¹
H_wall_loss_at_NG = k_wall_current * nH_NG_avg

print(f"From optimization balance analysis:")
print(f"  Chemical H production: {H_chem_production:.2e} cm⁻³/s")
print(f"  Wall loss at nH = {nH_NG_avg:.2e}: {H_wall_loss_at_NG:.2e} cm⁻³/s")
print(f"  Required drift source: {H_wall_loss_at_NG - H_chem_production:.2e} cm⁻³/s")
print()

print("=" * 80)
print("RECOMMENDATIONS")
print("=" * 80)
print()

H_drift_recommended = max(H_wall_loss_at_NG - H_chem_production, 0)
print(f"1. Set H_drift_gain = {H_drift_recommended:.2e} cm⁻³/s")
print(f"   (Based on steady-state balance at nH = {nH_NG_avg:.2e} cm⁻³)")
print()
print(f"2. This is {H_drift_current/H_drift_recommended:.1f}× lower than current value ({H_drift_current:.2e})")
print()
print(f"3. Update target H density to {nH_NG_avg:.2e} cm⁻³")
print(f"   (matches the NG plateau from experimental profile)")
