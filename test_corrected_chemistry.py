#!/usr/bin/env python3
"""
Test the corrected C2 chemistry to see the impact of Baulch fixes
"""

import numpy as np
from scipy.integrate import solve_ivp
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

# Set conditions
PRESSURE_MTORR = 400
TGAS_K = 570
NE = 2.3e9
TE_EV = 1.3

def pressure_to_density(pressure_mTorr, T_K):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

n_total = pressure_to_density(PRESSURE_MTORR, TGAS_K)

params = {
    'P': PRESSURE_MTORR,
    'Tgas': TGAS_K,
    'Te': TE_EV,
    'ne': NE,
    'E_field': 150,
    'L_discharge': 0.45,
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
    },
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4',
                'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3',
                'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3',
                'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'C3H5', 'HMinus', 'C2H5Plus',
                'CH2Plus', 'C4H', 'H2Plus', 'C2H6Plus', 'C2H2Star', 'C3H6'],
    'stick_coeffs': {
        'ArPlus': 1.0, 'CH4Plus': 1.0, 'CH3Plus': 1.0, 'CH5Plus': 1.0,
        'ArHPlus': 1.0, 'CH3Minus': 1.0, 'C2': 0.01, 'CH': 0.001, 'H': 0.0,
        'C': 0.01, 'CH2': 0.001, 'CH3': 0.001, 'C2H': 0.001, 'C2H2': 0.001,
        'C2H3': 0.001, 'C2H4': 0.001, 'C2H5': 0.001, 'C2H6': 0.001, 'C3H2': 0.001,
        'CHPlus': 1.0, 'C3H': 0.001, 'C4H2': 0.001, 'C3H3': 0.001, 'C3H4': 0.001,
        'C3': 0.001, 'C2H4Plus': 1.0, 'C2H3Plus': 1.0, 'C2HPlus': 1.0,
        'C3H5': 0.001, 'HMinus': 1.0, 'C2H5Plus': 1.0, 'CH2Plus': 1.0,
        'C4H': 0.001, 'H3Plus': 1.0, 'H2Plus': 1.0, 'C2H6Plus': 1.0
    },
    'n_Ar': n_total * 0.97,
    'n_CH4': n_total * 0.03,
}

print("="*80)
print("TESTING CORRECTED C2 CHEMISTRY")
print("="*80)
print()
print("CORRECTIONS APPLIED:")
print("1. ❌ REMOVED: C2 + H → CH + C (k set to 0)")
print("2. ❌ REMOVED: CH + CH2 → C2 + H2 + H (k set to 0)")
print("3. ✓ FIXED: CH2 + CH2 → C2H2 + H2 (was C2 + H2 + H2)")
print("4. ✓ FIXED: CH + CH → C2H2 (was C2 + H2)")
print("5. ✓ FIXED: C2H2 + C → C2: rate 1.0e-10 → 2.0e-10")
print("6. ✓ FIXED: CH + H → C + H2: rate 1.2e-10 → 2.0e-10")
print("7. ✓ FIXED: H + C2H2 → C2: Arrhenius form (was 7.5e9× too high)")
print()

# Build reactions (need to import define_rates)
from define_rates import define_rates
k = define_rates(params)
params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Initialize
species = params['species']
y0 = np.ones(len(species)) * 1e3
y0[species.index('Ar')] = n_total * 0.97
y0[species.index('CH4')] = n_total * 0.03
y0[species.index('e')] = NE

# Run simulation
ode_func = PlasmaODE_Optimized(params)
ode_func.H_drift_gain = 5.7e16

print("Solving for steady state...")
sol = solve_ivp(ode_func, (0, 500), y0, method='BDF',
               rtol=1e-7, atol=1e-9, max_step=1.0)

if not sol.success:
    print(f"❌ Solver failed: {sol.message}")
    exit(1)

# Get final densities
y_final = sol.y[:, -1]

# Extract key species
n_H = y_final[species.index('H')]
n_CH = y_final[species.index('CH')]
n_C2 = y_final[species.index('C2')]

print()
print("="*80)
print("RESULTS")
print("="*80)
print()
print("SPECIES DENSITIES:")
print(f"  H  = {n_H:.2e} cm⁻³  (target: 2.30e14, ratio: {n_H/2.3e14:.3f})")
print(f"  CH = {n_CH:.2e} cm⁻³  (target: 1.34e9, ratio: {n_CH/1.34e9:.3f})")
print(f"  C2 = {n_C2:.2e} cm⁻³  (target: 5.60e11, ratio: {n_C2/5.6e11:.3f})")
print()

# Comparison with old
old_c2 = 2.75e8
print(f"Old C2 (with wrong rates): {old_c2:.2e} cm⁻³")
print(f"New C2 (with corrections):  {n_C2:.2e} cm⁻³")
print(f"Improvement factor: {n_C2/old_c2:.1f}×")
print()

# Calculate errors
targets = {'H': 2.3e14, 'CH': 1.34e9, 'C2': 5.6e11}
results = {'H': n_H, 'CH': n_CH, 'C2': n_C2}

errors = []
for sp in ['H', 'CH', 'C2']:
    error_pct = 100 * abs(results[sp] - targets[sp]) / targets[sp]
    errors.append(error_pct)
    status = "✓" if error_pct < 50 else "✗"
    print(f"{status} {sp:4s}: error = {error_pct:6.1f}%")

rms_error = np.sqrt(np.mean(np.array(errors)**2))
print()
print(f"RMS Error: {rms_error:.1f}%")
print()

if n_C2 > old_c2:
    print(f"✓ C2 increased as expected!")
else:
    print(f"✗ Warning: C2 did not increase")

print("="*80)
