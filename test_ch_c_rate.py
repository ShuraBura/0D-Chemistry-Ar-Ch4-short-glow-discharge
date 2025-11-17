#!/usr/bin/env python3
"""
Test hypothesis: CH + C → C2 rate constants are underestimated

Current rates: k ~ 1-1.2×10⁻¹⁰ cm³/s (no citations!)
These contribute ~38% of C2 production

Will test multiplying by: 1×, 10×, 100×, 1000×
"""

import numpy as np
from scipy.integrate import solve_ivp
from build_reactions import build_reactions
from define_rates import define_rates
from odefun_optimized import PlasmaODE_Optimized

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

print("="*80)
print("TESTING CH + C → C2 RATE CONSTANT SENSITIVITY")
print("="*80)
print()
print("Hypothesis: CH + C → C2 rate constants might be underestimated")
print("Current values: k ~ 1-1.2×10⁻¹⁰ cm³/s (NO CITATIONS)")
print("These contribute ~38% of total C2 production")
print()
print("Testing multipliers: 1×, 10×, 100×, 1000×")
print("="*80)
print()

multipliers = [1, 10, 100, 1000]
results = []

for mult in multipliers:
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
            'ArHPlus': 1.0, 'CH3Minus': 1.0, 'C2': 0.001, 'CH': 0.001, 'H': 0.0,
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

    # Build reactions
    k = define_rates(params)

    # Modify CH + C rates
    k['CH_C_C2_H_cm3_7_9'] *= mult
    k['CH_C_C2_H2_cm3_7_24'] *= mult
    k['CH_C_C2_H_cm3_7_43'] *= mult
    k['C_CH_C2_H_cm3_7_4'] *= mult  # Also C + CH (same reaction)

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

    print(f"Testing multiplier = {mult}×...")

    try:
        sol = solve_ivp(ode_func, (0, 500), y0, method='BDF',
                       rtol=1e-7, atol=1e-9, max_step=1.0)

        if not sol.success:
            print(f"  ❌ Solver failed: {sol.message}")
            print()
            continue

        # Get final densities
        y_final = sol.y[:, -1]
        n_H = y_final[species.index('H')]
        n_CH = y_final[species.index('CH')]
        n_C2 = y_final[species.index('C2')]

        # Targets
        target_H = 2.3e14
        target_CH = 1.34e9
        target_C2 = 5.6e11

        # Errors
        error_H = abs(n_H - target_H) / target_H * 100
        error_CH = abs(n_CH - target_CH) / target_CH * 100
        error_C2 = abs(n_C2 - target_C2) / target_C2 * 100
        rms_error = np.sqrt((error_H**2 + error_CH**2 + error_C2**2) / 3)

        results.append({
            'mult': mult,
            'H': n_H,
            'CH': n_CH,
            'C2': n_C2,
            'rms_error': rms_error
        })

        print(f"  H  = {n_H:.2e} cm⁻³ (error: {error_H:.1f}%)")
        print(f"  CH = {n_CH:.2e} cm⁻³ (error: {error_CH:.1f}%)")
        print(f"  C2 = {n_C2:.2e} cm⁻³ (error: {error_C2:.1f}%)")
        print(f"  RMS error: {rms_error:.1f}%")
        print()

    except Exception as e:
        print(f"  ❌ Error: {e}")
        print()

print("="*80)
print("SUMMARY")
print("="*80)
print()

if results:
    print(f"{'Multiplier':<12} {'C2 (cm⁻³)':<15} {'vs Target':<15} {'RMS Error':<15}")
    print("-"*80)

    for r in results:
        gap = 5.6e11 / r['C2']
        print(f"{r['mult']:<12}× {r['C2']:<15.2e} {gap:<15.1f}× {r['rms_error']:<15.1f}%")

    print()
    print("Target: C2 = 5.6×10¹¹ cm⁻³")
    print()

    # Find best
    best = min(results, key=lambda x: x['rms_error'])
    print(f"✓ Best multiplier: {best['mult']}× (RMS error: {best['rms_error']:.1f}%)")

    # Check if any get close to target
    for r in results:
        if abs(r['C2'] - 5.6e11) / 5.6e11 < 0.5:  # Within factor of 2
            print(f"✓✓✓ Multiplier {r['mult']}× gets C2 within factor of 2 of target!")

print()
print("="*80)
