#!/usr/bin/env python3
"""
Test Te sensitivity to verify electron impact rate predictions

Based on CRITICAL_TE_ELECTRON_IMPACT_ANALYSIS.md:
- At Te = 1.3 eV: C2 = 4.43×10⁸ cm⁻³ (baseline)
- At Te = 3.0 eV: Predicted 76× increase → C2 ≈ 3.4×10¹⁰ cm⁻³

This script tests Te = 1.3, 1.5, 2.0, 2.5, 3.0, 4.0 eV
"""

import numpy as np
from scipy.integrate import solve_ivp
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

# Set conditions
PRESSURE_MTORR = 400
TGAS_K = 570
NE = 2.3e9

# Te values to test
TE_VALUES = [1.3, 1.5, 2.0, 2.5, 3.0, 4.0]

def pressure_to_density(pressure_mTorr, T_K):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

n_total = pressure_to_density(PRESSURE_MTORR, TGAS_K)

# Experimental targets
targets = {
    'H': 2.3e14,
    'CH': 1.34e9,
    'C2': 5.6e11
}

print("="*80)
print("Te SENSITIVITY STUDY")
print("="*80)
print()
print("Testing electron impact rate hypothesis:")
print("  - Baseline (Te=1.3 eV): C2 = 4.43×10⁸ cm⁻³")
print("  - Prediction (Te=3.0 eV): C2 ≈ 3.4×10¹⁰ cm⁻³ (76× increase)")
print()
print("="*80)
print()

results = []

for TE_EV in TE_VALUES:
    print(f"Testing Te = {TE_EV} eV...")

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

    try:
        sol = solve_ivp(ode_func, (0, 500), y0, method='BDF',
                       rtol=1e-7, atol=1e-9, max_step=1.0)

        if not sol.success:
            print(f"  ❌ Solver failed: {sol.message}")
            continue

        # Get final densities
        y_final = sol.y[:, -1]
        n_H = y_final[species.index('H')]
        n_CH = y_final[species.index('CH')]
        n_C2 = y_final[species.index('C2')]

        # Calculate errors
        error_H = abs(n_H - targets['H']) / targets['H'] * 100
        error_CH = abs(n_CH - targets['CH']) / targets['CH'] * 100
        error_C2 = abs(n_C2 - targets['C2']) / targets['C2'] * 100
        rms_error = np.sqrt((error_H**2 + error_CH**2 + error_C2**2) / 3)

        results.append({
            'Te': TE_EV,
            'H': n_H,
            'CH': n_CH,
            'C2': n_C2,
            'error_H': error_H,
            'error_CH': error_CH,
            'error_C2': error_C2,
            'rms_error': rms_error
        })

        print(f"  ✓ H  = {n_H:.2e} cm⁻³ (error: {error_H:.1f}%)")
        print(f"  ✓ CH = {n_CH:.2e} cm⁻³ (error: {error_CH:.1f}%)")
        print(f"  ✓ C2 = {n_C2:.2e} cm⁻³ (error: {error_C2:.1f}%)")
        print(f"  RMS error: {rms_error:.1f}%")
        print()

    except Exception as e:
        print(f"  ❌ Error: {e}")
        print()

print("="*80)
print("RESULTS SUMMARY")
print("="*80)
print()

# Create results table
print(f"{'Te (eV)':<8} {'H (cm⁻³)':<12} {'CH (cm⁻³)':<12} {'C2 (cm⁻³)':<12} {'C2 vs 1.3eV':<12} {'RMS Error':<10}")
print("-"*80)

baseline_C2 = results[0]['C2'] if results else 1.0

for r in results:
    c2_ratio = r['C2'] / baseline_C2
    print(f"{r['Te']:<8.1f} {r['H']:<12.2e} {r['CH']:<12.2e} {r['C2']:<12.2e} {c2_ratio:<12.1f}× {r['rms_error']:<10.1f}%")

print()
print("="*80)
print("ANALYSIS")
print("="*80)
print()

# Find best Te
if results:
    best = min(results, key=lambda x: x['rms_error'])
    print(f"✓ Best Te: {best['Te']} eV (RMS error: {best['rms_error']:.1f}%)")
    print()

    # Compare with predictions
    print("Prediction vs. Reality:")
    print(f"  Baseline (Te=1.3 eV): C2 = {results[0]['C2']:.2e} cm⁻³")

    te3_result = [r for r in results if r['Te'] == 3.0]
    if te3_result:
        r = te3_result[0]
        predicted_ratio = 76.8  # From analysis
        actual_ratio = r['C2'] / results[0]['C2']
        print(f"  Te=3.0 eV:")
        print(f"    Predicted: 76× increase → {results[0]['C2']*76:.2e} cm⁻³")
        print(f"    Actual:    {actual_ratio:.1f}× increase → {r['C2']:.2e} cm⁻³")
        print(f"    Match: {actual_ratio/predicted_ratio*100:.1f}% of predicted")

    print()

    # Check if any Te gives good match
    print("Distance from experimental targets:")
    for r in results:
        gap_H = targets['H'] / r['H']
        gap_CH = targets['CH'] / r['CH']
        gap_C2 = targets['C2'] / r['C2']
        print(f"  Te={r['Te']} eV: H {gap_H:.1f}×, CH {gap_CH:.1f}×, C2 {gap_C2:.1f}×")

print()
print("="*80)
