#!/usr/bin/env python3
"""
Pressure scan: 340-500 mTorr

At T = 400 K (typical discharge temperature):
- 340 mTorr → n_total = 8.21e15 cm⁻³
- 400 mTorr → n_total = 9.66e15 cm⁻³ (baseline)
- 500 mTorr → n_total = 1.21e16 cm⁻³

Gas mixture: 85% Ar, 15% CH4 (constant)

This tests how pressure affects:
- Collision frequency (∝ pressure)
- 2-body vs 3-body reaction balance
- Wall loss vs volume reaction rates
- CH, C2, H densities
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


# Experimental targets
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}


def pressure_to_density(pressure_mTorr, T_K=400):
    """
    Convert pressure to number density.

    n = P / (kB * T)

    Parameters:
    -----------
    pressure_mTorr : float
        Pressure in mTorr
    T_K : float
        Temperature in Kelvin (default 400 K for discharge)

    Returns:
    --------
    float : number density in cm⁻³
    """
    # Constants
    kB = 1.38064852e-23  # J/K
    Torr_to_Pa = 133.322

    # Convert pressure to Pa
    P_Pa = pressure_mTorr * 1e-3 * Torr_to_Pa

    # Calculate density in m⁻³
    n_m3 = P_Pa / (kB * T_K)

    # Convert to cm⁻³
    n_cm3 = n_m3 * 1e-6

    return n_cm3


def run_simulation_at_pressure(pressure_mTorr, Te=1.0, ne=3.3e9, E_field=50.0, verbose=False):
    """
    Run simulation at specified pressure.

    Returns:
    --------
    dict : {'H': value, 'CH': value, 'C2': value, 'C2H2': value, ...} or None if failed
    """
    try:
        # Calculate total neutral density
        n_total = pressure_to_density(pressure_mTorr)

        # Setup parameters
        params = {
            'E_field': E_field,
            'L_discharge': 0.45,
            'ne': ne,
            'Te': Te,
            'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                        'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                        'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                        'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                        'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
            'mobilities': {
                'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
                'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
                'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
                'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
            }
        }

        # Load rates with temperature dependence
        k = define_rates(params)
        db = get_complete_rate_database()

        # Ensure all rates within literature bounds
        for name, rate_db in db.items():
            if name in k:
                if k[name] < rate_db.min:
                    k[name] = rate_db.min
                elif k[name] > rate_db.max:
                    k[name] = rate_db.max

        params['k'] = k
        params['R'], params['tags'] = build_reactions(params)

        # Initial conditions
        species = params['species']
        ns = len(species)
        y0 = np.ones(ns) * 1e3

        def set_density(name, value):
            try:
                idx = species.index(name)
                y0[idx] = value
            except ValueError:
                pass

        # Set neutral gas densities based on pressure and mixture
        n_Ar = 0.85 * n_total
        n_CH4 = 0.15 * n_total

        set_density('e', ne)
        set_density('Ar', n_Ar)
        set_density('CH4', n_CH4)
        set_density('ArPlus', 1e7)
        set_density('CH4Plus', 1e5)
        set_density('CH3Plus', 1e5)
        set_density('CH5Plus', 1e3)
        set_density('ArHPlus', 5e5)
        set_density('CH3Minus', 5e4)
        set_density('H2', 1e12)
        set_density('ArStar', 5e6)
        set_density('H', 1e11)
        set_density('C2', 5e7)
        set_density('CH', 5e4)
        set_density('C2H4', 5e7)
        set_density('C2H6', 1e6)
        set_density('CH2', 1e11)
        set_density('C2H2', 1e12)
        set_density('C2H5', 1e6)
        set_density('CH3', 5e7)
        set_density('C', 5e7)

        # Run simulation
        if verbose:
            print(f"  Running at {pressure_mTorr:.0f} mTorr (n_total = {n_total:.2e} cm⁻³)...")

        ode_func = PlasmaODE_Optimized(params)
        sol = solve_ivp(
            ode_func,
            (0, 100),
            y0,
            method='BDF',
            rtol=1e-5,
            atol=1e-6,
            max_step=10.0
        )

        if not sol.success:
            if verbose:
                print(f"    Solver failed: {sol.message}")
            return None

        # Extract final densities
        y_final = sol.y[:, -1]

        def get_density(name):
            try:
                idx = species.index(name)
                return y_final[idx]
            except ValueError:
                return 0.0

        results = {
            'pressure': pressure_mTorr,
            'n_total': n_total,
            'H': get_density('H'),
            'CH': get_density('CH'),
            'C2': get_density('C2'),
            'C2H2': get_density('C2H2'),
            'CH3': get_density('CH3'),
            'CH2': get_density('CH2'),
            'C': get_density('C'),
            'Ar': get_density('Ar'),
            'CH4': get_density('CH4'),
        }

        if verbose:
            print(f"    ✓ Success")

        return results

    except Exception as e:
        if verbose:
            print(f"    Simulation failed: {e}")
        return None


def main():
    print("=" * 80)
    print(" PRESSURE SCAN: 340-500 mTorr")
    print("=" * 80)
    print("\nGas mixture: 85% Ar, 15% CH4 (constant)")
    print("Temperature: 400 K (assumed)")
    print("\nPressure range:")
    print("  340 mTorr → n_total = 8.21e15 cm⁻³")
    print("  400 mTorr → n_total = 9.66e15 cm⁻³ (baseline)")
    print("  500 mTorr → n_total = 1.21e16 cm⁻³")

    # Use baseline parameters from previous optimizations
    # Or use optimized parameters from Te-Ne-E optimization
    Te = 1.0      # eV - baseline
    ne = 3.3e9    # cm⁻³ - baseline
    E_field = 50.0  # V/cm - baseline

    print(f"\nFixed parameters:")
    print(f"  Te: {Te} eV")
    print(f"  Ne: {ne:.2e} cm⁻³")
    print(f"  E-field: {E_field} V/cm")

    # Pressure scan points
    pressures = np.linspace(340, 500, 15)

    print(f"\n\nScanning {len(pressures)} pressure points from 340 to 500 mTorr...")
    print("=" * 80)

    results_list = []

    for pressure in pressures:
        result = run_simulation_at_pressure(pressure, Te, ne, E_field, verbose=True)
        if result is not None:
            results_list.append(result)

    if not results_list:
        print("\n❌ All simulations failed!")
        return

    print(f"\n✓ Successfully completed {len(results_list)}/{len(pressures)} simulations")

    # Extract data for plotting
    pressures_success = [r['pressure'] for r in results_list]
    H_densities = [r['H'] for r in results_list]
    CH_densities = [r['CH'] for r in results_list]
    C2_densities = [r['C2'] for r in results_list]
    C2H2_densities = [r['C2H2'] for r in results_list]
    CH3_densities = [r['CH3'] for r in results_list]

    # Calculate ratios vs targets
    H_ratios = [h / TARGETS['H'] for h in H_densities]
    CH_ratios = [ch / TARGETS['CH'] for ch in CH_densities]
    C2_ratios = [c2 / TARGETS['C2'] for c2 in C2_densities]

    # Print summary table
    print("\n" + "=" * 80)
    print(" RESULTS SUMMARY")
    print("=" * 80)
    print(f"\n{'P (mTorr)':<12} {'H (cm⁻³)':<14} {'H ratio':<10} {'CH (cm⁻³)':<14} {'CH ratio':<10} {'C2 (cm⁻³)':<14} {'C2 ratio':<10}")
    print("-" * 100)

    for i, p in enumerate(pressures_success):
        print(f"{p:<12.0f} {H_densities[i]:<14.2e} {H_ratios[i]:<10.2f} "
              f"{CH_densities[i]:<14.2e} {CH_ratios[i]:<10.1f} "
              f"{C2_densities[i]:<14.2e} {C2_ratios[i]:<10.2f}")

    # Find best pressure for each species
    best_H_idx = min(range(len(H_ratios)), key=lambda i: abs(H_ratios[i] - 1.0))
    best_CH_idx = min(range(len(CH_ratios)), key=lambda i: abs(CH_ratios[i] - 1.0))
    best_C2_idx = min(range(len(C2_ratios)), key=lambda i: abs(C2_ratios[i] - 1.0))

    print("\n" + "=" * 80)
    print(" BEST PRESSURE FOR EACH TARGET")
    print("=" * 80)
    print(f"\nH:  {pressures_success[best_H_idx]:.0f} mTorr (ratio: {H_ratios[best_H_idx]:.2f}x)")
    print(f"CH: {pressures_success[best_CH_idx]:.0f} mTorr (ratio: {CH_ratios[best_CH_idx]:.1f}x)")
    print(f"C2: {pressures_success[best_C2_idx]:.0f} mTorr (ratio: {C2_ratios[best_C2_idx]:.2f}x)")

    # Create plots
    print("\nGenerating plots...")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Absolute densities
    ax1.plot(pressures_success, H_densities, 'o-', label='H', linewidth=2)
    ax1.axhline(TARGETS['H'], color='blue', linestyle='--', alpha=0.5, label='H target')
    ax1.plot(pressures_success, CH_densities, 's-', label='CH', linewidth=2)
    ax1.axhline(TARGETS['CH'], color='orange', linestyle='--', alpha=0.5, label='CH target')
    ax1.plot(pressures_success, C2_densities, '^-', label='C2', linewidth=2)
    ax1.axhline(TARGETS['C2'], color='green', linestyle='--', alpha=0.5, label='C2 target')
    ax1.set_xlabel('Pressure (mTorr)', fontsize=12)
    ax1.set_ylabel('Density (cm⁻³)', fontsize=12)
    ax1.set_title('Species Densities vs Pressure', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')

    # Plot 2: Ratios vs target
    ax2.plot(pressures_success, H_ratios, 'o-', label='H / target', linewidth=2)
    ax2.plot(pressures_success, CH_ratios, 's-', label='CH / target', linewidth=2)
    ax2.plot(pressures_success, C2_ratios, '^-', label='C2 / target', linewidth=2)
    ax2.axhline(1.0, color='black', linestyle='--', alpha=0.5, label='Perfect match')
    ax2.set_xlabel('Pressure (mTorr)', fontsize=12)
    ax2.set_ylabel('Ratio vs Target', fontsize=12)
    ax2.set_title('Target Matching vs Pressure', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')

    # Plot 3: C2H2 and CH3 (precursors)
    ax3.plot(pressures_success, C2H2_densities, 'o-', label='C2H2', linewidth=2)
    ax3.plot(pressures_success, CH3_densities, 's-', label='CH3', linewidth=2)
    ax3.set_xlabel('Pressure (mTorr)', fontsize=12)
    ax3.set_ylabel('Density (cm⁻³)', fontsize=12)
    ax3.set_title('Key Precursor Densities', fontsize=14, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale('log')

    # Plot 4: Error metric
    # Calculate weighted error (same as optimization)
    errors = []
    for i in range(len(pressures_success)):
        error = (
            1.0 * (H_ratios[i] - 1.0)**2 +
            20.0 * (CH_ratios[i] - 1.0)**2 +
            3.0 * (C2_ratios[i] - 1.0)**2
        )
        errors.append(error)

    ax4.plot(pressures_success, errors, 'o-', linewidth=2, color='red')
    best_error_idx = np.argmin(errors)
    ax4.plot(pressures_success[best_error_idx], errors[best_error_idx],
             'g*', markersize=20, label=f'Best: {pressures_success[best_error_idx]:.0f} mTorr')
    ax4.set_xlabel('Pressure (mTorr)', fontsize=12)
    ax4.set_ylabel('Weighted Error', fontsize=12)
    ax4.set_title('Overall Error vs Pressure', fontsize=14, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_yscale('log')

    plt.tight_layout()
    plt.savefig('pressure_scan_results.png', dpi=150, bbox_inches='tight')
    print("✓ Plot saved to: pressure_scan_results.png")

    # Final recommendations
    print("\n" + "=" * 80)
    print(" PRESSURE EFFECT ANALYSIS")
    print("=" * 80)

    # Calculate trends
    H_trend = (H_densities[-1] - H_densities[0]) / H_densities[0] * 100
    CH_trend = (CH_densities[-1] - CH_densities[0]) / CH_densities[0] * 100
    C2_trend = (C2_densities[-1] - C2_densities[0]) / C2_densities[0] * 100

    print(f"\nDensity changes from 340 to 500 mTorr:")
    print(f"  H:  {H_trend:+.1f}%")
    print(f"  CH: {CH_trend:+.1f}%")
    print(f"  C2: {C2_trend:+.1f}%")

    print(f"\nBest overall pressure: {pressures_success[best_error_idx]:.0f} mTorr")
    print(f"  Error: {errors[best_error_idx]:.1f} (baseline at 400 mTorr: {errors[min(range(len(pressures_success)), key=lambda i: abs(pressures_success[i]-400))]:.1f})")

    if errors[best_error_idx] < errors[min(range(len(pressures_success)), key=lambda i: abs(pressures_success[i]-400))]:
        improvement = (1 - errors[best_error_idx] / errors[min(range(len(pressures_success)), key=lambda i: abs(pressures_success[i]-400))]) * 100
        print(f"  → {improvement:.1f}% improvement over baseline!")
    else:
        print(f"  → No significant improvement over baseline")

    print("\n" + "=" * 80)


if __name__ == '__main__':
    main()
