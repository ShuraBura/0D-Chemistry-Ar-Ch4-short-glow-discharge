#!/usr/bin/env python3
"""
Ar/CH4 Plasma Simulation - Optimized for Cathode Glow (CG) Region
Targets: TALIF measurements at cathode glow
"""

import numpy as np
from scipy.integrate import solve_ivp
import time
import matplotlib.pyplot as plt

from define_rates import define_rates
from build_reactions import build_reactions
from odefun import PlasmaODE


# TALIF Experimental Targets (Cathode Glow)
TARGETS = {
    'H': 8.57e15,    # cm⁻³
    'CH': 2.75e8,    # cm⁻³
    'C2': 1.12e11,   # cm⁻³
}

# Derived ratios (more robust than absolute values)
TARGET_RATIOS = {
    'H/CH': TARGETS['H'] / TARGETS['CH'],      # = 3.1e7
    'C2/CH': TARGETS['C2'] / TARGETS['CH'],    # = 408
    'H/C2': TARGETS['H'] / TARGETS['C2'],      # = 7.6e4
}


def calculate_initial_densities(params):
    """Calculate initial species densities."""
    species = params['species']
    ns = len(species)
    y0 = np.ones(ns) * 1e3  # Default minimum density

    # Helper function to set density
    def set_density(name, value):
        try:
            idx = species.index(name)
            y0[idx] = value
        except ValueError:
            pass

    # Set specific initial densities
    set_density('e', params['ne'])
    set_density('Ar', 0.85 * 9.66e15)
    set_density('CH4', 0.15 * 9.66e15)
    set_density('ArPlus', 5e8)      # Lower for CG
    set_density('CH4Plus', 1e6)
    set_density('CH3Plus', 1e6)
    set_density('CH5Plus', 1e4)
    set_density('ArHPlus', 1e6)
    set_density('CH3Minus', 5e5)
    set_density('H2', 1e12)
    set_density('ArStar', 1e7)
    set_density('H', 1e15)          # Start near target
    set_density('C2', 1e11)         # Start near target
    set_density('CH', 1e8)          # Start near target
    set_density('C2H4', 5e7)
    set_density('C2H6', 1e6)
    set_density('CH2', 1e11)
    set_density('C2H2', 1e12)
    set_density('C2H5', 1e6)
    set_density('CH3', 5e7)
    set_density('C', 5e7)

    return y0


def print_comparison(sol, params):
    """Print comparison with experimental targets."""
    species = params['species']

    print("\n" + "=" * 70)
    print("COMPARISON WITH TALIF TARGETS (Cathode Glow)")
    print("=" * 70)

    # Get final densities
    results = {}
    for name in ['H', 'CH', 'C2']:
        try:
            idx = species.index(name)
            results[name] = sol.y[idx, -1]
        except ValueError:
            results[name] = 0.0

    # Print absolute comparison
    print("\nAbsolute Densities:")
    print(f"{'Species':<10} {'Target (cm⁻³)':<20} {'Model (cm⁻³)':<20} {'Ratio (M/T)':<15}")
    print("-" * 70)
    for name in ['H', 'CH', 'C2']:
        target = TARGETS[name]
        model = results[name]
        ratio = model / target if target > 0 else 0
        print(f"{name:<10} {target:>18.2e}   {model:>18.2e}   {ratio:>13.2f}x")

    # Print ratio comparison (more robust!)
    print("\nDensity Ratios (More Robust):")
    print(f"{'Ratio':<15} {'Target':<20} {'Model':<20} {'Agreement':<15}")
    print("-" * 70)

    model_ratios = {
        'H/CH': results['H'] / results['CH'] if results['CH'] > 0 else 0,
        'C2/CH': results['C2'] / results['CH'] if results['CH'] > 0 else 0,
        'H/C2': results['H'] / results['C2'] if results['C2'] > 0 else 0,
    }

    for ratio_name in ['H/CH', 'C2/CH', 'H/C2']:
        target = TARGET_RATIOS[ratio_name]
        model = model_ratios[ratio_name]
        agreement = model / target if target > 0 else 0
        print(f"{ratio_name:<15} {target:>18.2e}   {model:>18.2e}   {agreement:>13.2f}x")

    # Overall assessment
    print("\n" + "=" * 70)
    print("ASSESSMENT:")
    print("=" * 70)

    # Check if within factor of 2-10 (good for plasma models!)
    h_ok = 0.1 < (results['H'] / TARGETS['H']) < 10
    ch_ok = 0.1 < (results['CH'] / TARGETS['CH']) < 10
    c2_ok = 0.1 < (results['C2'] / TARGETS['C2']) < 10

    print(f"H  density: {'✓ Good' if h_ok else '✗ Needs tuning'}")
    print(f"CH density: {'✓ Good' if ch_ok else '✗ Needs tuning'}")
    print(f"C2 density: {'✓ Good' if c2_ok else '✗ Needs tuning'}")

    if h_ok and ch_ok and c2_ok:
        print("\n🎉 MODEL VALIDATED against TALIF targets!")
    else:
        print("\n⚠️  Model needs parameter tuning")
        print("    Try adjusting: Te (0.5-2.5 eV), E-field (50-120 V/cm), or ne")


def main():
    """Main simulation function optimized for CG."""
    print("=" * 70)
    print("Ar/CH4 Plasma Simulation - CATHODE GLOW (CG) Region")
    print("=" * 70)
    print(f"\nTALIF Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm⁻³")
    print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
    print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")
    print(f"\nTarget Ratios:")
    print(f"  H/CH:  {TARGET_RATIOS['H/CH']:.2e}")
    print(f"  C2/CH: {TARGET_RATIOS['C2/CH']:.2e}")

    # Simulation parameters (CG-optimized)
    params = {
        'P': 0.4,                # Pressure in Torr
        'Tg': 400,               # Gas temperature in K
        'n_tot': 9.66e15,        # Total density in cm^-3
        'ne': 5e9,               # Electron density in cm^-3 (LOW for CG)
        'Te': 1.5,               # Electron temperature in eV (effective for non-thermal EEDF)
        'E_field': 80,           # V/cm (HIGH for cathode region)
        'L_discharge': 0.45,     # cm
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                    'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                    'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                    'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'C2H2Star'],
        'ion_species': ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH3Minus', 'H3Plus', 'CHPlus',
                        'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus'],
        'mobilities': {
            'ArPlus': 3057.28,
            'CH4Plus': 6432,
            'CH3Plus': 4949.6,
            'CH5Plus': 4761.6,
            'ArHPlus': 2969.6,
            'CH2Plus': 4949.6,
            'C2H5Plus': 4949.6,
            'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6,
            'C2HPlus': 5000,
            'H3Plus': 5000,
            'CHPlus': 5000,
            'CH3Minus': 3000,
            'HMinus': 3000
        }
    }

    print(f"\nModel Parameters (CG Region):")
    print(f"  Pressure:    {params['P']} Torr")
    print(f"  Tgas:        {params['Tg']} K")
    print(f"  ne:          {params['ne']:.1e} cm⁻³ (LOW - cathode glow)")
    print(f"  Te:          {params['Te']} eV (effective, non-thermal EEDF)")
    print(f"  E-field:     {params['E_field']} V/cm (HIGH - cathode)")
    print(f"  L_discharge: {params['L_discharge']} cm")

    # Build reaction network
    print("\nBuilding reaction network...")
    params['k'] = define_rates(params)
    params['R'], params['tags'] = build_reactions(params)
    print(f"  Number of species: {len(params['species'])}")
    print(f"  Number of reactions: {len(params['R'])}")

    # Initial conditions
    print("\nSetting initial conditions...")
    y0 = calculate_initial_densities(params)

    # Time span
    t_span = (0, 100)  # 0 to 100 seconds
    t_eval = np.logspace(-3, 2, 50)  # Logarithmic time points

    # Create ODE function
    ode_func = PlasmaODE(params)

    # Solver options
    print("\nStarting simulation...")
    print(f"  Time span: {t_span[0]} to {t_span[1]} s")
    print(f"  Method: BDF (stiff solver)")

    start_time = time.time()

    # Solve ODE system
    sol = solve_ivp(
        ode_func,
        t_span,
        y0,
        method='BDF',
        t_eval=t_eval,
        rtol=1e-7,
        atol=1e-8,
    )

    elapsed_time = time.time() - start_time

    print(f"\nSimulation completed in {elapsed_time:.2f} seconds")
    print(f"  Status: {'Success' if sol.success else 'Failed'}")
    print(f"  Number of function evaluations: {sol.nfev}")

    # Compare with targets
    print_comparison(sol, params)

    # Plot results
    print("\nGenerating plots...")
    plot_results(sol, params)

    print("\nDone! Check comparison above and plots in plasma_cg_results.png")


def plot_results(sol, params):
    """Plot time evolution with target lines."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Key radicals with target lines
    ax = axes[0, 0]
    species_plot = ['H', 'CH', 'C2']
    colors = ['blue', 'red', 'green']

    for species_name, color in zip(species_plot, colors):
        try:
            idx = params['species'].index(species_name)
            ax.loglog(sol.t, sol.y[idx, :], 'o-', label=f'{species_name} (model)',
                     color=color, markersize=3, alpha=0.7)
            # Add target line
            target = TARGETS[species_name]
            ax.axhline(target, color=color, linestyle='--', linewidth=2,
                      label=f'{species_name} (TALIF target)', alpha=0.8)
        except (ValueError, KeyError):
            pass

    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Density (cm⁻³)', fontsize=12)
    ax.set_title('Key Radicals vs TALIF Targets', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Ratios
    ax = axes[0, 1]
    try:
        idx_h = params['species'].index('H')
        idx_ch = params['species'].index('CH')
        idx_c2 = params['species'].index('C2')

        ratio_h_ch = sol.y[idx_h, :] / (sol.y[idx_ch, :] + 1e-30)
        ratio_c2_ch = sol.y[idx_c2, :] / (sol.y[idx_ch, :] + 1e-30)

        ax.loglog(sol.t, ratio_h_ch, 'o-', label='H/CH (model)', markersize=3)
        ax.axhline(TARGET_RATIOS['H/CH'], color='blue', linestyle='--', linewidth=2,
                  label='H/CH (target)')

        ax.loglog(sol.t, ratio_c2_ch, 's-', label='C2/CH (model)', markersize=3)
        ax.axhline(TARGET_RATIOS['C2/CH'], color='red', linestyle='--', linewidth=2,
                  label='C2/CH (target)')

        ax.set_xlabel('Time (s)', fontsize=12)
        ax.set_ylabel('Ratio', fontsize=12)
        ax.set_title('Density Ratios vs Targets', fontsize=13, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
    except ValueError:
        pass

    # Stable molecules
    key_molecules = ['H2', 'CH4', 'C2H2', 'C2H4', 'C2H6']
    ax = axes[1, 0]
    for species_name in key_molecules:
        try:
            idx = params['species'].index(species_name)
            ax.loglog(sol.t, sol.y[idx, :], 'o-', label=species_name, markersize=3)
        except ValueError:
            pass
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Density (cm⁻³)', fontsize=12)
    ax.set_title('Stable Molecules', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Major ions
    key_ions = ['ArPlus', 'CH3Plus', 'CH5Plus', 'e']
    ax = axes[1, 1]
    for species_name in key_ions:
        try:
            idx = params['species'].index(species_name)
            ax.loglog(sol.t, sol.y[idx, :], 'o-', label=species_name, markersize=3)
        except ValueError:
            pass
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Density (cm⁻³)', fontsize=12)
    ax.set_title('Major Ions', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('plasma_cg_results.png', dpi=150)
    print("  Plot saved to: plasma_cg_results.png")


if __name__ == '__main__':
    main()
