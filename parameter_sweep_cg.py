#!/usr/bin/env python3
"""
Parameter Sweep for CG Region
Test different Te (0.5-10 eV) and ne (1e7-1e10 cm⁻³) to match TALIF targets
"""

import numpy as np
from scipy.integrate import solve_ivp
import time
import matplotlib.pyplot as plt
from multiprocessing import Pool
import itertools

from define_rates import define_rates
from build_reactions import build_reactions
from odefun import PlasmaODE


# TALIF Experimental Targets (Cathode Glow)
TARGETS = {
    'H': 8.57e15,    # cm⁻³
    'CH': 2.75e8,    # cm⁻³
    'C2': 1.12e11,   # cm⁻³
}


def calculate_initial_densities(params):
    """Calculate initial species densities."""
    species = params['species']
    ns = len(species)
    y0 = np.ones(ns) * 1e3

    def set_density(name, value):
        try:
            idx = species.index(name)
            y0[idx] = value
        except ValueError:
            pass

    set_density('e', params['ne'])
    set_density('Ar', 0.85 * 9.66e15)
    set_density('CH4', 0.15 * 9.66e15)
    set_density('ArPlus', params['ne'] * 0.5)  # Scale with ne
    set_density('CH4Plus', params['ne'] * 0.05)
    set_density('CH3Plus', params['ne'] * 0.05)
    set_density('CH5Plus', params['ne'] * 0.01)
    set_density('ArHPlus', params['ne'] * 0.1)
    set_density('CH3Minus', params['ne'] * 0.01)
    set_density('H2', 1e12)
    set_density('ArStar', params['ne'] * 0.5)
    set_density('H', 1e15)
    set_density('C2', 1e11)
    set_density('CH', 1e8)
    set_density('C2H4', 5e7)
    set_density('C2H6', 1e6)
    set_density('CH2', 1e11)
    set_density('C2H2', 1e12)
    set_density('C2H5', 1e6)
    set_density('CH3', 5e7)
    set_density('C', 5e7)

    return y0


def run_single_simulation(params_tuple):
    """Run a single simulation with given Te and ne."""
    Te, ne, run_id = params_tuple

    try:
        # Base parameters
        params = {
            'P': 0.4,
            'Tg': 400,
            'n_tot': 9.66e15,
            'ne': ne,
            'Te': Te,
            'E_field': 80,
            'L_discharge': 0.45,
            'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                        'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                        'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                        'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                        'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'C2H2Star'],
            'ion_species': ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH3Minus',
                            'H3Plus', 'CHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus',
                            'HMinus', 'C2HPlus'],
            'mobilities': {
                'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
                'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
                'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
                'CH3Minus': 3000, 'HMinus': 3000
            }
        }

        # Build reaction network
        params['k'] = define_rates(params)
        params['R'], params['tags'] = build_reactions(params)

        # Initial conditions
        y0 = calculate_initial_densities(params)

        # Time span - shorter for sweep
        t_span = (0, 50)
        t_eval = [50]  # Only final time

        # Create ODE function
        ode_func = PlasmaODE(params)

        # Solve
        sol = solve_ivp(
            ode_func,
            t_span,
            y0,
            method='BDF',
            t_eval=t_eval,
            rtol=1e-6,
            atol=1e-7,
            max_step=1.0
        )

        if not sol.success:
            return None

        # Extract results
        species = params['species']
        results = {}
        for name in ['H', 'CH', 'C2']:
            try:
                idx = species.index(name)
                results[name] = sol.y[idx, -1]
            except ValueError:
                results[name] = 1e-30

        # Calculate error metric (log-scale, more appropriate for plasmas)
        error_h = np.abs(np.log10(results['H'] / TARGETS['H']))
        error_ch = np.abs(np.log10(results['CH'] / TARGETS['CH']))
        error_c2 = np.abs(np.log10(results['C2'] / TARGETS['C2']))

        total_error = error_h + error_ch + error_c2

        # Calculate ratios
        ratio_h_ch = results['H'] / (results['CH'] + 1e-30)
        ratio_c2_ch = results['C2'] / (results['CH'] + 1e-30)

        return {
            'Te': Te,
            'ne': ne,
            'H': results['H'],
            'CH': results['CH'],
            'C2': results['C2'],
            'error': total_error,
            'error_h': error_h,
            'error_ch': error_ch,
            'error_c2': error_c2,
            'ratio_h_ch': ratio_h_ch,
            'ratio_c2_ch': ratio_c2_ch,
            'success': True
        }

    except Exception as e:
        print(f"  Run {run_id}: Te={Te:.1f}, ne={ne:.1e} FAILED: {str(e)[:50]}")
        return None


def main():
    """Run parameter sweep."""
    print("=" * 70)
    print("PARAMETER SWEEP FOR CATHODE GLOW (CG)")
    print("=" * 70)
    print("\nTALIF Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm⁻³")
    print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
    print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")

    # Parameter ranges
    Te_values = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0]  # eV
    ne_values = [1e7, 5e7, 1e8, 5e8, 1e9, 5e9, 1e10]      # cm⁻³

    print(f"\nSweep Parameters:")
    print(f"  Te range: {min(Te_values)} - {max(Te_values)} eV ({len(Te_values)} values)")
    print(f"  ne range: {min(ne_values):.0e} - {max(ne_values):.0e} cm⁻³ ({len(ne_values)} values)")
    print(f"  Total runs: {len(Te_values) * len(ne_values)}")

    # Create parameter combinations
    param_combinations = [(Te, ne, i) for i, (Te, ne) in
                         enumerate(itertools.product(Te_values, ne_values))]

    print("\nStarting parameter sweep...")
    print("(This may take several minutes)")

    start_time = time.time()

    # Run simulations (sequential for simplicity)
    results = []
    total = len(param_combinations)

    for i, params in enumerate(param_combinations):
        if (i + 1) % 10 == 0 or i == 0:
            print(f"  Progress: {i+1}/{total} ({100*(i+1)/total:.0f}%)")

        result = run_single_simulation(params)
        if result is not None:
            results.append(result)

    elapsed = time.time() - start_time
    print(f"\nCompleted {len(results)}/{total} successful runs in {elapsed:.1f} seconds")

    if len(results) == 0:
        print("ERROR: No successful runs!")
        return

    # Find best parameters
    results_sorted = sorted(results, key=lambda x: x['error'])
    best = results_sorted[0]

    print("\n" + "=" * 70)
    print("BEST FIT PARAMETERS")
    print("=" * 70)
    print(f"Te:  {best['Te']:.2f} eV")
    print(f"ne:  {best['ne']:.2e} cm⁻³")
    print(f"\nResults:")
    print(f"  H:  {best['H']:.2e} cm⁻³  (target: {TARGETS['H']:.2e}, factor: {best['H']/TARGETS['H']:.2f}x)")
    print(f"  CH: {best['CH']:.2e} cm⁻³  (target: {TARGETS['CH']:.2e}, factor: {best['CH']/TARGETS['CH']:.2f}x)")
    print(f"  C2: {best['C2']:.2e} cm⁻³  (target: {TARGETS['C2']:.2e}, factor: {best['C2']/TARGETS['C2']:.2f}x)")
    print(f"\nTotal error: {best['error']:.3f} (log-scale)")

    # Show top 5
    print("\n" + "=" * 70)
    print("TOP 5 PARAMETER SETS")
    print("=" * 70)
    print(f"{'Rank':<5} {'Te (eV)':<10} {'ne (cm⁻³)':<12} {'H factor':<10} {'CH factor':<10} {'C2 factor':<10} {'Error':<8}")
    print("-" * 70)

    for i, r in enumerate(results_sorted[:5]):
        h_factor = r['H'] / TARGETS['H']
        ch_factor = r['CH'] / TARGETS['CH']
        c2_factor = r['C2'] / TARGETS['C2']
        print(f"{i+1:<5} {r['Te']:<10.1f} {r['ne']:<12.1e} {h_factor:<10.2f} {ch_factor:<10.2f} {c2_factor:<10.2f} {r['error']:<8.3f}")

    # Create visualization
    print("\nGenerating parameter sweep plots...")
    plot_sweep_results(results, Te_values, ne_values)

    # Save results
    print("\nSaving results to parameter_sweep_results.txt...")
    with open('parameter_sweep_results.txt', 'w') as f:
        f.write("Parameter Sweep Results for CG Region\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"BEST FIT:\n")
        f.write(f"  Te = {best['Te']:.2f} eV\n")
        f.write(f"  ne = {best['ne']:.2e} cm⁻³\n")
        f.write(f"\nRESULTS:\n")
        f.write(f"  H:  {best['H']:.2e} cm⁻³ (factor: {best['H']/TARGETS['H']:.2f}x)\n")
        f.write(f"  CH: {best['CH']:.2e} cm⁻³ (factor: {best['CH']/TARGETS['CH']:.2f}x)\n")
        f.write(f"  C2: {best['C2']:.2e} cm⁻³ (factor: {best['C2']/TARGETS['C2']:.2f}x)\n")
        f.write(f"\nAll results sorted by error:\n")
        f.write(f"{'Te (eV)':<10} {'ne (cm⁻³)':<12} {'H':<12} {'CH':<12} {'C2':<12} {'Error':<8}\n")
        for r in results_sorted:
            f.write(f"{r['Te']:<10.1f} {r['ne']:<12.1e} {r['H']:<12.2e} {r['CH']:<12.2e} {r['C2']:<12.2e} {r['error']:<8.3f}\n")

    print("\nDone! Check:")
    print("  - parameter_sweep_results.txt")
    print("  - parameter_sweep_heatmap.png")


def plot_sweep_results(results, Te_values, ne_values):
    """Create heatmap visualizations."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Create grids for heatmaps
    Te_grid, ne_grid = np.meshgrid(Te_values, ne_values)
    error_grid = np.full_like(Te_grid, np.nan)

    # Fill error grid
    for r in results:
        i_te = Te_values.index(r['Te'])
        i_ne = ne_values.index(r['ne'])
        error_grid[i_ne, i_te] = r['error']

    # Total error heatmap
    ax = axes[0, 0]
    im = ax.pcolormesh(Te_grid, ne_grid, error_grid, shading='auto', cmap='viridis_r')
    ax.set_xlabel('Te (eV)', fontsize=12)
    ax.set_ylabel('ne (cm⁻³)', fontsize=12)
    ax.set_yscale('log')
    ax.set_title('Total Error (lower is better)', fontsize=13, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Log-scale error')
    ax.grid(True, alpha=0.3)

    # H error
    ax = axes[0, 1]
    h_error_grid = np.full_like(Te_grid, np.nan)
    for r in results:
        i_te = Te_values.index(r['Te'])
        i_ne = ne_values.index(r['ne'])
        h_error_grid[i_ne, i_te] = r['error_h']
    im = ax.pcolormesh(Te_grid, ne_grid, h_error_grid, shading='auto', cmap='Reds')
    ax.set_xlabel('Te (eV)', fontsize=12)
    ax.set_ylabel('ne (cm⁻³)', fontsize=12)
    ax.set_yscale('log')
    ax.set_title('H Error', fontsize=13, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Log error')
    ax.grid(True, alpha=0.3)

    # CH error
    ax = axes[1, 0]
    ch_error_grid = np.full_like(Te_grid, np.nan)
    for r in results:
        i_te = Te_values.index(r['Te'])
        i_ne = ne_values.index(r['ne'])
        ch_error_grid[i_ne, i_te] = r['error_ch']
    im = ax.pcolormesh(Te_grid, ne_grid, ch_error_grid, shading='auto', cmap='Blues')
    ax.set_xlabel('Te (eV)', fontsize=12)
    ax.set_ylabel('ne (cm⁻³)', fontsize=12)
    ax.set_yscale('log')
    ax.set_title('CH Error', fontsize=13, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Log error')
    ax.grid(True, alpha=0.3)

    # C2 error
    ax = axes[1, 1]
    c2_error_grid = np.full_like(Te_grid, np.nan)
    for r in results:
        i_te = Te_values.index(r['Te'])
        i_ne = ne_values.index(r['ne'])
        c2_error_grid[i_ne, i_te] = r['error_c2']
    im = ax.pcolormesh(Te_grid, ne_grid, c2_error_grid, shading='auto', cmap='Greens')
    ax.set_xlabel('Te (eV)', fontsize=12)
    ax.set_ylabel('ne (cm⁻³)', fontsize=12)
    ax.set_yscale('log')
    ax.set_title('C2 Error', fontsize=13, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Log error')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('parameter_sweep_heatmap.png', dpi=150)
    print("  Saved: parameter_sweep_heatmap.png")


if __name__ == '__main__':
    main()
