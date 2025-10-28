#!/usr/bin/env python3
"""
Focused Parameter Sweep for CG with Experimental Constraints

Hard constraints:
- ne ≤ 1e9 cm⁻³
- E ~ 500 V/cm
- CG: 0-2 mm from copper cathode
- Tg: 570-800 K
- L_diff: 0-1 mm
"""

import numpy as np
from scipy.integrate import solve_ivp
import time
import itertools
import json

from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from odefun import PlasmaODE


# TALIF Targets (CG)
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

    ne = params['ne']
    set_density('e', ne)
    set_density('Ar', 0.85 * 9.66e15)
    set_density('CH4', 0.15 * 9.66e15)
    set_density('ArPlus', ne * 0.5)
    set_density('CH4Plus', ne * 0.05)
    set_density('CH3Plus', ne * 0.05)
    set_density('CH5Plus', ne * 0.01)
    set_density('ArHPlus', ne * 0.1)
    set_density('CH3Minus', ne * 0.01)
    set_density('H2', 5e12)  # Higher initial H2
    set_density('ArStar', ne * 0.5)
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


def run_simulation(params_dict):
    """Run a single simulation."""
    try:
        # Base parameters
        params = {
            'P': 0.4,
            'n_tot': 9.66e15,
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

        # Add variable parameters
        params.update(params_dict)

        # Build reactions with tunable rates
        params['k'] = define_rates_tunable(params)
        params['R'], params['tags'] = build_reactions(params)

        # Initial conditions
        y0 = calculate_initial_densities(params)

        # Shorter time for sweep
        t_span = (0, 30)
        t_eval = [30]

        # Solve
        ode_func = PlasmaODE(params)
        sol = solve_ivp(
            ode_func,
            t_span,
            y0,
            method='BDF',
            t_eval=t_eval,
            rtol=1e-5,
            atol=1e-6,
            max_step=0.5
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

        # Calculate error (log-scale)
        error_h = np.abs(np.log10(results['H'] / TARGETS['H']))
        error_ch = np.abs(np.log10(results['CH'] / TARGETS['CH']))
        error_c2 = np.abs(np.log10(results['C2'] / TARGETS['C2']))
        total_error = error_h + error_ch + error_c2

        return {
            **params_dict,
            'H': results['H'],
            'CH': results['CH'],
            'C2': results['C2'],
            'error': total_error,
            'error_h': error_h,
            'error_ch': error_ch,
            'error_c2': error_c2,
            'success': True
        }

    except Exception as e:
        print(f"  ERROR: {str(e)[:80]}")
        return None


def main():
    """Run focused parameter sweep with experimental constraints."""
    print("=" * 80)
    print("COMPREHENSIVE PARAMETER SWEEP - CG with Experimental Constraints")
    print("=" * 80)

    print("\nExperimental Constraints:")
    print("  ne ≤ 1e9 cm⁻³  (cannot go higher in CG)")
    print("  E: 400-800 V/cm (gradient in cathode fall, now TUNABLE!)")
    print("  Cathode: Copper (clean initially, contaminated during operation)")
    print("  Tg: 570-800 K (higher than assumed!)")
    print("  L_diff: 0-1 mm")
    print("  Spatial averaging: 0-1 mm from cathode (CRITICAL!)")

    print("\nTALIF Targets (CG):")
    print(f"  H:  {TARGETS['H']:.2e} cm⁻³")
    print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
    print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")

    # COMPREHENSIVE parameter ranges (with E-field sweep!)
    param_grid = {
        'ne': [3e8, 5e8, 8e8, 1e9],              # cm⁻³ (MAX 1e9!)
        'Te': [3.0, 5.0, 7.0],                   # eV (high, non-thermal)
        'E_field': [400, 500, 600, 800],         # V/cm (TUNABLE, gradient in cathode fall)
        'Tgas': [570, 700],                      # K (from measurement)
        'L_diff': [0.1],                         # cm (fixed at 1 mm)
        'gamma_H': [0.001, 0.01, 0.05, 0.1],     # Copper: clean→contaminated
    }

    print("\nParameter Ranges:")
    for param, values in param_grid.items():
        print(f"  {param}: {values}")

    # Generate all combinations
    keys = list(param_grid.keys())
    values = list(param_grid.values())
    combinations = list(itertools.product(*values))

    total_runs = len(combinations)
    print(f"\nTotal runs: {total_runs}")
    if total_runs > 200:
        print(f"Expected runtime: ~30-60 minutes (comprehensive sweep)")
    elif total_runs > 100:
        print(f"Expected runtime: ~15-30 minutes")
    else:
        print(f"Expected runtime: ~5-15 minutes")

    # Run sweep
    print("\nStarting sweep...")
    start_time = time.time()

    results = []
    for i, combo in enumerate(combinations):
        params_dict = dict(zip(keys, combo))

        if (i + 1) % 5 == 0 or i == 0:
            print(f"  Progress: {i+1}/{total_runs} ({100*(i+1)/total_runs:.0f}%) "
                  f"[Te={params_dict['Te']:.1f}, ne={params_dict['ne']:.1e}, "
                  f"γ_H={params_dict['gamma_H']:.3f}]")

        result = run_simulation(params_dict)
        if result is not None:
            results.append(result)

    elapsed = time.time() - start_time
    print(f"\nCompleted {len(results)}/{total_runs} runs in {elapsed:.1f} s")
    print(f"Average: {elapsed/max(len(results),1):.2f} s/run")

    if len(results) == 0:
        print("ERROR: No successful runs!")
        return

    # Sort by error
    results_sorted = sorted(results, key=lambda x: x['error'])
    best = results_sorted[0]

    # Print results
    print("\n" + "=" * 80)
    print("BEST FIT PARAMETERS")
    print("=" * 80)
    print(f"  ne:      {best['ne']:.2e} cm⁻³")
    print(f"  Te:      {best['Te']:.1f} eV")
    print(f"  E-field: {best['E_field']:.0f} V/cm")
    print(f"  Tgas:    {best['Tgas']:.0f} K")
    print(f"  L_diff:  {best['L_diff']:.3f} cm")
    print(f"  γ_H:     {best['gamma_H']:.4f} (copper recombination)")

    print("\nResults vs Targets:")
    print(f"  {'Species':<8} {'Target':<15} {'Model':<15} {'Factor':<10} {'Log Error':<10}")
    print("  " + "-" * 60)
    for name in ['H', 'CH', 'C2']:
        target = TARGETS[name]
        model = best[name]
        factor = model / target
        log_err = best[f'error_{name.lower()}']
        print(f"  {name:<8} {target:<15.2e} {model:<15.2e} {factor:<10.2f} {log_err:<10.3f}")

    print(f"\nTotal error: {best['error']:.3f}")

    # Check if any species is close
    h_close = 0.1 < (best['H'] / TARGETS['H']) < 10
    ch_close = 0.1 < (best['CH'] / TARGETS['CH']) < 10
    c2_close = 0.1 < (best['C2'] / TARGETS['C2']) < 10

    print("\nAssessment:")
    print(f"  H:  {'✓ Within factor 10' if h_close else '✗ Outside factor 10'}")
    print(f"  CH: {'✓ Within factor 10' if ch_close else '✗ Outside factor 10'}")
    print(f"  C2: {'✓ Within factor 10' if c2_close else '✗ Outside factor 10'}")

    # Show top 10
    print("\n" + "=" * 80)
    print("TOP 10 PARAMETER SETS")
    print("=" * 80)
    print(f"{'#':<3} {'ne':>10} {'Te':>6} {'Tg':>5} {'γ_H':>7} {'H_factor':>9} {'CH_factor':>10} {'C2_factor':>10} {'Error':>7}")
    print("-" * 80)

    for i, r in enumerate(results_sorted[:10]):
        h_factor = r['H'] / TARGETS['H']
        ch_factor = r['CH'] / TARGETS['CH']
        c2_factor = r['C2'] / TARGETS['C2']
        print(f"{i+1:<3} {r['ne']:>10.1e} {r['Te']:>6.1f} {r['Tgas']:>5.0f} "
              f"{r['gamma_H']:>7.4f} {h_factor:>9.2f} {ch_factor:>10.2f} {c2_factor:>10.2f} {r['error']:>7.3f}")

    # Save results
    print("\nSaving results...")
    with open('cg_constrained_sweep_results.json', 'w') as f:
        json.dump({
            'best': best,
            'all_results': results_sorted[:50]  # Top 50
        }, f, indent=2)

    print("\nFiles created:")
    print("  - cg_constrained_sweep_results.json")
    print("\nDone!")


if __name__ == '__main__':
    main()
