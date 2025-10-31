#!/usr/bin/env python3
"""
ULTRA-REDUCED Parallel CSB Parameter Sweep with H Diffusion Influx

Grid: 12 runs total (2×1×2×3)
Expected time: 30-60 minutes with 12 parallel workers

Focus: Find optimal H_diffusion_influx value quickly
"""

import numpy as np
from scipy.integrate import solve_ivp
import time
import itertools
import json
import multiprocessing as mp
from pathlib import Path
from datetime import datetime
from functools import partial

from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from odefun import PlasmaODE

# CSB Targets
TARGETS = {
    'H': 6.35e14,
    'CH': 9.27e8,
    'C2': 5.56e11,
}

class PlasmaODE_CSB(PlasmaODE):
    """PlasmaODE with custom H diffusion influx for CSB."""
    def __init__(self, params, H_diffusion_influx):
        super().__init__(params)
        self.H_drift_gain = H_diffusion_influx


def calculate_initial_densities(params):
    """Calculate initial species densities for CSB."""
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
    n_tot = params['P'] * 9.66e15

    set_density('e', ne)
    set_density('Ar', 0.85 * n_tot)
    set_density('CH4', 0.15 * n_tot)
    set_density('H', 1e14)
    set_density('ArPlus', ne * 0.5)
    set_density('CH4Plus', ne * 0.05)
    set_density('CH3Plus', ne * 0.05)
    set_density('CH5Plus', ne * 0.01)
    set_density('ArHPlus', ne * 0.1)
    set_density('CH3Minus', ne * 0.01)
    set_density('H2', 1e13)
    set_density('ArStar', ne * 0.3)
    set_density('CH', 5e8)
    set_density('C2', 5e11)
    set_density('CH2', 1e11)
    set_density('CH3', 1e8)
    set_density('C2H2', 1e12)
    set_density('C2H4', 5e7)
    set_density('C', 1e8)

    return y0


def run_single_simulation(param_values, param_keys):
    """Run a single simulation with given parameters."""
    params_dict = dict(zip(param_keys, param_values))

    try:
        params = {
            'P': 0.4,
            'n_tot': 9.66e15,
            'L_discharge': 0.45,
            'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                        'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                        'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                        'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                        'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
            'ion_species': ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH3Minus',
                            'H3Plus', 'CHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus',
                            'HMinus', 'C2HPlus', 'H2Plus'],
            'mobilities': {
                'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
                'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
                'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
                'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
            },
            'scale_e_impact': 1.0
        }

        params.update(params_dict)
        H_diffusion_influx = params_dict['H_diffusion_influx']

        params['k'] = define_rates_tunable(params)
        params['R'], params['tags'] = build_reactions(params)

        y0 = calculate_initial_densities(params)

        # Shorter integration for speed
        t_span = (0, 0.5)  # 0.5 seconds only
        t_eval = [0.5]

        t0 = time.time()
        ode = PlasmaODE_CSB(params, H_diffusion_influx)
        sol = solve_ivp(
            ode,
            t_span,
            y0,
            method='BDF',
            t_eval=t_eval,
            rtol=1e-5,
            atol=1e-6,
            max_step=0.1
        )
        runtime = time.time() - t0

        if sol.success:
            species = params['species']
            y_final = sol.y[:, -1]

            def get_density(name):
                try:
                    idx = species.index(name)
                    return y_final[idx]
                except ValueError:
                    return 0.0

            nH = get_density('H')
            nCH = get_density('CH')
            nC2 = get_density('C2')
            nH2 = get_density('H2')
            ne_final = get_density('e')

            err_H = abs(np.log10(nH / TARGETS['H'])) if nH > 0 else 999
            err_CH = abs(np.log10(nCH / TARGETS['CH'])) if nCH > 0 else 999
            err_C2 = abs(np.log10(nC2 / TARGETS['C2'])) if nC2 > 0 else 999
            total_error = err_H + err_CH + err_C2

            return {
                'params': params_dict,
                'H': nH,
                'CH': nCH,
                'C2': nC2,
                'H2': nH2,
                'ne_final': ne_final,
                'error_H': err_H,
                'error_CH': err_CH,
                'error_C2': err_C2,
                'total_error': total_error,
                'runtime': runtime,
                'success': True,
                'message': sol.message,
                'nfev': sol.nfev
            }
        else:
            return {
                'params': params_dict,
                'success': False,
                'message': sol.message,
                'runtime': runtime
            }

    except Exception as e:
        return {
            'params': params_dict,
            'success': False,
            'message': f"Exception: {str(e)}",
            'runtime': 0
        }


def main():
    """Run ultra-reduced parallel CSB sweep."""
    print("="*80)
    print("ULTRA-REDUCED PARALLEL CSB SWEEP (Fast)")
    print("="*80)

    print(f"\nCSB Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm⁻³")
    print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
    print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")

    # ULTRA-REDUCED grid: 12 runs
    param_grid = {
        'ne': [1e9, 5e9],                         # 2 values
        'Te': [1.0],                              # 1 value (fixed)
        'E_field': [100, 200],                    # 2 values
        'Tgas': [400],                            # 1 value (fixed)
        'H_diffusion_influx': [1e18, 1e19, 1e20], # 3 values (KEY!)
        'L_diff': [0.1],                          # 1 value (fixed)
        'gamma_H': [0.01],                        # 1 value (fixed)
    }

    print(f"\nParameter Grid (ULTRA-REDUCED):")
    for param, values in param_grid.items():
        print(f"  {param}: {values}")

    keys = list(param_grid.keys())
    values = list(param_grid.values())
    combinations = list(itertools.product(*values))

    total_runs = len(combinations)
    print(f"\nTotal runs: {total_runs}")

    # Use all but 1 core, or max 12 workers (one per run)
    n_cores = mp.cpu_count()
    n_workers = min(max(1, n_cores - 1), total_runs)
    print(f"CPU cores available: {n_cores}")
    print(f"Workers to use: {n_workers}")

    est_time_per_run = 300  # 5 min estimate
    est_total_parallel = (total_runs / n_workers) * est_time_per_run / 60
    print(f"\nEstimated time: {est_total_parallel:.1f} minutes")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"csb_ultra_reduced_{timestamp}.json"

    run_func = partial(run_single_simulation, param_keys=keys)

    print(f"\nStarting ultra-reduced sweep...")
    print(f"Output: {output_file}")
    print()

    start_time = time.time()
    results = []

    try:
        with mp.Pool(processes=n_workers) as pool:
            for i, result in enumerate(pool.imap_unordered(run_func, combinations), 1):
                results.append(result)

                elapsed = time.time() - start_time
                avg_time = elapsed / i
                eta = avg_time * (total_runs - i)

                status = "✓" if result['success'] else "✗"
                err_str = f"err={result.get('total_error', 999):.2f}" if result['success'] else "FAILED"

                print(f"[{i}/{total_runs}] {status} {err_str} | "
                      f"Runtime: {result.get('runtime', 0):.1f}s | "
                      f"ETA: {eta/60:.1f} min")

                # Save after each completion
                save_data = {
                    'timestamp': timestamp,
                    'targets': TARGETS,
                    'param_grid': param_grid,
                    'total_runs': total_runs,
                    'completed_runs': len(results),
                    'results': results,
                }
                with open(output_file, 'w') as f:
                    json.dump(save_data, f, indent=2)

    except KeyboardInterrupt:
        print("\n\n⚠️  Interrupted. Saving partial results...")

    total_time = time.time() - start_time

    print(f"\n{'='*80}")
    print("SWEEP COMPLETE!")
    print(f"{'='*80}")
    print(f"Total runtime: {total_time/60:.1f} minutes")
    print(f"Successful: {sum(1 for r in results if r['success'])}/{len(results)}")
    print(f"Results: {output_file}")

    # Find best result
    successful = [r for r in results if r['success']]
    if successful:
        best = min(successful, key=lambda x: x['total_error'])
        print(f"\n{'='*80}")
        print("BEST FIT:")
        print(f"{'='*80}")
        print(f"Parameters:")
        for k, v in best['params'].items():
            print(f"  {k}: {v}")
        print(f"\nDensities:")
        print(f"  H:  {best['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {best['H']/TARGETS['H']:.2f})")
        print(f"  CH: {best['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {best['CH']/TARGETS['CH']:.2f})")
        print(f"  C2: {best['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {best['C2']/TARGETS['C2']:.2f})")
        print(f"\nTotal error: {best['total_error']:.3f}")

if __name__ == '__main__':
    mp.freeze_support()
    main()
