#!/usr/bin/env python3
"""
SEQUENTIAL Ultra-Reduced CSB Sweep (No Multiprocessing)

Runs 12 simulations one at a time to avoid multiprocessing issues.
Expected time: ~1-2 hours (but will actually complete!)
"""

import numpy as np
from scipy.integrate import solve_ivp
import time
import itertools
import json
from datetime import datetime

from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from odefun import PlasmaODE

TARGETS = {
    'H': 6.35e14,
    'CH': 9.27e8,
    'C2': 5.56e11,
}

class PlasmaODE_CSB(PlasmaODE):
    def __init__(self, params, H_diffusion_influx):
        super().__init__(params)
        self.H_drift_gain = H_diffusion_influx

print("="*80)
print("SEQUENTIAL ULTRA-REDUCED CSB SWEEP (No Parallel)")
print("="*80)
print(f"\nCSB Targets:")
print(f"  H:  {TARGETS['H']:.2e} cm⁻³")
print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")

# Ultra-reduced grid
param_grid = {
    'ne': [1e9, 5e9],
    'Te': [1.0],
    'E_field': [100, 200],
    'Tgas': [400],
    'H_diffusion_influx': [1e18, 1e19, 1e20],
    'L_diff': [0.1],
    'gamma_H': [0.01],
}

print(f"\nParameter Grid:")
for param, values in param_grid.items():
    print(f"  {param}: {values}")

keys = list(param_grid.keys())
values = list(param_grid.values())
combinations = list(itertools.product(*values))

total_runs = len(combinations)
print(f"\nTotal runs: {total_runs}")
print(f"Running SEQUENTIALLY (no parallel)")
print(f"Estimated time: ~60-120 minutes\n")

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_file = f"csb_sequential_{timestamp}.json"

start_time = time.time()
results = []

for i, param_values in enumerate(combinations, 1):
    params_dict = dict(zip(keys, param_values))

    print(f"\n[{i}/{total_runs}] Running: ne={params_dict['ne']:.1e}, E={params_dict['E_field']}, H_flux={params_dict['H_diffusion_influx']:.1e}")

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

        # Run simulation
        t0 = time.time()
        ode = PlasmaODE_CSB(params, H_diffusion_influx)
        sol = solve_ivp(ode, (0, 0.5), y0, method='BDF', t_eval=[0.5], rtol=1e-5, atol=1e-6, max_step=0.1)
        runtime = time.time() - t0

        if sol.success:
            def get_density(name):
                try:
                    idx = species.index(name)
                    return sol.y[idx, -1]
                except ValueError:
                    return 0.0

            nH = get_density('H')
            nCH = get_density('CH')
            nC2 = get_density('C2')

            err_H = abs(np.log10(nH / TARGETS['H'])) if nH > 0 else 999
            err_CH = abs(np.log10(nCH / TARGETS['CH'])) if nCH > 0 else 999
            err_C2 = abs(np.log10(nC2 / TARGETS['C2'])) if nC2 > 0 else 999
            total_error = err_H + err_CH + err_C2

            result = {
                'params': params_dict,
                'H': nH,
                'CH': nCH,
                'C2': nC2,
                'error_H': err_H,
                'error_CH': err_CH,
                'error_C2': err_C2,
                'total_error': total_error,
                'runtime': runtime,
                'success': True
            }

            print(f"  ✓ SUCCESS! Runtime: {runtime:.1f}s")
            print(f"    H={nH:.2e} CH={nCH:.2e} C2={nC2:.2e}")
            print(f"    Error: {total_error:.3f}")

        else:
            result = {
                'params': params_dict,
                'success': False,
                'message': sol.message,
                'runtime': runtime
            }
            print(f"  ✗ FAILED: {sol.message}")

    except Exception as e:
        result = {
            'params': params_dict,
            'success': False,
            'message': str(e),
            'runtime': 0
        }
        print(f"  ✗ EXCEPTION: {e}")

    results.append(result)

    # Save after each run
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

    # Progress
    elapsed = time.time() - start_time
    avg_time = elapsed / i
    eta = avg_time * (total_runs - i)
    print(f"  Progress: {i}/{total_runs} | Avg: {avg_time/60:.1f} min/run | ETA: {eta/60:.1f} min")

total_time = time.time() - start_time

print(f"\n{'='*80}")
print("SWEEP COMPLETE!")
print(f"{'='*80}")
print(f"Total runtime: {total_time/60:.1f} minutes")
print(f"Successful: {sum(1 for r in results if r['success'])}/{len(results)}")
print(f"Results: {output_file}")

# Find best
successful = [r for r in results if r['success']]
if successful:
    best = min(successful, key=lambda x: x['total_error'])
    print(f"\n{'='*80}")
    print("BEST FIT:")
    print(f"{'='*80}")
    for k, v in best['params'].items():
        print(f"  {k}: {v}")
    print(f"\nDensities:")
    print(f"  H:  {best['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {best['H']/TARGETS['H']:.2f})")
    print(f"  CH: {best['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {best['CH']/TARGETS['CH']:.2f})")
    print(f"  C2: {best['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {best['C2']/TARGETS['C2']:.2f})")
    print(f"\nTotal error: {best['total_error']:.3f}")
