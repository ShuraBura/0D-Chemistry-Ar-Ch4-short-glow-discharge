#!/usr/bin/env python3
"""
Single CSB Simulation - Overnight Test
Purpose: Get accurate timing for one full simulation
"""

import numpy as np
from scipy.integrate import solve_ivp
import time
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
print("SINGLE CSB SIMULATION - OVERNIGHT TEST")
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)

# Test parameters
params_dict = {
    'ne': 1e9,
    'Te': 1.0,
    'E_field': 100,
    'Tgas': 400,
    'H_diffusion_influx': 1e19,  # Medium estimate
    'L_diff': 0.1,
    'gamma_H': 0.01,
}

print(f"\nTest Parameters:")
for k, v in params_dict.items():
    print(f"  {k}: {v}")

print(f"\nTargets:")
print(f"  H:  {TARGETS['H']:.2e} cm⁻³")
print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")

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

print("\n[1/4] Defining rates...")
t0 = time.time()
params['k'] = define_rates_tunable(params)
t_rates = time.time() - t0
print(f"  Done in {t_rates:.2f}s")

print("[2/4] Building reactions...")
t0 = time.time()
params['R'], params['tags'] = build_reactions(params)
t_build = time.time() - t0
print(f"  Done in {t_build:.2f}s")
print(f"  Species: {len(params['species'])}, Reactions: {len(params['R'])}")

print("[3/4] Setting initial conditions...")
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
print("  Initial conditions set")

print("\n[4/4] Running ODE integration...")
print("  Integration time: 0 to 1.0 seconds")
print("  Method: BDF (stiff)")
print("  This may take 30+ minutes...")
print(f"  Started: {datetime.now().strftime('%H:%M:%S')}")

t0 = time.time()
ode = PlasmaODE_CSB(params, H_diffusion_influx)

try:
    sol = solve_ivp(
        ode,
        (0, 1.0),
        y0,
        method='BDF',
        t_eval=[1.0],
        rtol=1e-5,
        atol=1e-6,
        max_step=0.1
    )
    runtime = time.time() - t0

    print(f"\n{'='*80}")
    print(f"INTEGRATION COMPLETE!")
    print(f"{'='*80}")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Runtime: {runtime:.1f}s ({runtime/60:.1f} minutes)")
    print(f"Success: {sol.success}")
    print(f"Message: {sol.message}")
    print(f"Function evaluations: {sol.nfev}")

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
        nH2 = get_density('H2')

        print(f"\n{'='*80}")
        print("RESULTS")
        print(f"{'='*80}")
        print(f"{'Species':<10} {'Model':<15} {'Target':<15} {'Ratio':<10}")
        print(f"{'-'*60}")
        print(f"{'H':<10} {nH:<15.2e} {TARGETS['H']:<15.2e} {nH/TARGETS['H']:<10.2f}")
        print(f"{'CH':<10} {nCH:<15.2e} {TARGETS['CH']:<15.2e} {nCH/TARGETS['CH']:<10.2f}")
        print(f"{'C2':<10} {nC2:<15.2e} {TARGETS['C2']:<15.2e} {nC2/TARGETS['C2']:<10.2f}")

        err_H = abs(np.log10(nH / TARGETS['H'])) if nH > 0 else 999
        err_CH = abs(np.log10(nCH / TARGETS['CH'])) if nCH > 0 else 999
        err_C2 = abs(np.log10(nC2 / TARGETS['C2'])) if nC2 > 0 else 999
        total_error = err_H + err_CH + err_C2

        print(f"\nLog errors:")
        print(f"  H:  {err_H:.3f}")
        print(f"  CH: {err_CH:.3f}")
        print(f"  C2: {err_C2:.3f}")
        print(f"  Total: {total_error:.3f}")

        # Save results
        result = {
            'timestamp': datetime.now().isoformat(),
            'params': params_dict,
            'targets': TARGETS,
            'results': {
                'H': nH,
                'CH': nCH,
                'C2': nC2,
                'H2': nH2,
            },
            'errors': {
                'H': err_H,
                'CH': err_CH,
                'C2': err_C2,
                'total': total_error
            },
            'timing': {
                'rates_s': t_rates,
                'build_s': t_build,
                'integration_s': runtime,
                'total_s': t_rates + t_build + runtime
            },
            'success': True,
            'message': sol.message,
            'nfev': sol.nfev
        }

        with open('single_test_result.json', 'w') as f:
            json.dump(result, f, indent=2)

        print(f"\n✓ Results saved to: single_test_result.json")

        print(f"\n{'='*80}")
        print("TIMING ESTIMATE FOR FULL SWEEP")
        print(f"{'='*80}")
        print(f"Single run: {runtime/60:.1f} minutes")
        print(f"12 runs:    {runtime*12/3600:.1f} hours")
        print(f"36 runs:    {runtime*36/3600:.1f} hours")

    else:
        print(f"\n✗ Integration FAILED: {sol.message}")

        result = {
            'timestamp': datetime.now().isoformat(),
            'params': params_dict,
            'success': False,
            'message': sol.message,
            'runtime': runtime
        }

        with open('single_test_result.json', 'w') as f:
            json.dump(result, f, indent=2)

except Exception as e:
    runtime = time.time() - t0
    print(f"\n✗ EXCEPTION after {runtime/60:.1f} minutes:")
    print(f"  {e}")

    result = {
        'timestamp': datetime.now().isoformat(),
        'params': params_dict,
        'success': False,
        'exception': str(e),
        'runtime': runtime
    }

    with open('single_test_result.json', 'w') as f:
        json.dump(result, f, indent=2)

print(f"\n{'='*80}")
print("TEST COMPLETE")
print(f"{'='*80}")
