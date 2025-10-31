#!/usr/bin/env python3
"""
Test single CSB simulation to diagnose multiprocessing issues
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from odefun import PlasmaODE

# PlasmaODE with custom H diffusion influx
class PlasmaODE_CSB(PlasmaODE):
    def __init__(self, params, H_diffusion_influx):
        super().__init__(params)
        self.H_drift_gain = H_diffusion_influx

print("Testing single CSB simulation...")

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
    'ne': 1e9,
    'Te': 0.8,
    'E_field': 50,
    'Tgas': 400,
    'L_diff': 0.1,
    'gamma_H': 0.01,
    'scale_e_impact': 1.0,
}

H_diffusion_influx = 1e19

print(f"Parameters: ne={params['ne']:.1e}, Te={params['Te']}, E={params['E_field']}, H_flux={H_diffusion_influx:.1e}")

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

n_tot = params['P'] * 9.66e15
set_density('e', params['ne'])
set_density('Ar', 0.85 * n_tot)
set_density('CH4', 0.15 * n_tot)
set_density('H', 1e14)
set_density('ArPlus', params['ne'] * 0.5)
set_density('CH', 5e8)
set_density('C2', 5e11)

print("Running integration (0.5 sec)...")
t0 = time.time()
ode = PlasmaODE_CSB(params, H_diffusion_influx)
sol = solve_ivp(ode, (0, 0.5), y0, method='BDF', t_eval=[0.5], rtol=1e-5, atol=1e-6, max_step=0.1)
runtime = time.time() - t0

print(f"Runtime: {runtime:.1f}s")
print(f"Success: {sol.success}")

if sol.success:
    idx_H = species.index('H')
    idx_CH = species.index('CH')
    idx_C2 = species.index('C2')
    print(f"H:  {sol.y[idx_H, -1]:.2e}")
    print(f"CH: {sol.y[idx_CH, -1]:.2e}")
    print(f"C2: {sol.y[idx_C2, -1]:.2e}")
    print("\n✓ Single simulation WORKS!")
else:
    print(f"✗ Failed: {sol.message}")
