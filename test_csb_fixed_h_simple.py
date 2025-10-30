#!/usr/bin/env python3
"""
CSB Validation with FIXED H density - SIMPLE VERSION
Just set H derivative to zero after computing reactions
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from odefun import PlasmaODE

# CSB Targets
TARGETS = {
    'H': 6.35e14,    # cm⁻³ (FIXED)
    'CH': 9.27e8,    # cm⁻³
    'C2': 5.56e11,   # cm⁻³
}

# Fixed H ODE wrapper
class PlasmaODE_FixedH(PlasmaODE):
    """Plasma ODE with FIXED H density."""
    def __init__(self, params, fixed_species=None):
        super().__init__(params)
        self.fixed_species = fixed_species or {}
        self.fixed_indices = {}
        for name in self.fixed_species:
            try:
                self.fixed_indices[name] = self.species.index(name)
            except ValueError:
                pass

    def __call__(self, t, y):
        dydt = super().__call__(t, y)
        # Set derivatives to zero for fixed species
        for name, idx in self.fixed_indices.items():
            dydt[idx] = 0.0
        return dydt

print("="*80)
print("CSB VALIDATION with FIXED H")
print("="*80)
print(f"\nCSB Targets:")
print(f"  H:  {TARGETS['H']:.2e} cm⁻³  [FIXED]")
print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")

# CSB parameters
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
    'ne': 5e9,          # cm⁻³
    'Te': 1.0,          # eV
    'E_field': 100,     # V/cm
    'Tgas': 400,        # K
    'L_diff': 0.1,      # cm
    'gamma_H': 0.01,
    'scale_e_impact': 1.0
}

print(f"\nCSB Parameters:")
print(f"  ne = {params['ne']:.1e} cm⁻³")
print(f"  Te = {params['Te']} eV")
print(f"  E = {params['E_field']} V/cm")
print(f"  H = {TARGETS['H']:.2e} cm⁻³ [FIXED]")

print("\nDefining rates...")
t0 = time.time()
params['k'] = define_rates_tunable(params)
print(f"  Done in {time.time()-t0:.2f}s")

print("Building reactions...")
t0 = time.time()
params['R'], params['tags'] = build_reactions(params)
print(f"  Done in {time.time()-t0:.2f}s")

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
set_density('H', TARGETS['H'])  # FIXED
set_density('ArPlus', params['ne'] * 0.5)
set_density('CH', 5e8)
set_density('C2', 5e11)

print(f"\nIntegrating (H FIXED)...")
t0 = time.time()
ode = PlasmaODE_FixedH(params, fixed_species={'H'})
sol = solve_ivp(
    ode,
    (0, 5),
    y0,
    method='BDF',
    t_eval=[5],
    rtol=1e-5,
    atol=1e-6,
    max_step=0.5
)
runtime = time.time() - t0

print(f"\n{'='*80}")
if sol.success:
    print(f"✓ SUCCESS! Runtime: {runtime:.2f}s")

    def get_density(name):
        try:
            idx = species.index(name)
            return sol.y[idx, -1]
        except ValueError:
            return 0.0

    nH = get_density('H')
    nCH = get_density('CH')
    nC2 = get_density('C2')

    print(f"\n{'Species':<10} {'Model':<15} {'Target':<15} {'Ratio':<10}")
    print(f"{'-'*60}")
    print(f"{'H (FIXED)':<10} {nH:<15.2e} {TARGETS['H']:<15.2e} {nH/TARGETS['H']:<10.2f}")
    print(f"{'CH':<10} {nCH:<15.2e} {TARGETS['CH']:<15.2e} {nCH/TARGETS['CH']:<10.2f}")
    print(f"{'C2':<10} {nC2:<15.2e} {TARGETS['C2']:<15.2e} {nC2/TARGETS['C2']:<10.2f}")

    err_CH = abs(np.log10(nCH / TARGETS['CH'])) if nCH > 0 else 999
    err_C2 = abs(np.log10(nC2 / TARGETS['C2'])) if nC2 > 0 else 999

    print(f"\nErrors (log10):")
    print(f"  CH: {err_CH:.3f}")
    print(f"  C2: {err_C2:.3f}")
    print(f"  Total: {err_CH + err_C2:.3f}")

    if err_CH + err_C2 < 0.5:
        print("\n🎉 EXCELLENT!")
    elif err_CH + err_C2 < 1.0:
        print("\n✓ GOOD!")

else:
    print(f"✗ FAILED: {sol.message}")

print(f"\n{'='*80}")
print(f"Ready for CSB parameter sweep with FIXED H!")
print(f"Sweep ranges:")
print(f"  ne: [1e9, 3e9, 5e9, 9e9] cm⁻³")
print(f"  Te: [0.5, 1.0, 1.5, 2.0] eV")
print(f"  E:  [20, 50, 100, 200, 300] V/cm")
print(f"{'='*80}")
