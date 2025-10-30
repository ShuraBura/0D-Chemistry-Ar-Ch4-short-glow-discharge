#!/usr/bin/env python3
"""
CSB Validation with FIXED H density
H is diffused from CG, not produced locally in CSB
This should be MUCH faster since H is not evolved dynamically

User-specified CSB parameters:
- E-field: 20-300 V/cm
- Te: 0.5-2 eV
- ne: 9e8-9e9 cm⁻³
- H: FIXED (diffused from CG)
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions

# CSB Targets
TARGETS = {
    'H': 6.35e14,    # cm⁻³ (FIXED - diffused from CG)
    'CH': 9.27e8,    # cm⁻³ (evolved)
    'C2': 5.56e11,   # cm⁻³ (evolved)
}

print("="*80)
print("CSB VALIDATION with FIXED H (Option B)")
print("="*80)
print(f"\nCSB Targets:")
print(f"  H:  {TARGETS['H']:.2e} cm⁻³  [FIXED - diffused from CG]")
print(f"  CH: {TARGETS['CH']:.2e} cm⁻³  [evolved]")
print(f"  C2: {TARGETS['C2']:.2e} cm⁻³  [evolved]")

# CSB parameters from user
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
    'ne': 5e9,          # cm⁻³ (user range: 9e8-9e9)
    'Te': 1.0,          # eV (user range: 0.5-2 eV)
    'E_field': 100,     # V/cm (user range: 20-300 V/cm)
    'Tgas': 400,        # K
    'L_diff': 0.1,      # cm
    'gamma_H': 0.01,
    'scale_e_impact': 1.0,
    'nH_fixed': TARGETS['H'],  # FIXED H density
}

print(f"\nCSB Parameters (User-specified ranges):")
print(f"  ne = {params['ne']:.1e} cm⁻³  (range: 9e8-9e9)")
print(f"  Te = {params['Te']} eV  (range: 0.5-2 eV)")
print(f"  E = {params['E_field']} V/cm  (range: 20-300 V/cm)")
print(f"  Tg = {params['Tgas']} K")
print(f"  H = {params['nH_fixed']:.2e} cm⁻³  [FIXED - from diffusion]")

print("\nDefining rates...")
t0 = time.time()
params['k'] = define_rates_tunable(params)
print(f"  Done in {time.time()-t0:.2f}s")

print("Building reactions...")
t0 = time.time()
params['R'], params['tags'] = build_reactions(params)
print(f"  Done in {time.time()-t0:.2f}s")
print(f"  Species: {len(params['species'])}, Reactions: {len(params['R'])}")

species = params['species']

# ODE function with FIXED H
class PlasmaODE_FixedH:
    def __init__(self, params):
        self.params = params
        self.species = params['species']
        self.R = params['R']
        self.k = params['k']
        self.n_tot = params['n_tot']
        self.nH_fixed = params['nH_fixed']

        # Find H index
        try:
            self.H_idx = self.species.index('H')
        except ValueError:
            self.H_idx = None

    def __call__(self, t, y):
        # Compute rates
        rates = np.zeros(len(self.species))
        for reaction in self.R:
            rate = self.k.get(reaction['k_name'], 0.0)
            for reactant in reaction['reactants']:
                idx = self.species.index(reactant)
                rate *= y[idx]

            for reactant in reaction['reactants']:
                idx = self.species.index(reactant)
                rates[idx] -= rate

            for product in reaction['products']:
                idx = self.species.index(product)
                rates[idx] += rate

        # SET H rate to zero (FIXED)
        if self.H_idx is not None:
            rates[self.H_idx] = 0.0

        return rates

# Initial densities
ns = len(species)
y0 = np.ones(ns) * 1e3

def set_density(name, value):
    try:
        idx = species.index(name)
        y0[idx] = value
    except ValueError:
        pass

n_tot = params['P'] * 9.66e15
n_Ar = 0.85 * n_tot
n_CH4 = 0.15 * n_tot

set_density('e', params['ne'])
set_density('Ar', n_Ar)
set_density('CH4', n_CH4)
set_density('H', TARGETS['H'])  # FIXED
set_density('ArPlus', params['ne'] * 0.5)
set_density('H2', 1e13)
set_density('CH', 5e8)
set_density('C2', 5e11)
set_density('CH3', 1e8)

print(f"\nInitial conditions (H FIXED at {TARGETS['H']:.2e}).")

# Time span
t_span = (0, 5)  # Even shorter - 5s
t_eval = [5]

print(f"\nIntegrating with FIXED H...")
print(f"  Time span: 0 to {t_span[1]:.1f} s")

t0 = time.time()
ode = PlasmaODE_FixedH(params)
sol = solve_ivp(
    ode,
    t_span,
    y0,
    method='BDF',
    t_eval=t_eval,
    rtol=1e-5,
    atol=1e-6,
    max_step=0.5
)
runtime = time.time() - t0

print(f"\n{'='*80}")
if sol.success:
    print(f"✓ SUCCESS! Runtime: {runtime:.2f}s")

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

    print(f"\n{'='*80}")
    print("RESULTS (H FIXED)")
    print(f"{'='*80}")
    print(f"{'Species':<10} {'Model':<15} {'Target':<15} {'Ratio':<10}")
    print(f"{'-'*80}")
    print(f"{'H':<10} {nH:<15.2e} {TARGETS['H']:<15.2e} {nH/TARGETS['H']:<10.2f}  [FIXED]")
    print(f"{'CH':<10} {nCH:<15.2e} {TARGETS['CH']:<15.2e} {nCH/TARGETS['CH']:<10.2f}")
    print(f"{'C2':<10} {nC2:<15.2e} {TARGETS['C2']:<15.2e} {nC2/TARGETS['C2']:<10.2f}")

    # Errors for CH and C2 only
    err_CH = abs(np.log10(nCH / TARGETS['CH'])) if nCH > 0 else 999
    err_C2 = abs(np.log10(nC2 / TARGETS['C2'])) if nC2 > 0 else 999
    total_error = err_CH + err_C2

    print(f"\nLog errors (CH + C2 only):")
    print(f"  CH: {err_CH:.3f}")
    print(f"  C2: {err_C2:.3f}")
    print(f"  Total: {total_error:.3f}")

    if total_error < 0.5:
        print("\n  🎉 EXCELLENT! CH and C₂ within factor 3!")
    elif total_error < 1.0:
        print("\n  ✓ GOOD! CH and C₂ within factor 10!")
    else:
        print("\n  ⚠️  Needs parameter tuning")

else:
    print(f"✗ FAILED! Runtime: {runtime:.2f}s")
    print(f"  {sol.message}")

print(f"\n{'='*80}")
print(f"Runtime with FIXED H: {runtime:.2f}s")
print(f"\nNext step: Run parameter sweep with FIXED H for CSB")
print(f"  Sweep ranges:")
print(f"    ne: [1e9, 3e9, 5e9, 9e9] cm⁻³")
print(f"    Te: [0.5, 1.0, 1.5, 2.0] eV")
print(f"    E:  [20, 50, 100, 200, 300] V/cm")
print(f"    Tg: [300, 400] K")
print(f"{'='*80}")
