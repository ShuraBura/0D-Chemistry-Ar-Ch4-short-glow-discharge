#!/usr/bin/env python3
"""
Test single simulation run to check timing
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from odefun import PlasmaODE

# Updated CG Targets
TARGETS = {
    'H': 8.58e15,    # cm⁻³
    'CH': 4.6e8,     # cm⁻³
    'C2': 1.44e11,   # cm⁻³
}

print("="*80)
print("SINGLE SIMULATION TEST - CG Region")
print("="*80)
print(f"\nTargets:")
print(f"  H:  {TARGETS['H']:.2e} cm⁻³")
print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")

# Test parameters
params = {
    'P': 0.4,           # Torr
    'n_tot': 9.66e15,   # cm⁻³
    'L_discharge': 0.45,# cm
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
    'ne': 5e8,          # cm⁻³
    'Te': 3.0,          # eV
    'E_field': 400,     # V/cm
    'Tgas': 570,        # K
    'L_diff': 0.1,      # cm
    'gamma_H': 0.01,    # wall recombination
    'scale_e_impact': 1.0
}

print(f"\nTest Parameters:")
print(f"  ne = {params['ne']:.1e} cm⁻³")
print(f"  Te = {params['Te']} eV")
print(f"  E = {params['E_field']} V/cm")
print(f"  Tg = {params['Tgas']} K")
print(f"  γ_H = {params['gamma_H']}")

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

# Initial densities
ns = len(species)
y0 = np.ones(ns) * 1e3

def set_density(name, value):
    try:
        idx = species.index(name)
        y0[idx] = value
    except ValueError:
        pass

n_tot = params['P'] * 9.66e15  # cm⁻³
n_Ar = 0.85 * n_tot
n_CH4 = 0.15 * n_tot

set_density('e', params['ne'])
set_density('Ar', n_Ar)
set_density('CH4', n_CH4)

print(f"\nInitial conditions set:")
print(f"  ne = {params['ne']:.1e} cm⁻³")
print(f"  n_Ar = {n_Ar:.2e} cm⁻³")
print(f"  n_CH4 = {n_CH4:.2e} cm⁻³")

# Time span
t_span = (0, 30)  # 0 to 30 s
t_eval = [30]  # Only final state

print(f"\nIntegrating ODE system...")
print(f"  Time span: 0 to {t_span[1]:.1f} s")
print(f"  Method: BDF (stiff), max_step=0.5")

t0 = time.time()
ode = PlasmaODE(params)
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
    print(f"  Status: {sol.message}")
    print(f"  Steps: {sol.nfev}")

    # Extract final densities
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

    print(f"\nFinal densities:")
    print(f"  H:  {nH:.2e} cm⁻³  (target: {TARGETS['H']:.2e})")
    print(f"  CH: {nCH:.2e} cm⁻³  (target: {TARGETS['CH']:.2e})")
    print(f"  C2: {nC2:.2e} cm⁻³  (target: {TARGETS['C2']:.2e})")

    # Calculate errors
    err_H = abs(np.log10(nH / TARGETS['H'])) if nH > 0 else 999
    err_CH = abs(np.log10(nCH / TARGETS['CH'])) if nCH > 0 else 999
    err_C2 = abs(np.log10(nC2 / TARGETS['C2'])) if nC2 > 0 else 999
    total_error = err_H + err_CH + err_C2

    print(f"\nLog errors:")
    print(f"  H:  {err_H:.3f}")
    print(f"  CH: {err_CH:.3f}")
    print(f"  C2: {err_C2:.3f}")
    print(f"  Total: {total_error:.3f}")

else:
    print(f"✗ FAILED! Runtime: {runtime:.2f}s")
    print(f"  Status: {sol.message}")

print(f"{'='*80}")
print(f"\nEstimated time for 108 runs: {runtime*108/60:.1f} minutes")
print(f"Estimated time for 384 runs: {runtime*384/60:.1f} minutes")
