#!/usr/bin/env python3
"""
Quick validation test with default L_diff = 0.1 cm
Tests best-guess parameters against measured targets
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from odefun import PlasmaODE

# Updated CG Targets (Measured 2025-10-30)
TARGETS = {
    'H': 9.58e15,    # cm⁻³ (measured spatial average)
    'CH': 4.29e8,    # cm⁻³ (measured spatial average)
    'C2': 1.50e11,   # cm⁻³ (measured spatial average)
}

print("="*80)
print("QUICK VALIDATION TEST - CG Region")
print("="*80)
print(f"\nTargets (measured spatial avg 0-1mm):")
print(f"  H:  {TARGETS['H']:.2e} cm⁻³")
print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")

# Test parameters - using DEFAULT L_diff for speed
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
    'ne': 8e8,          # cm⁻³ (middle-high range)
    'Te': 5.0,          # eV (middle range for non-thermal)
    'E_field': 600,     # V/cm (middle of 400-800 range)
    'Tgas': 570,        # K (measured)
    'L_diff': 0.1,      # cm (DEFAULT for speed - measured was 0.057)
    'gamma_H': 0.01,    # wall recombination (middle range)
    'scale_e_impact': 1.0
}

print(f"\nTest Parameters (Best Guess):")
print(f"  ne = {params['ne']:.1e} cm⁻³")
print(f"  Te = {params['Te']} eV")
print(f"  E = {params['E_field']} V/cm")
print(f"  Tg = {params['Tgas']} K")
print(f"  L_diff = {params['L_diff']} cm (using DEFAULT for speed)")
print(f"  γ_H = {params['gamma_H']}")
print(f"\nNote: Measured L_diff_H = 0.057 cm, but using 0.1 cm for faster integration")

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

# Add reasonable initial guesses for key species
set_density('ArPlus', params['ne'] * 0.5)
set_density('H', 1e15)
set_density('H2', 5e12)
set_density('CH', 1e8)
set_density('C2', 1e11)
set_density('CH3', 5e7)
set_density('CH4Plus', params['ne'] * 0.05)

print(f"\nInitial conditions set.")

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
    nH2 = get_density('H2')
    ne_final = get_density('e')

    print(f"\n{'='*80}")
    print("FINAL DENSITIES vs TARGETS")
    print(f"{'='*80}")
    print(f"{'Species':<10} {'Model':<15} {'Target':<15} {'Ratio':<10} {'Error':<10}")
    print(f"{'-'*80}")

    ratio_H = nH / TARGETS['H']
    ratio_CH = nCH / TARGETS['CH']
    ratio_C2 = nC2 / TARGETS['C2']

    err_H = abs(np.log10(ratio_H)) if nH > 0 else 999
    err_CH = abs(np.log10(ratio_CH)) if nCH > 0 else 999
    err_C2 = abs(np.log10(ratio_C2)) if nC2 > 0 else 999

    print(f"{'H':<10} {nH:<15.2e} {TARGETS['H']:<15.2e} {ratio_H:<10.2f} {err_H:<10.3f}")
    print(f"{'CH':<10} {nCH:<15.2e} {TARGETS['CH']:<15.2e} {ratio_CH:<10.2f} {err_CH:<10.3f}")
    print(f"{'C2':<10} {nC2:<15.2e} {TARGETS['C2']:<15.2e} {ratio_C2:<10.2f} {err_C2:<10.3f}")
    print(f"{'-'*80}")

    total_error = err_H + err_CH + err_C2
    print(f"{'Total log error:':<40} {total_error:.3f}")

    print(f"\n{'='*80}")
    print("OTHER KEY SPECIES")
    print(f"{'='*80}")
    print(f"  H₂: {nH2:.2e} cm⁻³")
    print(f"  ne: {ne_final:.2e} cm⁻³ (input: {params['ne']:.2e})")

    print(f"\n{'='*80}")
    print("KEY RATIOS")
    print(f"{'='*80}")
    print(f"  H/CH = {nH/nCH:.2e} (target: {TARGETS['H']/TARGETS['CH']:.2e})")
    print(f"  C₂/CH = {nC2/nCH:.2e} (target: {TARGETS['C2']/TARGETS['CH']:.2e})")
    print(f"  H/C₂ = {nH/nC2:.2e} (target: {TARGETS['H']/TARGETS['C2']:.2e})")

else:
    print(f"✗ FAILED! Runtime: {runtime:.2f}s")
    print(f"  Status: {sol.message}")

print(f"\n{'='*80}")
print(f"Estimated sweep times (if all runs ~ {runtime:.1f}s):")
print(f"  Reduced (216 runs):       {runtime*216/60:.1f} minutes")
print(f"  Comprehensive (768 runs): {runtime*768/60:.1f} minutes")
print(f"{'='*80}")
