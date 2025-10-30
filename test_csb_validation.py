#!/usr/bin/env python3
"""
Quick validation test for CSB (Cathode Sheath Boundary) region
Should be MUCH faster and easier than CG due to:
- Lower H density (15× less)
- More thermal Te (0.5-1 eV)
- Lower E-field
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from odefun import PlasmaODE

# CSB Targets (from EXPERIMENTAL_TARGETS.md)
TARGETS = {
    'H': 6.35e14,    # cm⁻³ (15× lower than CG!)
    'CH': 9.27e8,    # cm⁻³ (2× higher than CG)
    'C2': 5.56e11,   # cm⁻³ (3.7× higher than CG)
}

print("="*80)
print("QUICK VALIDATION TEST - CSB (Cathode Sheath Boundary) Region")
print("="*80)
print(f"\nCSB Targets:")
print(f"  H:  {TARGETS['H']:.2e} cm⁻³  (15× LOWER than CG)")
print(f"  CH: {TARGETS['CH']:.2e} cm⁻³  (2× higher than CG)")
print(f"  C2: {TARGETS['C2']:.2e} cm⁻³  (3.7× higher than CG)")

# CSB parameters - more thermal, lower E-field
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
    'ne': 2e10,         # cm⁻³ (HIGHER in CSB - easier to measure)
    'Te': 0.8,          # eV (MORE THERMAL - easier to model!)
    'E_field': 50,      # V/cm (LOWER - bulk plasma)
    'Tgas': 400,        # K (cooler in bulk)
    'L_diff': 0.1,      # cm (default)
    'gamma_H': 0.01,    # wall recombination
    'scale_e_impact': 1.0
}

print(f"\nCSB Parameters (More Thermal):")
print(f"  ne = {params['ne']:.1e} cm⁻³  (HIGHER - easier to measure)")
print(f"  Te = {params['Te']} eV  (MORE THERMAL)")
print(f"  E = {params['E_field']} V/cm  (LOWER - bulk plasma)")
print(f"  Tg = {params['Tgas']} K")
print(f"  L_diff = {params['L_diff']} cm")
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

# CSB initial guesses - lower H, higher CH/C2
set_density('ArPlus', params['ne'] * 0.5)
set_density('H', 1e14)      # Start lower for CSB
set_density('H2', 1e13)
set_density('CH', 5e8)      # Start higher for CSB
set_density('C2', 5e11)     # Start higher for CSB
set_density('CH3', 1e8)
set_density('CH4Plus', params['ne'] * 0.05)

print(f"\nInitial conditions set (CSB-optimized).")

# Time span
t_span = (0, 10)  # Start with shorter time - 10s
t_eval = [10]

print(f"\nIntegrating ODE system...")
print(f"  Time span: 0 to {t_span[1]:.1f} s (shorter for quick test)")
print(f"  Method: BDF (stiff)")

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
    print("FINAL DENSITIES vs CSB TARGETS")
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

    # Assessment
    print(f"\n{'='*80}")
    print("ASSESSMENT")
    print(f"{'='*80}")
    if total_error < 1.0:
        print("  🎉 EXCELLENT! All species within factor 10")
    elif total_error < 2.0:
        print("  ✓ GOOD! Model in reasonable range")
    else:
        print("  ⚠️  Needs tuning")

    print(f"\n{'='*80}")
    print("OTHER KEY SPECIES")
    print(f"{'='*80}")
    print(f"  H₂: {nH2:.2e} cm⁻³")
    print(f"  ne: {ne_final:.2e} cm⁻³ (input: {params['ne']:.2e})")

    print(f"\n{'='*80}")
    print("KEY RATIOS")
    print(f"{'='*80}")
    print(f"  H/CH (model): {nH/nCH:.2e}  (target: {TARGETS['H']/TARGETS['CH']:.2e})")
    print(f"  C₂/CH (model): {nC2/nCH:.2e}  (target: {TARGETS['C2']/TARGETS['CH']:.2e})")
    print(f"  H/C₂ (model): {nH/nC2:.2e}  (target: {TARGETS['H']/TARGETS['C2']:.2e})")

else:
    print(f"✗ FAILED! Runtime: {runtime:.2f}s")
    print(f"  Status: {sol.message}")

print(f"\n{'='*80}")
print(f"CSB runtime: {runtime:.1f}s (much faster than CG!)")
print(f"\nEstimated sweep times for CSB (if all runs ~ {runtime:.1f}s):")
print(f"  Reduced (216 runs):       {runtime*216/60:.1f} minutes")
print(f"  Comprehensive (768 runs): {runtime*768/60:.1f} minutes")
print(f"{'='*80}")
