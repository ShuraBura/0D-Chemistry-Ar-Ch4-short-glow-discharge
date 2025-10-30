#!/usr/bin/env python3
"""
CSB Validation with H DIFFUSION INFLUX from CG

H balance in CSB:
  dH/dt = (local chemistry production) - (local losses) + (diffusion influx from CG)

The diffusion influx is a SOURCE TERM representing H flowing from CG → CSB.

User-specified CSB parameters:
- E-field: 20-300 V/cm
- Te: 0.5-2 eV
- ne: 9e8-9e9 cm⁻³
- H_diffusion_influx: TUNABLE parameter (cm⁻³/s)
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from odefun import PlasmaODE

# CSB Targets
TARGETS = {
    'H': 6.35e14,    # cm⁻³ (from diffusion + chemistry)
    'CH': 9.27e8,    # cm⁻³
    'C2': 5.56e11,   # cm⁻³
}

print("="*80)
print("CSB VALIDATION with H DIFFUSION INFLUX")
print("="*80)
print(f"\nCSB Targets:")
print(f"  H:  {TARGETS['H']:.2e} cm⁻³  (chemistry + diffusion from CG)")
print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")

# Estimate H diffusion influx from CG
# Simple estimate: flux = D * ∇n ≈ D * (n_CG - n_CSB) / L
# For H: D ~ 300 cm²/s, n_CG ~ 9.58e15, n_CSB ~ 6.35e14, L ~ 0.1 cm
# Flux ≈ 300 * (9.58e15 - 6.35e14) / 0.1 ≈ 2.8e19 cm⁻³/s
# But this is very rough - treat as tunable parameter!

H_diffusion_influx_estimates = {
    'low': 1e18,      # cm⁻³/s (weak diffusion)
    'medium': 1e19,   # cm⁻³/s (moderate)
    'high': 1e20,     # cm⁻³/s (strong)
}

# Start with medium estimate
H_diffusion_influx = H_diffusion_influx_estimates['medium']

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
    'ne': 5e9,          # cm⁻³ (user range: 9e8-9e9)
    'Te': 1.0,          # eV (user range: 0.5-2 eV)
    'E_field': 100,     # V/cm (user range: 20-300 V/cm)
    'Tgas': 400,        # K
    'L_diff': 0.1,      # cm
    'gamma_H': 0.01,
    'scale_e_impact': 1.0,
}

print(f"\nCSB Parameters:")
print(f"  ne = {params['ne']:.1e} cm⁻³  (range: 9e8-9e9)")
print(f"  Te = {params['Te']} eV  (range: 0.5-2 eV)")
print(f"  E = {params['E_field']} V/cm  (range: 20-300 V/cm)")
print(f"  Tg = {params['Tgas']} K")
print(f"  H_diffusion_influx = {H_diffusion_influx:.2e} cm⁻³/s")

print("\nDefining rates...")
t0 = time.time()
params['k'] = define_rates_tunable(params)
print(f"  Done in {time.time()-t0:.2f}s")

print("Building reactions...")
t0 = time.time()
params['R'], params['tags'] = build_reactions(params)
print(f"  Done in {time.time()-t0:.2f}s")

species = params['species']

# Modify PlasmaODE to use custom H_drift_gain
class PlasmaODE_CSB(PlasmaODE):
    """PlasmaODE with custom H diffusion influx for CSB."""
    def __init__(self, params, H_diffusion_influx):
        super().__init__(params)
        self.H_drift_gain = H_diffusion_influx  # Override default
        print(f"  H diffusion influx set to: {self.H_drift_gain:.2e} cm⁻³/s")

# Initial densities - start LOWER for CSB
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
set_density('H', 1e14)      # Start LOWER than target (will build up)
set_density('ArPlus', params['ne'] * 0.5)
set_density('H2', 1e13)
set_density('CH', 5e8)
set_density('C2', 5e11)
set_density('CH3', 1e8)

print(f"\nInitial conditions (low H, will build up from diffusion).")

# Shorter integration time for quick test
t_span = (0, 1)  # Just 1 second!
t_eval = [1]

print(f"\nIntegrating CSB with H diffusion...")
print(f"  Time span: 0 to {t_span[1]:.1f} s (short test)")
print(f"  Method: BDF (stiff)")

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
    max_step=0.1  # Smaller max_step for stability
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
    nH2 = get_density('H2')

    print(f"\n{'='*80}")
    print("RESULTS (after 1 second)")
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

    if total_error < 1.0:
        print("\n✓ GOOD! Within factor 10")
    elif total_error < 2.0:
        print("\n~ Reasonable, needs tuning")
    else:
        print("\n⚠️  Needs significant tuning")

    print(f"\nOther species:")
    print(f"  H₂: {nH2:.2e} cm⁻³")

else:
    print(f"✗ FAILED! Runtime: {runtime:.2f}s")
    print(f"  {sol.message}")

print(f"\n{'='*80}")
print(f"CSB with H diffusion runtime: {runtime:.2f}s")
print(f"\nNext: Sweep H_diffusion_influx parameter:")
print(f"  Low:    {H_diffusion_influx_estimates['low']:.1e} cm⁻³/s")
print(f"  Medium: {H_diffusion_influx_estimates['medium']:.1e} cm⁻³/s")
print(f"  High:   {H_diffusion_influx_estimates['high']:.1e} cm⁻³/s")
print(f"\nPlus sweep: ne, Te, E-field")
print(f"{'='*80}")
