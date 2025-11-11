#!/usr/bin/env python3
"""
Test that the CRITICAL FIX for C2 production is working:
1. Verify C2-producing reactions are force-included in optimization
2. Test if these reactions can actually produce C2 with current conditions
"""

import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '.')

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


def pressure_to_density(pressure_mTorr, T_K=400):
    """Convert pressure to number density."""
    kB = 1.38064852e-23  # J/K
    Torr_to_Pa = 133.322
    P_Pa = pressure_mTorr * 1e-3 * Torr_to_Pa
    n_m3 = P_Pa / (kB * T_K)
    n_cm3 = n_m3 * 1e-6
    return n_cm3


print("=" * 80)
print("TEST 1: Verify force-included C2 reactions are in tunable set")
print("=" * 80)
print()

# Import the select_tunable_rates function from the fixed optimizer
import importlib.util
spec = importlib.util.spec_from_file_location(
    "optimizer",
    "./optimize_charge_balanced_500mTorr_NO_TIMEOUT.py"
)
optimizer_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(optimizer_module)

selected_rates = optimizer_module.select_tunable_rates()

c2_critical = [
    'CH_CH_C2_H2_cm3_5_4',
    'C_CH_C2_H_cm3_7_4',
    'CH_C_C2_H_cm3_7_9',
    'e_C2H2_C2_H2_cm3_1_16',
    'C2HPlus_e_C2_H_cm3_6_18',
    'C2H_H_C2_H2_cm3_7_47',
]

print(f"Total selected rates: {len(selected_rates)}")
print()

all_included = True
for rxn in c2_critical:
    if rxn in selected_rates:
        print(f"✓ {rxn} is in tunable set")
    else:
        print(f"✗ {rxn} is MISSING from tunable set!")
        all_included = False

print()
if all_included:
    print("SUCCESS: All C2-producing reactions are force-included!")
else:
    print("FAILURE: Some C2-producing reactions are missing!")
    sys.exit(1)

print()
print("=" * 80)
print("TEST 2: Check C2 production with best parameters at steady state")
print("=" * 80)
print()

# Load best result
import json
with open('optimization_results_charge_balanced/best_f49.5.json', 'r') as f:
    best = json.load(f)

# Setup parameters (from optimizer params_base)
params = {
    'P': 500.0,
    'Te': best['Te'],
    'ne': best['Ne'],
    'E_field': best['E_field'],
    'T_gas': 400.0,
    'L_discharge': 0.45,
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
    },
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar',
                'C2H4', 'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C',
                'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H',
                'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus',
                'H2Plus', 'C2H2Star'],
}

n_total = pressure_to_density(500.0)
k = define_rates(params)
params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Get species list
species = params['species']
ns = len(species)

# Find indices
idx_C2 = species.index('C2')
idx_CH = species.index('CH')
idx_C = species.index('C')
idx_C2H2 = species.index('C2H2')
idx_C2HPlus = species.index('C2HPlus')
idx_C2H = species.index('C2H')
idx_H = species.index('H')
idx_e = species.index('e')

print(f"Using best result parameters:")
print(f"  Te = {best['Te']:.3f} eV")
print(f"  Ne = {best['Ne']:.2e} cm⁻³")
print(f"  E = {best['E_field']:.1f} V/cm")
print()

# Set initial conditions from best result
y0 = np.array([best['all_densities'].get(s, 1e3) for s in species])

print("Initial conditions (precursors for C2 production):")
print(f"  C2    = {y0[idx_C2]:.2e} cm⁻³")
print(f"  CH    = {y0[idx_CH]:.2e} cm⁻³")
print(f"  C     = {y0[idx_C]:.2e} cm⁻³")
print(f"  C2H2  = {y0[idx_C2H2]:.2e} cm⁻³")
print(f"  C2H+  = {y0[idx_C2HPlus]:.2e} cm⁻³")
print(f"  C2H   = {y0[idx_C2H]:.2e} cm⁻³")
print(f"  H     = {y0[idx_H]:.2e} cm⁻³")
print(f"  e     = {y0[idx_e]:.2e} cm⁻³")
print()

# Calculate reaction rates for C2 production
print("C2 PRODUCTION RATES (with default rate coefficients):")
print()

# CH + CH → C2 + H2
k_CH_CH = k.get('CH_CH_C2_H2_cm3_5_4', 0)
rate_CH_CH = k_CH_CH * y0[idx_CH] * y0[idx_CH]
print(f"1. CH + CH → C2 + H2")
print(f"   k = {k_CH_CH:.2e} cm³/s")
print(f"   rate = {rate_CH_CH:.2e} cm⁻³/s")
print()

# e + C2H2 → C2 + H2
k_e_C2H2 = k.get('e_C2H2_C2_H2_cm3_1_16', 0)
rate_e_C2H2 = k_e_C2H2 * y0[idx_e] * y0[idx_C2H2]
print(f"2. e + C2H2 → C2 + H2")
print(f"   k = {k_e_C2H2:.2e} cm³/s")
print(f"   rate = {rate_e_C2H2:.2e} cm⁻³/s")
print()

# C2H+ + e → C2 + H
k_C2HPlus_e = k.get('C2HPlus_e_C2_H_cm3_6_18', 0)
rate_C2HPlus_e = k_C2HPlus_e * y0[idx_C2HPlus] * y0[idx_e]
print(f"3. C2H+ + e → C2 + H")
print(f"   k = {k_C2HPlus_e:.2e} cm³/s")
print(f"   rate = {rate_C2HPlus_e:.2e} cm⁻³/s")
print()

# C + CH → C2 + H
k_C_CH = k.get('C_CH_C2_H_cm3_7_4', 0)
rate_C_CH = k_C_CH * y0[idx_C] * y0[idx_CH]
print(f"4. C + CH → C2 + H")
print(f"   k = {k_C_CH:.2e} cm³/s")
print(f"   rate = {rate_C_CH:.2e} cm⁻³/s")
print()

total_C2_production = rate_CH_CH + rate_e_C2H2 + rate_C2HPlus_e + rate_C_CH
print(f"TOTAL C2 PRODUCTION: {total_C2_production:.2e} cm⁻³/s")
print()

# C2 wall loss
if 'loss_C2_11_3' in k:
    k_C2_wall = k['loss_C2_11_3']
    rate_C2_wall = k_C2_wall * y0[idx_C2]
    print(f"C2 WALL LOSS: {rate_C2_wall:.2e} cm⁻³/s (k = {k_C2_wall:.2e} s⁻¹)")
    print()

    if total_C2_production > 0:
        balance = total_C2_production / rate_C2_wall
        print(f"Production/Loss ratio: {balance:.2e}")
        if balance < 1.0:
            print("⚠ WARNING: C2 loss exceeds production at these densities!")
        print()

print("=" * 80)
print("DIAGNOSIS")
print("=" * 80)
print()

if total_C2_production < 1e8:
    print("❌ C2 production is too low!")
    print()
    print("REASONS:")
    if y0[idx_CH] < 1e9:
        print(f"  • CH density too low: {y0[idx_CH]:.2e} cm⁻³ (need ~1e9)")
    if y0[idx_C2HPlus] < 1e6:
        print(f"  • C2H+ density too low: {y0[idx_C2HPlus]:.2e} cm⁻³")
    if k_CH_CH < 1e-11:
        print(f"  • CH+CH rate too slow: {k_CH_CH:.2e} cm³/s")
    print()
    print("RECOMMENDATIONS:")
    print("  1. Boost CH density (currently only {:.2e})".format(y0[idx_CH]))
    print("  2. Boost C2H+ production to enable ion recombination pathway")
    print("  3. Increase C2-producing reaction rates in optimizer")
else:
    print(f"✓ C2 production looks viable: {total_C2_production:.2e} cm⁻³/s")
    print()
    print("The CRITICAL FIX should help when optimizer runs!")
