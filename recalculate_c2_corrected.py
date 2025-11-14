#!/usr/bin/env python3
"""
Recalculate C2 pathways with CORRECTED chemistry
After fixing all the Baulch-identified errors
"""

import numpy as np
from scipy.integrate import solve_ivp
from build_reactions import build_reactions
from define_rates import define_rates
from odefun_optimized import PlasmaODE_Optimized

# Set conditions
PRESSURE_MTORR = 400
TGAS_K = 570
NE = 2.3e9
TE_EV = 1.3

def pressure_to_density(pressure_mTorr, T_K):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

# Build params dict
n_total = pressure_to_density(PRESSURE_MTORR, TGAS_K)

params = {
    'P': PRESSURE_MTORR,
    'Tgas': TGAS_K,
    'Te': TE_EV,
    'ne': NE,
    'E_field': 150,
    'L_discharge': 0.45,
    'mobilities': {
        'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
        'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
        'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
        'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
    },
    'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4',
                'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3',
                'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3',
                'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'C3H5', 'HMinus', 'C2H5Plus',
                'CH2Plus', 'C4H', 'H2Plus', 'C2H6Plus'],
    'stick_coeffs': {
        'ArPlus': 1.0, 'CH4Plus': 1.0, 'CH3Plus': 1.0, 'CH5Plus': 1.0,
        'ArHPlus': 1.0, 'CH3Minus': 1.0, 'C2': 0.01, 'CH': 0.001, 'H': 0.0,
        'C': 0.01, 'CH2': 0.001, 'CH3': 0.001, 'C2H': 0.001, 'C2H2': 0.001,
        'C2H3': 0.001, 'C2H4': 0.001, 'C2H5': 0.001, 'C2H6': 0.001, 'C3H2': 0.001,
        'CHPlus': 1.0, 'C3H': 0.001, 'C4H2': 0.001, 'C3H3': 0.001, 'C3H4': 0.001,
        'C3': 0.001, 'C2H4Plus': 1.0, 'C2H3Plus': 1.0, 'C2HPlus': 1.0,
        'C3H5': 0.001, 'HMinus': 1.0, 'C2H5Plus': 1.0, 'CH2Plus': 1.0,
        'C4H': 0.001, 'H3Plus': 1.0, 'H2Plus': 1.0, 'C2H6Plus': 1.0
    },
    'n_Ar': n_total * 0.97,
    'n_CH4': n_total * 0.03,
    'H_drift_gain': 5.7e16
}

species_names = params['species']

# Build reactions with corrected rates
reactions, rate_names = build_reactions(params)
k = define_rates(params)

print("="*80)
print("C2 PATHWAY RECALCULATION WITH CORRECTED CHEMISTRY")
print("="*80)
print()
print("CORRECTIONS APPLIED:")
print("1. ❌ REMOVED: C2 + H → CH + C (was 93% of C2 destruction)")
print("2. ❌ REMOVED: CH + CH2 → C2 + H2 + H (was 21% of C2 production)")
print("3. ✓ FIXED: CH2 + CH2 → C2H2 + H2 (was incorrectly C2 + H2 + H2)")
print("4. ✓ FIXED: CH + CH → C2H2 (was incorrectly C2 + H2)")
print("5. ✓ FIXED: C2H2 + C → C2 rate: 1.0e-10 → 2.0e-10")
print("6. ✓ FIXED: CH + H → C + H2 rate: 1.2e-10 → 2.0e-10")
print("7. ✓ FIXED: H + C2H2 → C2 Arrhenius form (was 7.5e9× too high)")
print()

# Check key rates
print("KEY RATE CONSTANTS AT T=570K:")
print(f"  C2 + H → CH + C:        k = {k['C2_H_CH_C_cm3_7_6']:.2e} cm³/s (SHOULD BE 0!)")
print(f"  CH + CH2 → C2:          k = {k['CH2_CH_C2_H2_H_cm3_7_26']:.2e} cm³/s (SHOULD BE 0!)")
print(f"  C2H2 + C → C2 + CH2:    k = {k['C2H2_C_C2_CH2_cm3_7_19']:.2e} cm³/s (corrected)")
print(f"  CH + H → C + H2:        k = {k['CH_H_C_H2_cm3_7_3']:.2e} cm³/s (corrected)")
print(f"  H + C2H2 → C2 + H2 + H: k = {k['C2H2_H_C2_H2_H_cm3_7_50']:.2e} cm³/s (corrected Arrhenius)")
print()

# Solve steady state
y0 = np.zeros(len(species_names))
y0[species_names.index('Ar')] = params['n_Ar']
y0[species_names.index('CH4')] = params['n_CH4']

print("Solving for steady state...")
sol = solve_ivp(
    lambda t, y: odefun(t, y, params, reactions),
    [0, 1.0],  # seconds
    y0,
    method='BDF',
    rtol=1e-8,
    atol=1e-10,
    dense_output=True
)

if not sol.success:
    print(f"❌ Solver failed: {sol.message}")
    exit(1)

# Get final densities
y_final = sol.y[:, -1]

# Extract key species
idx_H = species_names.index('H')
idx_CH = species_names.index('CH')
idx_C2 = species_names.index('C2')
idx_C2H2 = species_names.index('C2H2')
idx_C = species_names.index('C')
idx_CH2 = species_names.index('CH2')

n_H = y_final[idx_H]
n_CH = y_final[idx_CH]
n_C2 = y_final[idx_C2]
n_C2H2 = y_final[idx_C2H2]
n_C = y_final[idx_C]
n_CH2 = y_final[idx_CH2]

print()
print("="*80)
print("RESULTS")
print("="*80)
print()
print("SPECIES DENSITIES:")
print(f"  H     = {n_H:.2e} cm⁻³  (target: 2.30e14, ratio: {n_H/2.3e14:.2f})")
print(f"  CH    = {n_CH:.2e} cm⁻³  (target: 1.34e9, ratio: {n_CH/1.34e9:.2f})")
print(f"  C2    = {n_C2:.2e} cm⁻³  (target: 5.60e11, ratio: {n_C2/5.6e11:.2f})")
print(f"  C2H2  = {n_C2H2:.2e} cm⁻³")
print(f"  C     = {n_C:.2e} cm⁻³")
print(f"  CH2   = {n_CH2:.2e} cm⁻³")
print()

# Calculate C2 production/destruction rates
print("="*80)
print("C2 PRODUCTION/DESTRUCTION ANALYSIS")
print("="*80)
print()

# Find all reactions that produce/consume C2
c2_production = []
c2_destruction = []

for rxn, name in zip(reactions, rate_names):
    # Check if C2 is produced
    c2_prod = rxn['products'].get(idx_C2, 0) - rxn['reactants'].get(idx_C2, 0)

    if c2_prod > 0:
        # Calculate rate
        rate = rxn['k']
        for sp_idx, nu in rxn['reactants'].items():
            rate *= y_final[sp_idx]**nu
        c2_production.append((name, rate * c2_prod, rate))
    elif c2_prod < 0:
        # Calculate rate
        rate = rxn['k']
        for sp_idx, nu in rxn['reactants'].items():
            rate *= y_final[sp_idx]**nu
        c2_destruction.append((name, rate * abs(c2_prod), rate))

# Sort by magnitude
c2_production.sort(key=lambda x: x[1], reverse=True)
c2_destruction.sort(key=lambda x: x[1], reverse=True)

# Calculate totals
total_prod = sum(r[1] for r in c2_production)
total_dest = sum(r[1] for r in c2_destruction)

print(f"Total C2 production:  {total_prod:.2e} cm⁻³/s")
print(f"Total C2 destruction: {total_dest:.2e} cm⁻³/s")
print(f"Net C2 rate:          {total_prod - total_dest:.2e} cm⁻³/s")
print()

print("TOP C2 PRODUCTION PATHWAYS:")
for i, (name, rate, k_eff) in enumerate(c2_production[:10], 1):
    pct = 100 * rate / total_prod if total_prod > 0 else 0
    print(f"{i:2d}. {name:45s} {rate:10.2e} cm⁻³/s ({pct:5.2f}%)")

print()
print("TOP C2 DESTRUCTION PATHWAYS:")
for i, (name, rate, k_eff) in enumerate(c2_destruction[:10], 1):
    pct = 100 * rate / total_dest if total_dest > 0 else 0
    print(f"{i:2d}. {name:45s} {rate:10.2e} cm⁻³/s ({pct:5.2f}%)")

print()
print("="*80)
print("COMPARISON WITH TARGETS")
print("="*80)
print()

# Calculate RMS error
targets = {
    'H': 2.3e14,
    'CH': 1.34e9,
    'C2': 5.6e11
}

results = {
    'H': n_H,
    'CH': n_CH,
    'C2': n_C2
}

errors = []
for sp in ['H', 'CH', 'C2']:
    error_pct = 100 * abs(results[sp] - targets[sp]) / targets[sp]
    errors.append(error_pct)
    status = "✓" if error_pct < 50 else "✗"
    print(f"{status} {sp:4s}: {results[sp]:.2e} vs {targets[sp]:.2e} (error: {error_pct:6.1f}%)")

rms_error = np.sqrt(np.mean(np.array(errors)**2))
print()
print(f"RMS Error: {rms_error:.1f}%")
print()

# Compare with old results
print("="*80)
print("COMPARISON WITH OLD (INCORRECT) CHEMISTRY")
print("="*80)
print()
print("Old C2 = 2.75×10⁸ cm⁻³ (with wrong rates)")
print(f"New C2 = {n_C2:.2e} cm⁻³ (with corrected rates)")
print(f"Improvement factor: {n_C2/2.75e8:.1f}×")
print()
print("Analysis:")
print("  • Removed fake C2 destruction (C2 + H → CH + C): ~93% decrease")
print("  • Removed fake C2 production (CH + CH2 → C2): ~21% decrease")
print("  • Fixed CH2 + CH2 products: no longer produces C2")
print("  • Fixed CH + CH → C2H2: no longer produces C2")
print("  • Fixed C2H2 + C rate: +100% (2× increase)")
print("  • Fixed H + C2H2 rate: ~-99.9% (Arrhenius form)")
print()
print("Expected net effect:")
print("  • Production: -37% (lost fake pathways) + 18% (C2H2+C fix) ≈ -19%")
print("  • Destruction: -93% (removed C2+H)")
print("  • Net C2 increase: ~10-15× predicted")
print()

if n_C2 > 2.75e8:
    print(f"✓ C2 increased as predicted! Factor: {n_C2/2.75e8:.1f}×")
else:
    print(f"✗ C2 decreased unexpectedly. Factor: {n_C2/2.75e8:.1f}×")

print()
print("="*80)
