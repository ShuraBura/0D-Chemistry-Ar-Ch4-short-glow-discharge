#!/usr/bin/env python3
"""
Recalculate C2 production pathways with CORRECTED H + C2H2 rate from Baulch

Baulch et al. (2005): k = 1.67e-14 * T^1.64 * exp(-15250/T) cm³/s
At T=570K, this gives k ~ 3.2e-24 cm³/s (not 1e-11!)
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from scipy.integrate import solve_ivp

# Parameters
PRESSURE_MTORR = 400
TGAS_K = 570
NE = 2.3e9
TE_EV = 1.5

def pressure_to_density(pressure_mTorr, T_K):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

def baulch_rate_h_c2h2(T):
    """Correct Baulch rate for H + C2H2 → C2 + H2 + H"""
    return 1.67e-14 * (T**1.64) * np.exp(-15250.0 / T)

print("="*80)
print("RECALCULATING C2 PATHWAYS WITH CORRECT BAULCH RATE")
print("="*80)

print(f"\nH + C2H2 → C2 + H2 + H rate correction:")
print(f"  Model uses:   k = 1.0e-11 cm³/s (constant)")
print(f"  Baulch gives: k = 1.67e-14 * T^1.64 * exp(-15250/T)")

T = TGAS_K
k_correct = baulch_rate_h_c2h2(T)
k_model = 1e-11

print(f"\n  At T = {T} K:")
print(f"    k_Baulch = {k_correct:.2e} cm³/s")
print(f"    k_model  = {k_model:.2e} cm³/s")
print(f"    Ratio: k_model / k_Baulch = {k_model/k_correct:.2e}")
print()
print(f"  >>> Model rate is {k_model/k_correct:.1e}× TOO HIGH! <<<")
print()

# Run simulation with CORRECTED rate
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
                'C3H5', 'C4H', 'C3H6', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
}

k = define_rates(params)

# Apply CORRECTED Baulch rate
k['C2H2_H_C2_H2_H_cm3_7_50'] = k_correct

params['k'] = k
params['R'], params['tags'] = build_reactions(params)

species = params['species']
y0 = np.ones(len(species)) * 1e3
y0[species.index('Ar')] = n_total * 0.97
y0[species.index('CH4')] = n_total * 0.03
y0[species.index('e')] = NE

print("="*80)
print("RUNNING SIMULATION WITH CORRECTED RATE")
print("="*80)

ode_func = PlasmaODE_Optimized(params)
ode_func.H_drift_gain = 5.7e16

sol = solve_ivp(ode_func, (0, 500), y0, method='BDF',
               rtol=1e-7, atol=1e-9, max_step=1.0)

if not sol.success:
    print("Simulation failed!")
    sys.exit(1)

y_final = sol.y[:, -1]
C2_idx = species.index('C2')

print(f"\nFinal C2 density: {y_final[C2_idx]:.2e} cm^-3")
print()

# Analyze C2 production pathways
production_reactions = []
destruction_reactions = []

for rxn_idx, reaction in enumerate(params['R']):
    c2_change = reaction.products[C2_idx] - reaction.reactants[C2_idx]
    if c2_change == 0:
        continue

    rate_constant = params['k'][params['tags'][rxn_idx]]
    rate = rate_constant
    reactant_species = []

    for sp_idx, stoich in enumerate(reaction.reactants):
        if stoich > 0:
            rate *= y_final[sp_idx] ** stoich
            reactant_species.append(f"{species[sp_idx]}^{int(stoich)}" if stoich > 1 else species[sp_idx])

    c2_rate = rate * c2_change

    reactants_str = " + ".join(reactant_species)
    product_species = []
    for sp_idx, stoich in enumerate(reaction.products):
        if stoich > 0:
            product_species.append(f"{species[sp_idx]}^{int(stoich)}" if stoich > 1 else species[sp_idx])
    products_str = " + ".join(product_species)
    reaction_str = f"{reactants_str} → {products_str}"

    if c2_change > 0:
        production_reactions.append({
            'reaction': reaction_str,
            'rate': c2_rate,
            'k': rate_constant,
            'tag': params['tags'][rxn_idx]
        })
    else:
        destruction_reactions.append({
            'reaction': reaction_str,
            'rate': abs(c2_rate),
            'k': rate_constant,
            'tag': params['tags'][rxn_idx]
        })

total_production = sum(r['rate'] for r in production_reactions)
total_destruction = sum(r['rate'] for r in destruction_reactions)

production_reactions.sort(key=lambda x: x['rate'], reverse=True)
destruction_reactions.sort(key=lambda x: x['rate'], reverse=True)

print("="*80)
print(f"C2 PRODUCTION WITH CORRECTED RATE (Total: {total_production:.2e} cm^-3 s^-1)")
print("="*80)

for i, rxn in enumerate(production_reactions[:10]):
    pct = 100 * rxn['rate'] / total_production
    print(f"\nRank {i+1}: {pct:6.2f}%")
    print(f"  {rxn['reaction']}")
    print(f"  Rate: {rxn['rate']:.2e} cm^-3 s^-1")
    print(f"  k = {rxn['k']:.2e} cm^3/s")
    print(f"  Tag: {rxn['tag']}")

print("\n" + "="*80)
print(f"C2 DESTRUCTION (Total: {total_destruction:.2e} cm^-3 s^-1)")
print("="*80)

for i, rxn in enumerate(destruction_reactions[:5]):
    pct = 100 * rxn['rate'] / total_destruction
    print(f"\nRank {i+1}: {pct:6.2f}%")
    print(f"  {rxn['reaction']}")
    print(f"  Rate: {rxn['rate']:.2e} cm^-3 s^-1")
    print(f"  k = {rxn['k']:.2e} cm^3/s")

print("\n" + "="*80)
print("COMPARISON: OLD vs NEW")
print("="*80)
print(f"\nWith WRONG rate (k=1e-11):")
print(f"  H + C2H2 contributed 95% of C2 production")
print(f"\nWith CORRECT rate (k={k_correct:.2e}):")
print(f"  Now see results above ^^^")
print()
print("This reveals the ACTUAL C2 production mechanism!")
print("="*80)
