#!/usr/bin/env python3
"""
Analyze C2 production pathways and rank by contribution
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from scipy.integrate import solve_ivp

# Corrected parameters
PRESSURE_MTORR = 400
TGAS_K = 570
NE = 2.3e9
TE_EV = 1.5

def pressure_to_density(pressure_mTorr, T_K):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

def run_simulation_with_analysis():
    """Run simulation and analyze C2 production"""

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
    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    species = params['species']
    y0 = np.ones(len(species)) * 1e3
    y0[species.index('Ar')] = n_total * 0.97
    y0[species.index('CH4')] = n_total * 0.03
    y0[species.index('e')] = NE

    # Run simulation
    ode_func = PlasmaODE_Optimized(params)
    ode_func.H_drift_gain = 5.7e16

    sol = solve_ivp(ode_func, (0, 500), y0, method='BDF',
                   rtol=1e-7, atol=1e-9, max_step=1.0)

    if not sol.success:
        print("Simulation failed!")
        return

    # Get final densities
    y_final = sol.y[:, -1]

    # Calculate reaction rates at steady state
    C2_idx = species.index('C2')

    print("="*80)
    print("C2 PRODUCTION PATHWAY ANALYSIS")
    print("="*80)
    print(f"\nSimulation conditions:")
    print(f"  ne = {NE:.2e} cm^-3")
    print(f"  Te = {TE_EV} eV")
    print(f"  P = {PRESSURE_MTORR} mTorr")
    print(f"  Tgas = {TGAS_K} K")
    print(f"\nFinal C2 density: {y_final[C2_idx]:.2e} cm^-3")
    print()

    # Analyze each reaction
    production_reactions = []
    destruction_reactions = []

    for rxn_idx, reaction in enumerate(params['R']):
        # Check if this reaction produces or destroys C2
        c2_change = reaction.products[C2_idx] - reaction.reactants[C2_idx]

        if c2_change == 0:
            continue  # C2 not involved

        # Calculate reaction rate
        rate_constant = params['k'][params['tags'][rxn_idx]]

        # Calculate reactant term: k * [A] * [B] * ...
        rate = rate_constant
        reactant_species = []
        for sp_idx, stoich in enumerate(reaction.reactants):
            if stoich > 0:
                rate *= y_final[sp_idx] ** stoich
                reactant_species.append(f"{species[sp_idx]}^{int(stoich)}" if stoich > 1 else species[sp_idx])

        # Rate of C2 change (cm^-3 s^-1)
        c2_rate = rate * c2_change

        # Build reaction string
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

    # Calculate totals
    total_production = sum(r['rate'] for r in production_reactions)
    total_destruction = sum(r['rate'] for r in destruction_reactions)

    # Sort by rate
    production_reactions.sort(key=lambda x: x['rate'], reverse=True)
    destruction_reactions.sort(key=lambda x: x['rate'], reverse=True)

    # Print production reactions
    print("="*80)
    print(f"C2 PRODUCTION REACTIONS (Total: {total_production:.2e} cm^-3 s^-1)")
    print("="*80)

    for i, rxn in enumerate(production_reactions[:10]):
        pct = 100 * rxn['rate'] / total_production
        print(f"\nRank {i+1}: {pct:6.2f}%")
        print(f"  {rxn['reaction']}")
        print(f"  Rate: {rxn['rate']:.2e} cm^-3 s^-1")
        print(f"  k = {rxn['k']:.2e} cm^3/s")
        print(f"  Tag: {rxn['tag']}")

    # Print destruction reactions
    print("\n" + "="*80)
    print(f"C2 DESTRUCTION REACTIONS (Total: {total_destruction:.2e} cm^-3 s^-1)")
    print("="*80)

    for i, rxn in enumerate(destruction_reactions[:10]):
        pct = 100 * rxn['rate'] / total_destruction
        print(f"\nRank {i+1}: {pct:6.2f}%")
        print(f"  {rxn['reaction']}")
        print(f"  Rate: {rxn['rate']:.2e} cm^-3 s^-1")
        print(f"  k = {rxn['k']:.2e} cm^3/s")
        print(f"  Tag: {rxn['tag']}")

    # Net balance
    print("\n" + "="*80)
    print("NET C2 BALANCE")
    print("="*80)
    print(f"  Total production:   {total_production:.2e} cm^-3 s^-1")
    print(f"  Total destruction:  {total_destruction:.2e} cm^-3 s^-1")
    print(f"  Net rate:           {total_production - total_destruction:.2e} cm^-3 s^-1")
    print(f"  Production/Destruction ratio: {total_production/total_destruction:.3f}")

if __name__ == '__main__':
    run_simulation_with_analysis()
