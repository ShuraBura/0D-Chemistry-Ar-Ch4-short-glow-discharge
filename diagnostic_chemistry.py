#!/usr/bin/env python3
"""
Detailed chemistry diagnostic showing:
1. All species densities
2. Reaction rates and contributions to target species (H, CH, C2)
3. Production/loss breakdown for each target species
"""

import numpy as np
from scipy.integrate import solve_ivp

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


def run_diagnostic(ne=3.3e9, E_field=50):
    """Run simulation and generate detailed chemistry diagnostics."""

    # Setup parameters
    params = {
        'E_field': E_field,
        'L_discharge': 0.45,
        'ne': ne,
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                    'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                    'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                    'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'C2H2Star'],
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'CH3Minus': 3000, 'HMinus': 3000
        }
    }

    # Load rates with literature corrections
    k = define_rates(params)
    db = get_complete_rate_database()

    # Ensure all rates within literature bounds
    for name, rate_db in db.items():
        if name in k:
            if k[name] < rate_db.min:
                k[name] = rate_db.min
            elif k[name] > rate_db.max:
                k[name] = rate_db.max

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    # Initial conditions
    species = params['species']
    ns = len(species)
    y0 = np.ones(ns) * 1e3

    def set_density(name, value):
        try:
            idx = species.index(name)
            y0[idx] = value
        except ValueError:
            pass

    set_density('e', ne)
    set_density('Ar', 0.85 * 9.66e15)
    set_density('CH4', 0.15 * 9.66e15)
    set_density('ArPlus', 1e7)
    set_density('CH4Plus', 1e5)
    set_density('CH3Plus', 1e5)
    set_density('CH5Plus', 1e3)
    set_density('ArHPlus', 5e5)
    set_density('CH3Minus', 5e4)
    set_density('H2', 1e12)
    set_density('ArStar', 5e6)
    set_density('H', 1e11)
    set_density('C2', 5e7)
    set_density('CH', 5e4)
    set_density('C2H4', 5e7)
    set_density('C2H6', 1e6)
    set_density('CH2', 1e11)
    set_density('C2H2', 1e12)
    set_density('C2H5', 1e6)
    set_density('CH3', 5e7)
    set_density('C', 5e7)

    print("=" * 80)
    print(" DETAILED CHEMISTRY DIAGNOSTIC")
    print("=" * 80)
    print(f"\nParameters:")
    print(f"  Ne: {ne:.2e} cm^-3")
    print(f"  E field: {E_field} V/cm")
    print(f"\nRunning simulation to t=100s...")

    # Run simulation
    ode_func = PlasmaODE_Optimized(params)
    sol = solve_ivp(
        ode_func,
        (0, 100),
        y0,
        method='BDF',
        rtol=1e-6,
        atol=1e-7
    )

    if not sol.success:
        print(f"ERROR: Solver failed - {sol.message}")
        return

    y_final = sol.y[:, -1]

    # ========== PART 1: ALL SPECIES DENSITIES ==========
    print("\n" + "=" * 80)
    print(" PART 1: ALL SPECIES DENSITIES (t=100s)")
    print("=" * 80)

    # Sort species by density
    species_densities = [(species[i], y_final[i]) for i in range(ns)]
    species_densities.sort(key=lambda x: x[1], reverse=True)

    print("\n{:20s} {:>15s} {:>15s}".format("Species", "Density (cm^-3)", "Fraction"))
    print("-" * 80)
    total_density = sum(y_final)
    for sp_name, dens in species_densities:
        fraction = dens / total_density
        print(f"{sp_name:20s} {dens:15.4e} {fraction:15.4e}")

    # ========== PART 2: REACTION RATE ANALYSIS ==========
    print("\n" + "=" * 80)
    print(" PART 2: REACTION RATES (cm^-3 s^-1)")
    print("=" * 80)

    # Compute all reaction rates at final state
    rate_constants = np.array([k[tag] for tag in params['tags']])
    reactions = params['R']

    reaction_rates = rate_constants.copy()
    for rxn_idx, reaction in enumerate(reactions):
        react_species = np.where(reaction.reactants > 0)[0]
        for sp_idx in react_species:
            coeff = reaction.reactants[sp_idx]
            reaction_rates[rxn_idx] *= y_final[sp_idx] ** coeff

    # Find dominant reactions
    sorted_rxn_idx = np.argsort(reaction_rates)[::-1]

    print("\nTop 20 fastest reactions:")
    print("{:5s} {:50s} {:>15s}".format("#", "Reaction Tag", "Rate (cm^-3 s^-1)"))
    print("-" * 80)
    for i, idx in enumerate(sorted_rxn_idx[:20], 1):
        tag = params['tags'][idx]
        rate = reaction_rates[idx]
        print(f"{i:5d} {tag:50s} {rate:15.4e}")

    # ========== PART 3: TARGET SPECIES ANALYSIS ==========
    target_species_names = ['H', 'CH', 'C2']

    for target_sp in target_species_names:
        try:
            sp_idx = species.index(target_sp)
        except ValueError:
            continue

        print("\n" + "=" * 80)
        print(f" PART 3: {target_sp} PRODUCTION AND LOSS ANALYSIS")
        print("=" * 80)

        # Categorize reactions by production/loss
        production_rxns = []
        loss_rxns = []

        for rxn_idx, reaction in enumerate(reactions):
            net_change = reaction.products[sp_idx] - reaction.reactants[sp_idx]

            if net_change > 0:
                # Net production
                contribution = net_change * reaction_rates[rxn_idx]
                production_rxns.append((rxn_idx, contribution, params['tags'][rxn_idx]))
            elif net_change < 0:
                # Net loss
                contribution = abs(net_change) * reaction_rates[rxn_idx]
                loss_rxns.append((rxn_idx, contribution, params['tags'][rxn_idx]))

        # Sort by contribution
        production_rxns.sort(key=lambda x: x[1], reverse=True)
        loss_rxns.sort(key=lambda x: x[1], reverse=True)

        total_production = sum(x[1] for x in production_rxns)
        total_loss = sum(x[1] for x in loss_rxns)

        print(f"\nCurrent {target_sp} density: {y_final[sp_idx]:.4e} cm^-3")
        print(f"Total production rate: {total_production:.4e} cm^-3 s^-1")
        print(f"Total loss rate:       {total_loss:.4e} cm^-3 s^-1")
        print(f"Net rate:              {total_production - total_loss:.4e} cm^-3 s^-1")

        if total_loss > 0:
            lifetime = y_final[sp_idx] / total_loss
            print(f"Lifetime (dens/loss):  {lifetime:.4e} seconds")

        # Top production pathways
        print(f"\nTop {target_sp} production pathways:")
        print("{:5s} {:50s} {:>15s} {:>10s}".format("#", "Reaction", "Rate (cm^-3 s^-1)", "% of Total"))
        print("-" * 90)
        for i, (idx, contrib, tag) in enumerate(production_rxns[:15], 1):
            pct = (contrib / total_production * 100) if total_production > 0 else 0
            print(f"{i:5d} {tag:50s} {contrib:15.4e} {pct:10.2f}%")

        # Top loss pathways
        print(f"\nTop {target_sp} loss pathways:")
        print("{:5s} {:50s} {:>15s} {:>10s}".format("#", "Reaction", "Rate (cm^-3 s^-1)", "% of Total"))
        print("-" * 90)
        for i, (idx, contrib, tag) in enumerate(loss_rxns[:15], 1):
            pct = (contrib / total_loss * 100) if total_loss > 0 else 0
            print(f"{i:5d} {tag:50s} {contrib:15.4e} {pct:10.2f}%")

    # ========== SUMMARY ==========
    print("\n" + "=" * 80)
    print(" SUMMARY")
    print("=" * 80)

    targets = {
        'H': 5.18e13,
        'CH': 1.0e9,
        'C2': 1.3e11
    }

    print("\n{:10s} {:>15s} {:>15s} {:>10s}".format("Species", "Current", "Target", "Ratio"))
    print("-" * 60)
    for sp_name, target_dens in targets.items():
        try:
            idx = species.index(sp_name)
            current = y_final[idx]
            ratio = current / target_dens
            print(f"{sp_name:10s} {current:15.4e} {target_dens:15.4e} {ratio:10.2f}x")
        except ValueError:
            print(f"{sp_name:10s} NOT FOUND")

    print("\n" + "=" * 80)


if __name__ == '__main__':
    import sys

    # Allow command-line arguments for Ne and E
    ne = 3.3e9
    E_field = 50

    if len(sys.argv) > 1:
        ne = float(sys.argv[1])
    if len(sys.argv) > 2:
        E_field = float(sys.argv[2])

    run_diagnostic(ne, E_field)
