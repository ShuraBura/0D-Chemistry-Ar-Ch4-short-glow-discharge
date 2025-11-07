#!/usr/bin/env python3
"""
Analyze CH chemistry to understand why it's 59x too high.

Strategy:
1. Run baseline simulation
2. Identify top CH production and loss pathways
3. Calculate production/loss rates
4. Identify which rates need tuning to reduce CH
"""

import numpy as np
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


def analyze_ch_chemistry():
    """Run simulation and analyze CH pathways."""

    print("=" * 80)
    print(" ANALYZING CH CHEMISTRY")
    print("=" * 80)
    print("\nTarget: CH = 1.0e9 cm⁻³")
    print("Current: CH = 6.06e10 cm⁻³ (59× too high!)")
    print("\nLet's find out why...\n")

    # Setup parameters - using baseline from optimization
    params = {
        'E_field': 50.0,
        'L_discharge': 0.45,
        'ne': 3.3e9,
        'Te': 1.0,  # Default electron temperature
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                    'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                    'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                    'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
        }
    }

    # Get rate constants
    k = define_rates(params)
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

    set_density('e', 3.3e9)
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

    # Run simulation
    print("Running simulation...")
    ode_func = PlasmaODE_Optimized(params)
    sol = solve_ivp(
        ode_func,
        (0, 100),
        y0,
        method='BDF',
        rtol=1e-5,
        atol=1e-6,
        max_step=10.0
    )

    if not sol.success:
        print("Simulation failed!")
        return

    y_final = sol.y[:, -1]

    def get_density(name):
        try:
            idx = species.index(name)
            return y_final[idx]
        except ValueError:
            return 0.0

    # Calculate reaction rates
    rate_constants = np.array([params['k'][tag] for tag in params['tags']])
    reactions = params['R']

    reaction_rates = rate_constants.copy()
    for rxn_idx, reaction in enumerate(reactions):
        react_species = np.where(reaction.reactants > 0)[0]
        for sp_idx in react_species:
            coeff = reaction.reactants[sp_idx]
            reaction_rates[rxn_idx] *= y_final[sp_idx] ** coeff

    # Analyze CH
    ch_idx = species.index('CH')
    ch_density = y_final[ch_idx]

    print(f"✓ Simulation complete\n")
    print(f"Final CH density: {ch_density:.2e} cm⁻³")
    print(f"Target CH density: 1.0e9 cm⁻³")
    print(f"Ratio: {ch_density/1.0e9:.1f}× too high\n")

    # Find CH production and loss reactions
    production_rxns = []
    loss_rxns = []

    for rxn_idx, reaction in enumerate(reactions):
        net_change = reaction.products[ch_idx] - reaction.reactants[ch_idx]

        if net_change > 0:
            contribution = net_change * reaction_rates[rxn_idx]
            production_rxns.append((params['tags'][rxn_idx], contribution, k[params['tags'][rxn_idx]]))
        elif net_change < 0:
            contribution = abs(net_change) * reaction_rates[rxn_idx]
            loss_rxns.append((params['tags'][rxn_idx], contribution, k[params['tags'][rxn_idx]]))

    production_rxns.sort(key=lambda x: x[1], reverse=True)
    loss_rxns.sort(key=lambda x: x[1], reverse=True)

    total_production = sum(x[1] for x in production_rxns)
    total_loss = sum(x[1] for x in loss_rxns)

    print("=" * 80)
    print(" CH PRODUCTION PATHWAYS")
    print("=" * 80)
    print(f"\nTotal CH production: {total_production:.2e} cm⁻³/s\n")
    print(f"{'Rank':<6} {'Reaction':<45} {'Rate (cm⁻³/s)':<18} {'%':<8} {'k (cm³/s)':<12}")
    print("-" * 95)

    for i, (tag, rate, k_val) in enumerate(production_rxns[:15], 1):
        pct = (rate / total_production) * 100
        print(f"{i:<6} {tag:<45} {rate:<18.2e} {pct:<8.1f} {k_val:<12.2e}")

    print("\n" + "=" * 80)
    print(" CH LOSS PATHWAYS")
    print("=" * 80)
    print(f"\nTotal CH loss: {total_loss:.2e} cm⁻³/s\n")
    print(f"{'Rank':<6} {'Reaction':<45} {'Rate (cm⁻³/s)':<18} {'%':<8} {'k (cm³/s)':<12}")
    print("-" * 95)

    for i, (tag, rate, k_val) in enumerate(loss_rxns[:15], 1):
        pct = (rate / total_loss) * 100
        print(f"{i:<6} {tag:<45} {rate:<18.2e} {pct:<8.1f} {k_val:<12.2e}")

    # Imbalance
    net_rate = total_production - total_loss
    print("\n" + "=" * 80)
    print(" NET BALANCE")
    print("=" * 80)
    print(f"\nProduction: {total_production:.2e} cm⁻³/s")
    print(f"Loss:       {total_loss:.2e} cm⁻³/s")
    print(f"Net:        {net_rate:.2e} cm⁻³/s ({net_rate/total_production*100:+.1f}%)")

    if abs(net_rate) < 0.01 * total_production:
        print("\n✓ CH is at steady state (production ≈ loss)")
    else:
        print(f"\n⚠ CH is {'accumulating' if net_rate > 0 else 'depleting'}")

    # Recommendations
    print("\n" + "=" * 80)
    print(" RECOMMENDATIONS TO REDUCE CH")
    print("=" * 80)

    print("\nTop 3 production pathways contribute {:.1f}% of CH production:".format(
        sum(x[1] for x in production_rxns[:3]) / total_production * 100
    ))

    for i, (tag, rate, k_val) in enumerate(production_rxns[:3], 1):
        pct = (rate / total_production) * 100
        print(f"  {i}. {tag:<45} ({pct:.1f}%)")

    print("\nStrategies:")
    print("  1. REDUCE these production rates (within literature bounds)")
    print("  2. INCREASE CH loss rates (especially top pathways)")
    print("  3. Adjust precursor densities (H, CH2, C2, etc.)")

    # Check key densities
    print("\n" + "=" * 80)
    print(" KEY SPECIES DENSITIES")
    print("=" * 80)
    print(f"\n{'Species':<10} {'Density (cm⁻³)':<20} {'Notes':<40}")
    print("-" * 70)

    key_species = ['H', 'CH', 'CH2', 'CH3', 'C', 'C2', 'C2H2']
    for sp in key_species:
        dens = get_density(sp)
        note = ""
        if sp == 'CH':
            note = "TARGET: reduce to 1.0e9"
        elif sp == 'H':
            note = "Major reactant in CH production"
        elif sp == 'CH2':
            note = "Major CH precursor"
        print(f"{sp:<10} {dens:<20.2e} {note:<40}")

    # Specific rate analysis
    print("\n" + "=" * 80)
    print(" SPECIFIC RATE RECOMMENDATIONS")
    print("=" * 80)

    # Get literature database
    from rate_database_complete import get_complete_rate_database
    db = get_complete_rate_database()

    print("\nFor the top 3 CH production reactions, here are the literature ranges:\n")

    for i, (tag, rate, k_val) in enumerate(production_rxns[:3], 1):
        if tag in db:
            rate_info = db[tag]
            current_val = k_val
            min_val = rate_info.min
            max_val = rate_info.max
            range_factor = max_val / min_val if min_val > 0 else 1.0

            print(f"{i}. {tag}")
            print(f"   Current k: {current_val:.2e} cm³/s")
            print(f"   Literature: [{min_val:.2e}, {max_val:.2e}] cm³/s")
            print(f"   Range: {range_factor:.1f}×")

            if current_val > min_val:
                reduction = min_val / current_val
                new_ch_production = rate * reduction
                print(f"   → If reduced to MIN: CH production would drop {reduction:.2f}× (to {new_ch_production:.2e} cm⁻³/s)")
            else:
                print(f"   → Already at minimum")
            print()

    print("\n" + "=" * 80)
    print(" CONCLUSION")
    print("=" * 80)
    print("\nTo reduce CH by 59×, we need to:")
    print("  1. Reduce top CH production rates (tune within literature bounds)")
    print("  2. Possibly reduce H density (affects CH2+H→CH and C2+H→CH)")
    print("  3. Increase CH loss pathways")
    print("  4. Or use higher Te to change rate balance (now that we have Te-dependent rates!)")
    print("\nNext step: Run multi-parameter optimization with these rates as tunables.\n")


if __name__ == '__main__':
    analyze_ch_chemistry()
