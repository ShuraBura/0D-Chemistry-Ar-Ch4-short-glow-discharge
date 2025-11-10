#!/usr/bin/env python3
"""
Analyze C2H2 and C2 chemistry.

User's insight: C2H2 is a big source for C2
If we reduce C2H2, we reduce C2, which then reduces CH (via C2+H→CH)

This could solve multiple problems at once!
"""

import numpy as np
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


def analyze_c2h2_and_c2():
    """Run simulation and analyze C2H2 and C2 pathways."""

    print("=" * 80)
    print(" ANALYZING C2H2 → C2 → CH PATHWAY")
    print("=" * 80)
    print("\nHypothesis: Reducing C2H2 will reduce C2, which will reduce CH")
    print("\nCurrent status:")
    print("  C2H2: Unknown")
    print("  C2:   9.15e11 cm⁻³ (7× too high, target: 1.30e11)")
    print("  CH:   6.06e10 cm⁻³ (60× too high, target: 1.00e9)")
    print("\nLet's trace the chemistry...\n")

    # Setup parameters - using baseline
    params = {
        'E_field': 50.0,
        'L_discharge': 0.45,
        'ne': 3.3e9,
        'Te': 1.0,
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

    print(f"✓ Simulation complete\n")

    # Get densities
    c2h2_density = get_density('C2H2')
    c2_density = get_density('C2')
    ch_density = get_density('CH')
    h_density = get_density('H')

    print("=" * 80)
    print(" CURRENT DENSITIES")
    print("=" * 80)
    print(f"\nC2H2: {c2h2_density:.2e} cm⁻³")
    print(f"C2:   {c2_density:.2e} cm⁻³ (target: 1.30e11, ratio: {c2_density/1.3e11:.2f}x)")
    print(f"CH:   {ch_density:.2e} cm⁻³ (target: 1.00e9, ratio: {ch_density/1.0e9:.2f}x)")
    print(f"H:    {h_density:.2e} cm⁻³")

    # Analyze C2H2
    c2h2_idx = species.index('C2H2')
    c2h2_production = []
    c2h2_loss = []

    for rxn_idx, reaction in enumerate(reactions):
        net_change = reaction.products[c2h2_idx] - reaction.reactants[c2h2_idx]

        if net_change > 0:
            contribution = net_change * reaction_rates[rxn_idx]
            c2h2_production.append((params['tags'][rxn_idx], contribution, k[params['tags'][rxn_idx]]))
        elif net_change < 0:
            contribution = abs(net_change) * reaction_rates[rxn_idx]
            c2h2_loss.append((params['tags'][rxn_idx], contribution, k[params['tags'][rxn_idx]]))

    c2h2_production.sort(key=lambda x: x[1], reverse=True)
    c2h2_loss.sort(key=lambda x: x[1], reverse=True)

    total_c2h2_prod = sum(x[1] for x in c2h2_production)
    total_c2h2_loss = sum(x[1] for x in c2h2_loss)

    print("\n" + "=" * 80)
    print(" C2H2 PRODUCTION PATHWAYS")
    print("=" * 80)
    print(f"\nTotal C2H2 production: {total_c2h2_prod:.2e} cm⁻³/s\n")
    print(f"{'Rank':<6} {'Reaction':<45} {'Rate (cm⁻³/s)':<18} {'%':<8} {'k (cm³/s)':<12}")
    print("-" * 95)

    for i, (tag, rate, k_val) in enumerate(c2h2_production[:10], 1):
        pct = (rate / total_c2h2_prod) * 100
        print(f"{i:<6} {tag:<45} {rate:<18.2e} {pct:<8.1f} {k_val:<12.2e}")

    print("\n" + "=" * 80)
    print(" C2H2 LOSS PATHWAYS")
    print("=" * 80)
    print(f"\nTotal C2H2 loss: {total_c2h2_loss:.2e} cm⁻³/s\n")
    print(f"{'Rank':<6} {'Reaction':<45} {'Rate (cm⁻³/s)':<18} {'%':<8} {'k (cm³/s)':<12}")
    print("-" * 95)

    for i, (tag, rate, k_val) in enumerate(c2h2_loss[:10], 1):
        pct = (rate / total_c2h2_loss) * 100
        print(f"{i:<6} {tag:<45} {rate:<18.2e} {pct:<8.1f} {k_val:<12.2e}")

    # Analyze C2
    c2_idx = species.index('C2')
    c2_production = []
    c2_loss = []

    for rxn_idx, reaction in enumerate(reactions):
        net_change = reaction.products[c2_idx] - reaction.reactants[c2_idx]

        if net_change > 0:
            contribution = net_change * reaction_rates[rxn_idx]
            c2_production.append((params['tags'][rxn_idx], contribution, k[params['tags'][rxn_idx]]))
        elif net_change < 0:
            contribution = abs(net_change) * reaction_rates[rxn_idx]
            c2_loss.append((params['tags'][rxn_idx], contribution, k[params['tags'][rxn_idx]]))

    c2_production.sort(key=lambda x: x[1], reverse=True)
    c2_loss.sort(key=lambda x: x[1], reverse=True)

    total_c2_prod = sum(x[1] for x in c2_production)
    total_c2_loss = sum(x[1] for x in c2_loss)

    print("\n" + "=" * 80)
    print(" C2 PRODUCTION PATHWAYS")
    print("=" * 80)
    print(f"\nTotal C2 production: {total_c2_prod:.2e} cm⁻³/s\n")
    print(f"{'Rank':<6} {'Reaction':<45} {'Rate (cm⁻³/s)':<18} {'%':<8} {'k (cm³/s)':<12}")
    print("-" * 95)

    for i, (tag, rate, k_val) in enumerate(c2_production[:10], 1):
        pct = (rate / total_c2_prod) * 100
        marker = " ← C2H2!" if "C2H2" in tag else ""
        print(f"{i:<6} {tag:<45} {rate:<18.2e} {pct:<8.1f} {k_val:<12.2e}{marker}")

    print("\n" + "=" * 80)
    print(" C2 LOSS PATHWAYS")
    print("=" * 80)
    print(f"\nTotal C2 loss: {total_c2_loss:.2e} cm⁻³/s\n")
    print(f"{'Rank':<6} {'Reaction':<45} {'Rate (cm⁻³/s)':<18} {'%':<8} {'k (cm³/s)':<12}")
    print("-" * 95)

    for i, (tag, rate, k_val) in enumerate(c2_loss[:10], 1):
        pct = (rate / total_c2_loss) * 100
        marker = " ← TO CH!" if "CH" in tag and "CH2" not in tag and "CH3" not in tag else ""
        print(f"{i:<6} {tag:<45} {rate:<18.2e} {pct:<8.1f} {k_val:<12.2e}{marker}")

    # Key pathway analysis
    print("\n" + "=" * 80)
    print(" KEY PATHWAY ANALYSIS")
    print("=" * 80)

    # Find C2H2 → C2 reactions
    c2h2_to_c2_rate = 0.0
    c2h2_to_c2_reactions = []

    for tag, rate, k_val in c2_production:
        if "C2H2" in tag:
            c2h2_to_c2_rate += rate
            c2h2_to_c2_reactions.append((tag, rate, k_val))

    # Find C2 → CH reaction
    c2_to_ch_rate = 0.0
    for tag, rate, k_val in c2_loss:
        if tag == 'C2_H_CH_C_cm3_7_6':
            c2_to_ch_rate = rate

    print(f"\n1. C2H2 → C2 pathway:")
    print(f"   Total rate: {c2h2_to_c2_rate:.2e} cm⁻³/s")
    print(f"   Fraction of C2 production: {c2h2_to_c2_rate/total_c2_prod*100:.1f}%")
    if c2h2_to_c2_reactions:
        print(f"   Main reaction: {c2h2_to_c2_reactions[0][0]}")

    print(f"\n2. C2 → CH pathway:")
    print(f"   C2 + H → CH + C rate: {c2_to_ch_rate:.2e} cm⁻³/s")
    print(f"   Fraction of C2 loss: {c2_to_ch_rate/total_c2_loss*100:.1f}%")

    print(f"\n3. Chain effect:")
    print(f"   If C2H2 reduced by 50%:")
    print(f"     → C2 production reduced by ~{c2h2_to_c2_rate/total_c2_prod*50:.1f}%")
    print(f"     → C2 density might drop to ~{c2_density * (1 - 0.5*c2h2_to_c2_rate/total_c2_prod):.2e}")
    print(f"     → CH production from C2+H would drop proportionally")

    # Literature ranges
    from rate_database_complete import get_complete_rate_database
    db = get_complete_rate_database()

    print("\n" + "=" * 80)
    print(" TUNABLE RATES TO REDUCE C2H2")
    print("=" * 80)
    print("\nTop C2H2 production reactions with tunable ranges:\n")

    for i, (tag, rate, k_val) in enumerate(c2h2_production[:5], 1):
        if tag in db:
            rate_info = db[tag]
            range_factor = rate_info.max / rate_info.min if rate_info.min > 0 else 1.0

            print(f"{i}. {tag}")
            print(f"   Current k: {k_val:.2e} cm³/s")
            print(f"   Literature: [{rate_info.min:.2e}, {rate_info.max:.2e}] cm³/s")
            print(f"   Range: {range_factor:.1f}×")

            if range_factor > 1.5:
                if k_val > rate_info.min * 1.1:
                    reduction = rate_info.min / k_val
                    new_rate = rate * reduction
                    print(f"   → If reduced to MIN: C2H2 production drops {reduction:.2f}× (to {new_rate:.2e} cm⁻³/s)")
                else:
                    print(f"   → Already near minimum")
            else:
                print(f"   → Limited tunability")
            print()

    print("\n" + "=" * 80)
    print(" RECOMMENDATIONS")
    print("=" * 80)
    print("\n1. REDUCE C2H2 production (top 3 reactions contribute {:.1f}%)".format(
        sum(x[1] for x in c2h2_production[:3]) / total_c2h2_prod * 100
    ))
    print("2. INCREASE C2H2 loss pathways")
    print("3. This will CASCADE to reduce:")
    print(f"   → C2 (via reduced C2H2→C2)")
    print(f"   → CH (via reduced C2+H→CH)")
    print("\nThis strategy attacks all three problems at once!")


if __name__ == '__main__':
    analyze_c2h2_and_c2()
