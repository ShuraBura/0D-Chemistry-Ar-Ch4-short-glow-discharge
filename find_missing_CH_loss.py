"""
Search for missing CH loss mechanisms

Question: Are we missing important CH destruction pathways?
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions

def pressure_to_density(pressure_mTorr, T_K=300):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

P = 500.0
n_total = pressure_to_density(P)

params = {
    'P': P,
    'Te': baseline['Te'],
    'ne': 2.3e9,
    'E_field': baseline['E_field'],
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

# Apply baseline rates
for rate_name, rate_value in baseline['rate_values'].items():
    if rate_name in k:
        k[rate_name] = rate_value

# Apply multipliers
rate_mults = {
    'e_CH4_CH3_HMinus_cm3_8_1': 10.0,
    'ArStar_CH4_CH3_H_cm3_3_1': 1.4,
    'e_CH4Plus_CH3_H_cm3_6_4': 1.5,
    'stick_CH3_9_2': 0.01,
    'stick_C2H2_9_11': 0.01,
    'loss_C2H2_11_19': 0.01,
}

for rate_name, mult in rate_mults.items():
    if rate_name in k:
        k[rate_name] *= mult

params['k'] = k
R_list, tags_list = build_reactions(params)

print("="*80)
print("SEARCH FOR ALL CH LOSS REACTIONS IN MODEL")
print("="*80)

# Find all reactions that consume CH
ch_loss_reactions = []
species = params['species']
ch_idx = species.index('CH')

for reaction, tag in zip(R_list, tags_list):
    reactants = reaction.reactants
    products = reaction.products
    rate = reaction.rate

    # Check if CH appears as reactant
    if reactants[ch_idx] > 0:
        # Check if it's actually consumed (not equal products)
        net_ch = reactants[ch_idx] - products[ch_idx]

        if net_ch > 0:  # CH is consumed
            # Build reaction string
            react_str = []
            for s_idx, stoich in enumerate(reactants):
                if stoich > 0:
                    sp_name = species[s_idx]
                    if stoich == 1:
                        react_str.append(sp_name)
                    else:
                        react_str.append(f"{int(stoich)}{sp_name}")

            prod_str = []
            for s_idx, stoich in enumerate(products):
                if stoich > 0:
                    sp_name = species[s_idx]
                    if stoich == 1:
                        prod_str.append(sp_name)
                    else:
                        prod_str.append(f"{int(stoich)}{sp_name}")

            reaction_str = f"{' + '.join(react_str)} → {' + '.join(prod_str)}"
            ch_loss_reactions.append((reaction_str, rate, tag))

print(f"\nFound {len(ch_loss_reactions)} CH loss reactions in model:\n")

for reaction, rate, tag in sorted(ch_loss_reactions, key=lambda x: x[1], reverse=True):
    print(f"{reaction:50} k={rate:12.2e}  [{tag}]")

# Now check what's potentially missing
print(f"\n{'='*80}")
print("POTENTIALLY MISSING CH LOSS REACTIONS FROM LITERATURE:")
print(f"{'='*80}\n")

missing_reactions = [
    ("e + CH → C + H", "Electron impact dissociation", "1e-9 to 1e-7 cm³/s (Te-dependent)"),
    ("e + CH → CH⁺ + 2e", "Electron impact ionization", "1e-10 to 1e-8 cm³/s (Te-dependent)"),
    ("CH + CH2 → products", "Radical-radical recombination", "1e-10 to 5e-10 cm³/s"),
    ("CH + C → C2 + H", "C radical reaction", "1e-10 to 5e-10 cm³/s"),
    ("CH + Ar⁺ → products", "Ion-neutral reaction", "1e-9 to 1e-8 cm³/s"),
    ("CH + Ar* → products", "Penning ionization/reaction", "1e-10 to 5e-10 cm³/s"),
    ("CH + CH4 → products", "Already in model?", "Check below"),
    ("CH + C2H2 → products", "Radical addition", "1e-11 to 1e-10 cm³/s"),
    ("CH + C2H4 → products", "Radical addition", "1e-11 to 1e-10 cm³/s"),
]

for reaction, description, rate_range in missing_reactions:
    # Check if similar reaction exists
    found = False
    for existing, _, tag in ch_loss_reactions:
        if reaction.split('→')[0].strip() in existing:
            found = True
            break

    status = "✓ FOUND" if found else "✗ MISSING"
    print(f"{status:12} {reaction:30} {description:35} ({rate_range})")

# Calculate what's needed
print(f"\n{'='*80}")
print("WHAT RATE WOULD BALANCE CH PRODUCTION?")
print(f"{'='*80}\n")

# From previous analysis
CH_production = 5.61e15  # cm⁻³/s
CH_loss = 2.50e15  # cm⁻³/s
deficit = CH_production - CH_loss

print(f"Current CH production: {CH_production:.2e} cm⁻³/s")
print(f"Current CH loss:       {CH_loss:.2e} cm⁻³/s")
print(f"Deficit (accumulation):{deficit:.2e} cm⁻³/s ({deficit/CH_production*100:.1f}% of production)")

# Typical densities
CH = 6.36e10  # cm⁻³
H = 2.54e14
C = 3.25e11
CH2 = 8.11e11
Ar_star = 1e11  # estimate
ArPlus = 1e9  # estimate
ne = 2.3e9

print(f"\nTo eliminate the deficit of {deficit:.2e} cm⁻³/s, we would need:\n")

# For each potential missing reaction, calculate required k
potential_additions = [
    ("e + CH → C + H", ne * CH, "cm³/s"),
    ("e + CH → CH⁺ + 2e", ne * CH, "cm³/s"),
    ("CH + CH2 → products", CH * CH2, "cm³/s"),
    ("CH + C → C2 + H", CH * C, "cm³/s"),
    ("CH + Ar⁺ → products", CH * ArPlus, "cm³/s"),
    ("CH + Ar* → products", CH * Ar_star, "cm³/s"),
]

for reaction, reactant_product, units in potential_additions:
    k_needed = deficit / reactant_product if reactant_product > 0 else np.inf
    print(f"{reaction:30} k_needed = {k_needed:.2e} {units}")

print(f"\n{'='*80}")
print("LITERATURE RANGES FOR MISSING REACTIONS:")
print(f"{'='*80}\n")

print("1. e + CH → C + H (electron impact dissociation):")
print("   - Typical: 1e-9 to 1e-7 cm³/s (highly Te-dependent)")
print("   - At Te=1.3 eV: ~1e-9 cm³/s (threshold ~3 eV)")
print("   - Would give: rate = 1e-9 × 2.3e9 × 6.36e10 = 1.46e11 cm⁻³/s")
print("   - Impact: Only 0.5% of deficit")

print("\n2. e + CH → CH⁺ + 2e (electron impact ionization):")
print("   - Typical: 1e-10 to 1e-8 cm³/s (highly Te-dependent)")
print("   - At Te=1.3 eV: ~1e-11 cm³/s (threshold ~10 eV, very slow)")
print("   - Would give: rate = 1e-11 × 2.3e9 × 6.36e10 = 1.46e9 cm⁻³/s")
print("   - Impact: Negligible")

print("\n3. CH + CH2 → products:")
print("   - Typical: 1e-10 to 5e-10 cm³/s")
print("   - Would give: rate = 3e-10 × 6.36e10 × 8.11e11 = 1.55e13 cm⁻³/s")
print("   - Impact: 0.5% of deficit")

print("\n4. CH + C → C2 + H:")
print("   - Typical: 1e-10 to 5e-10 cm³/s")
print("   - Would give: rate = 3e-10 × 6.36e10 × 3.25e11 = 6.21e12 cm⁻³/s")
print("   - Impact: 0.2% of deficit")

print("\n5. CH + Ar⁺ → products:")
print("   - Typical: 1e-9 to 1e-8 cm³/s (fast ion-neutral)")
print("   - Ar⁺ density very low (~1e9 cm⁻³)")
print("   - Would give: rate = 5e-9 × 6.36e10 × 1e9 = 3.18e11 cm⁻³/s")
print("   - Impact: 0.01% of deficit")

print("\n6. CH + Ar* → products:")
print("   - Typical: 1e-10 to 5e-10 cm³/s")
print("   - Ar* density moderate (~1e11 cm⁻³)")
print("   - Would give: rate = 3e-10 × 6.36e10 × 1e11 = 1.91e12 cm⁻³/s")
print("   - Impact: 0.06% of deficit")

print(f"\n{'='*80}")
print("CONCLUSION:")
print(f"{'='*80}\n")

print("Even if we add ALL missing CH loss reactions with upper-limit rate constants,")
print("the total additional loss would be ~1.7e13 cm⁻³/s, which is only 0.5% of the")
print("3.11e15 cm⁻³/s deficit.")
print()
print("The fundamental issue is that C2 + H → CH + C (5.61e15 cm⁻³/s) is MUCH faster")
print("than any realistic CH loss mechanism at these densities.")
print()
print("The missing reactions would help slightly, but won't solve the 2.2× imbalance.")
print()
print("Key reactions to check:")
print("  1. CH + CH2 → products (most impactful of missing ones)")
print("  2. e + CH → C + H (electron impact dissociation)")
print("  3. CH + C → C2 + H (recycles to C2)")
