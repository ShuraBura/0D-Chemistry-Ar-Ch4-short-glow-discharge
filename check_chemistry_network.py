"""
Comprehensive chemistry network audit
"""

import json
from define_rates import define_rates
from build_reactions import build_reactions

with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

params = {
    'P': 500.0,
    'Te': 1.0,
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
params['k'] = k
R, tags = build_reactions(params)

species = params['species']

print("="*80)
print("CHEMISTRY NETWORK AUDIT")
print("="*80)

# Count reactions by type
print(f"\nTotal reactions: {len(R)}")
print(f"Total rate constants: {len(k)}")

# Analyze for H, CH, C2
for target in ['H', 'CH', 'C2']:
    idx = species.index(target)
    production = sum(1 for r in R if r[idx] > 0)
    loss = sum(1 for r in R if r[idx] < 0)
    print(f"\n{target}: {production} production, {loss} loss reactions")

# List key reactions involving our target species
print(f"\n{'='*80}")
print("KEY REACTIONS FOR TARGET SPECIES")
print(f"{'='*80}")

def find_reactions_with_species(R, tags, species, target_species, max_show=15):
    """Find reactions involving target species"""
    results = []
    target_idx = species.index(target_species)
    
    for i, (r, tag) in enumerate(zip(R, tags)):
        if abs(r[target_idx]) > 1e-10:
            # Build reaction string
            reactants = []
            products = []
            for j, coef in enumerate(r[:len(species)]):
                if coef < 0:
                    reactants.append(f"{species[j]}" if abs(coef) == 1 else f"{int(abs(coef))}{species[j]}")
                elif coef > 0:
                    products.append(f"{species[j]}" if coef == 1 else f"{int(coef)}{species[j]}")
            
            rxn = f"{' + '.join(reactants)} → {' + '.join(products)}"
            rate = r[len(species)]
            sign = "PROD" if r[target_idx] > 0 else "LOSS"
            results.append((rxn, rate, sign, tag))
    
    return results

for target in ['H', 'CH', 'C2']:
    print(f"\n{target} reactions:")
    reactions = find_reactions_with_species(R, tags, species, target, max_show=20)
    
    prod = [(r, k, t) for r, k, s, t in reactions if s == "PROD"]
    loss_rxns = [(r, k, t) for r, k, s, t in reactions if s == "LOSS"]
    
    print(f"  Production (showing top 10 by rate):")
    for rxn, rate, tag in sorted(prod, key=lambda x: x[1], reverse=True)[:10]:
        print(f"    k={rate:.2e}: {rxn}")
    
    print(f"  Loss (showing top 10 by rate):")
    for rxn, rate, tag in sorted(loss_rxns, key=lambda x: x[1], reverse=True)[:10]:
        print(f"    k={rate:.2e}: {rxn}")

print(f"\n{'='*80}")
print("POTENTIALLY MISSING REACTIONS TO CHECK:")
print(f"{'='*80}")

missing = """
1. H + H + M → H2 + M (three-body recombination)
   - Critical for H loss at high H densities
   - k ~ 1e-32 cm^6/s at room temp
   
2. C + C + M → C2 + M (three-body C2 formation)
   - Alternative C2 production pathway
   - k ~ 1e-33 cm^6/s
   
3. CH + H2 → CH2 + H
   - CH loss via H2 abstraction
   - k ~ 1e-12 to 1e-11 cm^3/s
   
4. CH + CH4 → C2H4 + H
   - CH coupling with CH4
   - k ~ 1e-11 cm^3/s at 400K
   
5. C2 + CH4 → products  
   - C2 reaction with CH4
   - Check if present
   
6. H wall recombination: H + wall → 1/2 H2
   - Heterogeneous H loss (γ ~ 0.01-0.1)
   
7. C wall sticking: C + wall → loss
   - Heterogeneous C loss

8. CH2 + H2 → CH3 + H
   - CH2 hydrogenation
"""

print(missing)

print(f"\n{'='*80}")
print("ACTION ITEMS:")
print(f"{'='*80}")
print("""
1. Check define_rates.py for each potentially missing reaction
2. For any missing, search literature for rate constants
3. Check conditions at which rate constants were measured
4. Consider whether H_drift_gain represents physical H flux from cathode
""")

