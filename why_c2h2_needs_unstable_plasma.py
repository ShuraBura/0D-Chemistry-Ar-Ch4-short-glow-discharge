#!/usr/bin/env python3
"""
CRITICAL ANALYSIS: Why can't C2H2 accumulate with stable plasma?

Compare C2H2 chemistry between:
1. Stable plasma: Ni/Ne=2.95, C2H2=3.54e+09 (LOW)
2. Unstable plasma: Ni/Ne=215, C2H2=4.09e+12 (HIGH)

Find the mechanism linking charge balance to C2H2 accumulation!
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

def pressure_to_density(pressure_mTorr, T_K=400):
    kB = 1.38064852e-23
    Torr_to_Pa = 133.322
    P_Pa = pressure_mTorr * 1e-3 * Torr_to_Pa
    n_m3 = P_Pa / (kB * T_K)
    return n_m3 * 1e-6

def analyze_c2h2_balance(filename, label):
    with open(filename, 'r') as f:
        data = json.load(f)

    params = {
        'P': 500.0, 'Te': data['Te'], 'ne': data['Ne'], 'E_field': data['E_field'],
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

    k = define_rates(params)
    for rate_name, rate_val in data.get('rate_values', {}).items():
        if rate_name in k:
            k[rate_name] = rate_val

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    species = params['species']
    y0 = np.array([data['all_densities'].get(s, 1e3) for s in species])

    # Analyze C2H2 balance
    R = params['R']
    tags = params['tags']
    idx_C2H2 = species.index('C2H2')

    c2h2_production = []
    c2h2_loss = []

    for i, rxn in enumerate(R):
        net_c2h2 = rxn.products[idx_C2H2] - rxn.reactants[idx_C2H2]
        if net_c2h2 != 0:
            rate_coeff = rxn.rate
            reactant_density = 1.0
            for j, stoich in enumerate(rxn.reactants):
                if stoich > 0:
                    reactant_density *= y0[j]**stoich

            reaction_rate = rate_coeff * reactant_density
            c2h2_rate = reaction_rate * net_c2h2

            if net_c2h2 > 0:
                c2h2_production.append({
                    'rate': c2h2_rate,
                    'tag': tags[i],
                    'k': rate_coeff,
                    'rxn': rxn
                })
            else:
                c2h2_loss.append({
                    'rate': -c2h2_rate,
                    'tag': tags[i],
                    'k': rate_coeff,
                    'rxn': rxn
                })

    total_prod = sum(cp['rate'] for cp in c2h2_production)
    total_loss = sum(cl['rate'] for cl in c2h2_loss)

    c2h2_production.sort(key=lambda x: x['rate'], reverse=True)
    c2h2_loss.sort(key=lambda x: x['rate'], reverse=True)

    return {
        'label': label,
        'C2H2': data['all_densities']['C2H2'],
        'Ni_over_Ne': data['Ni_over_Ne'],
        'Te': data['Te'],
        'E_field': data['E_field'],
        'Ne': data['Ne'],
        'H': data['all_densities']['H'],
        'CH3': data['all_densities']['CH3'],
        'total_prod': total_prod,
        'total_loss': total_loss,
        'top_prod': c2h2_production[:3],
        'top_loss': c2h2_loss[:3],
        'species_densities': {sp: y0[species.index(sp)] for sp in ['e', 'H', 'CH3', 'C2H4']},
    }

print('=' * 80)
print('WHY CAN\'T C2H2 ACCUMULATE WITH STABLE PLASMA?')
print('=' * 80)
print()

# Analyze both states
stable = analyze_c2h2_balance('optimization_results_fixed_ionization/best_f22.8.json', 'STABLE')
unstable = analyze_c2h2_balance('optimization_results_comprehensive_1e12/best_f13946912.1.json', 'UNSTABLE')

print('STATE 1: STABLE PLASMA (Correct Physics)')
print('-' * 80)
print(f'  C2H2:     {stable["C2H2"]:.2e} cm⁻³ (LOW!)')
print(f'  Ni/Ne:    {stable["Ni_over_Ne"]:.2f} ✓ (stable)')
print(f'  Te:       {stable["Te"]:.2f} eV')
print(f'  E-field:  {stable["E_field"]:.1f} V/cm')
print(f'  Ne:       {stable["Ne"]:.2e} cm⁻³')
print(f'  H:        {stable["H"]:.2e} cm⁻³')
print(f'  CH3:      {stable["CH3"]:.2e} cm⁻³')
print()
print(f'  C2H2 production: {stable["total_prod"]:.2e} cm⁻³/s')
print(f'  C2H2 loss:       {stable["total_loss"]:.2e} cm⁻³/s')
print(f'  Production/Loss: {stable["total_prod"]/stable["total_loss"]:.2f}×')
print()

print('STATE 2: UNSTABLE PLASMA (Unphysical Ionization)')
print('-' * 80)
print(f'  C2H2:     {unstable["C2H2"]:.2e} cm⁻³ (HIGH!)')
print(f'  Ni/Ne:    {unstable["Ni_over_Ne"]:.1f} ✗ (unstable)')
print(f'  Te:       {unstable["Te"]:.2f} eV')
print(f'  E-field:  {unstable["E_field"]:.1f} V/cm')
print(f'  Ne:       {unstable["Ne"]:.2e} cm⁻³')
print(f'  H:        {unstable["H"]:.2e} cm⁻³')
print(f'  CH3:      {unstable["CH3"]:.2e} cm⁻³')
print()
print(f'  C2H2 production: {unstable["total_prod"]:.2e} cm⁻³/s')
print(f'  C2H2 loss:       {unstable["total_loss"]:.2e} cm⁻³/s')
print(f'  Production/Loss: {unstable["total_prod"]/unstable["total_loss"]:.2f}×')
print()

print('=' * 80)
print('KEY COMPARISON')
print('=' * 80)
print()

print('C2H2 Production:')
print(f'  Stable:   {stable["total_prod"]:.2e} cm⁻³/s')
print(f'  Unstable: {unstable["total_prod"]:.2e} cm⁻³/s')
print(f'  Ratio:    {unstable["total_prod"]/stable["total_prod"]:.1f}× more in unstable!')
print()

print('C2H2 Loss:')
print(f'  Stable:   {stable["total_loss"]:.2e} cm⁻³/s')
print(f'  Unstable: {unstable["total_loss"]:.2e} cm⁻³/s')
print(f'  Ratio:    {unstable["total_loss"]/stable["total_loss"]:.1f}× more in unstable')
print()

print('Species Densities:')
print(f'  H:   {stable["H"]:.2e} (stable) vs {unstable["H"]:.2e} (unstable) - {unstable["H"]/stable["H"]:.1f}× higher!')
print(f'  CH3: {stable["CH3"]:.2e} (stable) vs {unstable["CH3"]:.2e} (unstable) - {unstable["CH3"]/stable["CH3"]:.1f}× higher!')
print()

print('=' * 80)
print('THE MECHANISM:')
print('=' * 80)
print()

# Calculate ratio of C2H2 levels
c2h2_ratio = unstable['C2H2'] / stable['C2H2']
h_ratio = unstable['H'] / stable['H']
ch3_ratio = unstable['CH3'] / stable['CH3']

print(f'Unstable state has:')
print(f'  - C2H2: {c2h2_ratio:.0f}× higher')
print(f'  - H:    {h_ratio:.1f}× higher')
print(f'  - CH3:  {ch3_ratio:.1f}× higher (KEY PRECURSOR for C2H2!)')
print()

print('Top C2H2 producer is typically 2CH3 → C2H2 + 2H2')
print(f'Production rate ∝ [CH3]²')
print()
print(f'If CH3 is {ch3_ratio:.1f}× higher, production should be {ch3_ratio**2:.1f}× higher')
print(f'Actual C2H2 production ratio: {unstable["total_prod"]/stable["total_prod"]:.1f}×')
print()

print('CONCLUSION:')
print(f'The unstable (high ionization) state has much higher CH3 density!')
print(f'This drives massive C2H2 production via 2CH3 → C2H2')
print()
print('When ionization rates are FIXED at correct values:')
print('  → Lower CH3 density')
print('  → Lower C2H2 production (∝ [CH3]²)')
print('  → C2H2 accumulates to only {:.2e} instead of {:.2e}'.format(stable["C2H2"], unstable["C2H2"]))
print()
print('THE COUPLING: High ionization → High CH3 → High C2H2 production → High C2H2')
