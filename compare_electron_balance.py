#!/usr/bin/env python3
"""
Compare electron balance between:
1. High C2H2 state (C2H2=4.09e+12, Ni/Ne=215) - UNSTABLE
2. Good plasma state (C2H2=4.81e+09, Ni/Ne=3.12) - STABLE

Find what's different!
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

def analyze_state(filename, label):
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

    # Calculate electron balance
    R = params['R']
    tags = params['tags']
    idx_e = species.index('e')

    electron_production = []
    electron_loss = []

    for i, rxn in enumerate(R):
        net_e = rxn.products[idx_e] - rxn.reactants[idx_e]
        if net_e != 0:
            rate_coeff = rxn.rate
            reactant_density = 1.0
            for j, stoich in enumerate(rxn.reactants):
                if stoich > 0:
                    reactant_density *= y0[j]**stoich

            reaction_rate = rate_coeff * reactant_density
            electron_rate = reaction_rate * net_e

            if net_e > 0:
                electron_production.append({
                    'rate': electron_rate,
                    'tag': tags[i],
                    'k': rate_coeff
                })
            else:
                electron_loss.append({
                    'rate': -electron_rate,
                    'tag': tags[i],
                    'k': rate_coeff
                })

    total_prod = sum(ep['rate'] for ep in electron_production)
    total_loss = sum(el['rate'] for el in electron_loss)

    # Top ionization reactions
    electron_production.sort(key=lambda x: x['rate'], reverse=True)
    ionization_rates = {}
    for ep in electron_production[:5]:
        ionization_rates[ep['tag']] = {'rate': ep['rate'], 'k': ep['k']}

    return {
        'label': label,
        'C2H2': data['all_densities']['C2H2'],
        'Ne': data['Ne'],
        'Ni': data['n_i_total'],
        'Ni_over_Ne': data['Ni_over_Ne'],
        'Te': data['Te'],
        'E_field': data['E_field'],
        'total_e_prod': total_prod,
        'total_e_loss': total_loss,
        'ionization_rates': ionization_rates,
        'Ar': data['all_densities']['Ar'],
        'ArStar': data['all_densities']['ArStar'],
        'CH4': data['all_densities']['CH4'],
    }

print('=' * 80)
print('COMPARING ELECTRON BALANCE: UNSTABLE vs STABLE PLASMA')
print('=' * 80)
print()

# Analyze both states
unstable = analyze_state('optimization_results_comprehensive_1e12/best_f13946912.1.json', 'UNSTABLE (High C2H2)')
stable = analyze_state('optimization_results_comprehensive_1e12/best_f70.3.json', 'STABLE (Low C2H2)')

print('STATE 1: UNSTABLE PLASMA (High C2H2, Poor Charge Balance)')
print('-' * 80)
print(f'  C2H2:     {unstable["C2H2"]:.2e} cm⁻³ (409% of target) ✓')
print(f'  Ne:       {unstable["Ne"]:.2e} cm⁻³')
print(f'  Ni:       {unstable["Ni"]:.2e} cm⁻³')
print(f'  Ni/Ne:    {unstable["Ni_over_Ne"]:.1f} ✗ (way too high!)')
print(f'  Te:       {unstable["Te"]:.2f} eV')
print(f'  E-field:  {unstable["E_field"]:.1f} V/cm')
print()
print(f'  Electron production: {unstable["total_e_prod"]:.2e} cm⁻³/s')
print(f'  Electron loss:       {unstable["total_e_loss"]:.2e} cm⁻³/s')
print(f'  Production/Loss:     {unstable["total_e_prod"]/unstable["total_e_loss"]:.1f}×')
print()

print('STATE 2: STABLE PLASMA (Low C2H2, Good Charge Balance)')
print('-' * 80)
print(f'  C2H2:     {stable["C2H2"]:.2e} cm⁻³ (0.5% of target) ✗')
print(f'  Ne:       {stable["Ne"]:.2e} cm⁻³')
print(f'  Ni:       {stable["Ni"]:.2e} cm⁻³')
print(f'  Ni/Ne:    {stable["Ni_over_Ne"]:.1f} ✓ (good!)')
print(f'  Te:       {stable["Te"]:.2f} eV')
print(f'  E-field:  {stable["E_field"]:.1f} V/cm')
print()
print(f'  Electron production: {stable["total_e_prod"]:.2e} cm⁻³/s')
print(f'  Electron loss:       {stable["total_e_loss"]:.2e} cm⁻³/s')
print(f'  Production/Loss:     {stable["total_e_prod"]/stable["total_e_loss"]:.1f}×')
print()

print('=' * 80)
print('KEY DIFFERENCES')
print('=' * 80)
print()

print('Plasma Parameters:')
print(f'  Te:      {unstable["Te"]:.3f} eV (unstable) vs {stable["Te"]:.3f} eV (stable) - Δ = {abs(unstable["Te"]-stable["Te"]):.3f} eV')
print(f'  E-field: {unstable["E_field"]:.1f} V/cm (unstable) vs {stable["E_field"]:.1f} V/cm (stable) - Δ = {abs(unstable["E_field"]-stable["E_field"]):.1f} V/cm')
print()

print('Neutral Densities:')
print(f'  Ar:      {unstable["Ar"]:.2e} (unstable) vs {stable["Ar"]:.2e} (stable)')
print(f'  Ar*:     {unstable["ArStar"]:.2e} (unstable) vs {stable["ArStar"]:.2e} (stable) - {unstable["ArStar"]/stable["ArStar"]:.1f}× difference!')
print(f'  CH4:     {unstable["CH4"]:.2e} (unstable) vs {stable["CH4"]:.2e} (stable)')
print()

print('Ionization Rate Constants (top reaction: e + Ar → Ar+ + 2e):')
for tag in unstable['ionization_rates']:
    if tag in stable['ionization_rates']:
        k_unstable = unstable['ionization_rates'][tag]['k']
        k_stable = stable['ionization_rates'][tag]['k']
        print(f'  {tag}:')
        print(f'    Unstable: k = {k_unstable:.2e}')
        print(f'    Stable:   k = {k_stable:.2e}')
        print(f'    Ratio:    {k_unstable/k_stable:.1f}×')
        print()

print('=' * 80)
print('HYPOTHESIS')
print('=' * 80)
print()
print('The unstable state has EXTREMELY high ionization rates (e + Ar → Ar+ + 2e)')
print('due to tuned rate constants that are {:.0f}× higher than stable state.'.format(
    unstable['ionization_rates']['e_Ar_ArPlus_cm3_2_3']['k'] / stable['ionization_rates']['e_Ar_ArPlus_cm3_2_3']['k']
))
print()
print('This creates a runaway situation:')
print('  1. High ionization rate → More ions created')
print('  2. More ions → Higher Ni')
print('  3. Electron-ion recombination can\'t keep up → Low Ne')
print('  4. Result: Ni/Ne >> 7 (unstable plasma)')
print()
print('SOLUTION: The ionization rate constants may need TIGHTER bounds')
print('to prevent optimizer from creating unphysical high-ionization states.')
