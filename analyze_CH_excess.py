"""
Analyze why CH is 6355% of target (63× too high)

Question: What's producing so much CH and why isn't it being consumed?
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

def pressure_to_density(pressure_mTorr, T_K=300):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

targets = {'H': 2.52e14, 'CH': 1.0e9, 'C2': 5.6e11}

def analyze_CH_chemistry(ne_value, use_new_pathways=False):
    """Run simulation and analyze CH production/loss"""
    print(f"\n{'='*80}")
    print(f"ANALYZING CH CHEMISTRY")
    print(f"{'='*80}")
    print(f"ne: {ne_value:.2e} cm⁻³")
    print(f"New CH3 pathways: {'YES' if use_new_pathways else 'NO'}")

    P = 500.0
    n_total = pressure_to_density(P)

    params = {
        'P': P,
        'Te': baseline['Te'],
        'ne': ne_value,
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

    # Apply baseline tuned rates
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    # Apply physically realistic multipliers
    rate_mults = {
        'e_CH4_CH3_HMinus_cm3_8_1': 10.0,
        'ArStar_CH4_CH3_H_cm3_3_1': 1.4,
        'e_CH4Plus_CH3_H_cm3_6_4': 1.5,
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
    }

    # Suppress new pathways if not wanted
    if not use_new_pathways:
        new_pathway_keys = [
            'e_C2H4_CH3_CH_cm3_7_66',
            'e_C2H6_CH3_CH3_cm3_7_67',
            'e_C2H5_CH3_CH2_cm3_7_68',
            'ArStar_C2H4_CH3_CH_cm3_7_69',
            'ArStar_C2H6_CH3_CH3_cm3_7_70',
            'H_C2H5_CH3_CH2_cm3_7_71',
            'CH2_CH2_CH3_CH_cm3_7_72',
            'C2H5Plus_e_CH3_CH2_cm3_7_73',
            'CH2_H_M_CH3_M_cm6_7_74',
            'ArPlus_CH4_CH3Plus_ArH_cm3_7_75',
        ]
        for key in new_pathway_keys:
            if key in k:
                k[key] = 0.0

    for rate_name, mult in rate_mults.items():
        if rate_name in k:
            k[rate_name] *= mult

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    species = params['species']
    y0 = np.ones(len(species)) * 1e3
    y0[species.index('Ar')] = n_total * 0.85
    y0[species.index('CH4')] = n_total * 0.15
    y0[species.index('e')] = ne_value

    try:
        ode_func = PlasmaODE_Optimized(params)
        sol = solve_ivp(ode_func, (0, 500), y0, method='BDF', rtol=1e-7, atol=1e-9, max_step=1.0)

        if not sol.success:
            print("✗ Integration failed!")
            return None

        y_final = sol.y[:, -1]

        # Get final densities
        densities = {}
        for sp in ['H', 'CH', 'C2', 'CH2', 'CH3', 'C', 'C2H2', 'C2H', 'C2H4', 'CHPlus']:
            densities[sp] = y_final[species.index(sp)]

        print(f"\nFinal densities:")
        print(f"  H:    {densities['H']:.2e} ({densities['H']/targets['H']*100:5.1f}%)")
        print(f"  CH:   {densities['CH']:.2e} ({densities['CH']/targets['CH']*100:5.1f}%)")
        print(f"  C2:   {densities['C2']:.2e} ({densities['C2']/targets['C2']*100:5.2f}%)")
        print(f"  CH2:  {densities['CH2']:.2e}")
        print(f"  CH3:  {densities['CH3']:.2e}")
        print(f"  C:    {densities['C']:.2e}")
        print(f"  C2H2: {densities['C2H2']:.2e}")

        # Calculate CH production rates
        print(f"\n{'='*80}")
        print("CH PRODUCTION PATHWAYS:")
        print(f"{'='*80}")

        productions = []

        # C2 + H → CH + C (MAJOR SUSPECT!)
        if 'C2_H_CH_C_cm3_7_6' in k:
            rate_C2_H = k['C2_H_CH_C_cm3_7_6'] * densities['C2'] * densities['H']
            productions.append(('C2 + H → CH + C', rate_C2_H, k['C2_H_CH_C_cm3_7_6']))
            print(f"C2 + H → CH + C:")
            print(f"  k = {k['C2_H_CH_C_cm3_7_6']:.2e} cm³/s")
            print(f"  [C2] = {densities['C2']:.2e}, [H] = {densities['H']:.2e}")
            print(f"  Rate = {rate_C2_H:.2e} cm⁻³/s")

        # C + CH2 → CH + CH
        if 'C_CH2_CH_CH_cm3_7_8' in k:
            rate_C_CH2 = k['C_CH2_CH_CH_cm3_7_8'] * densities['C'] * densities['CH2']
            productions.append(('C + CH2 → CH + CH', rate_C_CH2, k['C_CH2_CH_CH_cm3_7_8']))
            print(f"\nC + CH2 → CH + CH:")
            print(f"  k = {k['C_CH2_CH_CH_cm3_7_8']:.2e} cm³/s")
            print(f"  Rate = {rate_C_CH2:.2e} cm⁻³/s")

        # H + C2H → CH + CH
        if 'H_C2H_CH_CH_cm3_7_7' in k:
            rate_H_C2H = k['H_C2H_CH_CH_cm3_7_7'] * densities['H'] * densities['C2H']
            productions.append(('H + C2H → CH + CH', rate_H_C2H, k['H_C2H_CH_CH_cm3_7_7']))
            print(f"\nH + C2H → CH + CH:")
            print(f"  k = {k['H_C2H_CH_CH_cm3_7_7']:.2e} cm³/s")
            print(f"  Rate = {rate_H_C2H:.2e} cm⁻³/s")

        # e + C2H4 → CH3 + CH (new pathway if enabled)
        if use_new_pathways and 'e_C2H4_CH3_CH_cm3_7_66' in k:
            rate_e_C2H4 = k['e_C2H4_CH3_CH_cm3_7_66'] * ne_value * densities['C2H4']
            productions.append(('e + C2H4 → CH3 + CH', rate_e_C2H4, k['e_C2H4_CH3_CH_cm3_7_66']))
            print(f"\ne + C2H4 → CH3 + CH:")
            print(f"  k = {k['e_C2H4_CH3_CH_cm3_7_66']:.2e} cm³/s")
            print(f"  Rate = {rate_e_C2H4:.2e} cm⁻³/s")

        # CH2 + CH2 → CH3 + CH (new pathway if enabled)
        if use_new_pathways and 'CH2_CH2_CH3_CH_cm3_7_72' in k:
            rate_CH2_CH2 = k['CH2_CH2_CH3_CH_cm3_7_72'] * densities['CH2']**2
            productions.append(('CH2 + CH2 → CH3 + CH', rate_CH2_CH2, k['CH2_CH2_CH3_CH_cm3_7_72']))
            print(f"\nCH2 + CH2 → CH3 + CH:")
            print(f"  k = {k['CH2_CH2_CH3_CH_cm3_7_72']:.2e} cm³/s")
            print(f"  Rate = {rate_CH2_CH2:.2e} cm⁻³/s")

        # Ar* + C2H4 → CH3 + CH (new pathway if enabled)
        if use_new_pathways and 'ArStar_C2H4_CH3_CH_cm3_7_69' in k:
            ArStar = y_final[species.index('ArStar')]
            rate_ArStar_C2H4 = k['ArStar_C2H4_CH3_CH_cm3_7_69'] * ArStar * densities['C2H4']
            productions.append(('Ar* + C2H4 → CH3 + CH', rate_ArStar_C2H4, k['ArStar_C2H4_CH3_CH_cm3_7_69']))
            print(f"\nAr* + C2H4 → CH3 + CH:")
            print(f"  k = {k['ArStar_C2H4_CH3_CH_cm3_7_69']:.2e} cm³/s")
            print(f"  Rate = {rate_ArStar_C2H4:.2e} cm⁻³/s")

        total_production = sum(rate for _, rate, _ in productions)

        # Calculate CH loss rates
        print(f"\n{'='*80}")
        print("CH LOSS PATHWAYS:")
        print(f"{'='*80}")

        losses = []

        # CH + CH → C2 + H2 (multiple possible keys)
        for key in ['CH_CH_C2_H2_cm3_5_4', 'CH_CH_C2_H2_cm3_7_9', 'CH_CH_C2_H2_cm3_7_44']:
            if key in k and k[key] > 0:
                rate_CH_CH = k[key] * densities['CH']**2
                losses.append((f'CH + CH → C2 + H2 ({key})', rate_CH_CH, k[key]))
                print(f"CH + CH → C2 + H2 ({key}):")
                print(f"  k = {k[key]:.2e} cm³/s")
                print(f"  [CH]² = {densities['CH']**2:.2e}")
                print(f"  Rate = {rate_CH_CH:.2e} cm⁻³/s")

        # CH + H → C + H2
        for key in ['CH_H_C_H2_cm3_7_3', 'CH_H_C_H2_cm3_7_10']:
            if key in k and k[key] > 0:
                rate_CH_H = k[key] * densities['CH'] * densities['H']
                losses.append((f'CH + H → C + H2 ({key})', rate_CH_H, k[key]))
                print(f"\nCH + H → C + H2 ({key}):")
                print(f"  k = {k[key]:.2e} cm³/s")
                print(f"  Rate = {rate_CH_H:.2e} cm⁻³/s")

        # Wall sticking
        for key in ['stick_CH_9_3', 'stick_CH_9_5']:
            if key in k and k[key] > 0:
                rate_stick = k[key] * densities['CH']
                losses.append((f'CH wall sticking ({key})', rate_stick, k[key]))
                print(f"\nCH wall sticking ({key}):")
                print(f"  k = {k[key]:.2e} s⁻¹")
                print(f"  Rate = {rate_stick:.2e} cm⁻³/s")

        # CH + CH3 → products (multiple channels)
        for key in ['CH_CH3_C2H4_cm3_7_5', 'CH_CH3_C2H3_H_cm3_7_10', 'CH_CH3_C2H3_H_cm3_7_11', 'CH_CH3_C2H2_H2_cm3_7_23']:
            if key in k and k[key] > 0:
                rate_CH_CH3 = k[key] * densities['CH'] * densities['CH3']
                losses.append((f'CH + CH3 ({key})', rate_CH_CH3, k[key]))
                print(f"\nCH + CH3 ({key}):")
                print(f"  k = {k[key]:.2e} cm³/s")
                print(f"  Rate = {rate_CH_CH3:.2e} cm⁻³/s")

        total_loss = sum(rate for _, rate, _ in losses)

        # Summary
        print(f"\n{'='*80}")
        print("SUMMARY:")
        print(f"{'='*80}")
        print(f"\nTotal CH production: {total_production:.2e} cm⁻³/s")
        print(f"Total CH loss:       {total_loss:.2e} cm⁻³/s")
        print(f"Net rate:            {total_production - total_loss:.2e} cm⁻³/s")

        print(f"\nCH production breakdown:")
        for name, rate, k_val in sorted(productions, key=lambda x: x[1], reverse=True):
            pct = (rate / total_production * 100) if total_production > 0 else 0
            print(f"  {name:30} {rate:>12.2e} ({pct:5.1f}%)")

        print(f"\nCH loss breakdown:")
        for name, rate, k_val in sorted(losses, key=lambda x: x[1], reverse=True):
            pct = (rate / total_loss * 100) if total_loss > 0 else 0
            print(f"  {name:30} {rate:>12.2e} ({pct:5.1f}%)")

        print(f"\n{'='*80}")
        print("ANALYSIS:")
        print(f"{'='*80}")

        # Find dominant production
        if productions:
            dom_prod_name, dom_prod_rate, _ = max(productions, key=lambda x: x[1])
            print(f"\nDominant CH production: {dom_prod_name}")
            print(f"  Rate: {dom_prod_rate:.2e} cm⁻³/s ({dom_prod_rate/total_production*100:.1f}%)")

        # Find dominant loss
        if losses:
            dom_loss_name, dom_loss_rate, _ = max(losses, key=lambda x: x[1])
            print(f"\nDominant CH loss: {dom_loss_name}")
            print(f"  Rate: {dom_loss_rate:.2e} cm⁻³/s ({dom_loss_rate/total_loss*100:.1f}%)")

        # Imbalance
        if total_production > 0 and total_loss > 0:
            imbalance = abs(total_production - total_loss) / max(total_production, total_loss) * 100
            print(f"\nProduction/Loss imbalance: {imbalance:.1f}%")

            if total_production > total_loss:
                print(f"  ⚠️  CH ACCUMULATING (production {total_production/total_loss:.1f}× loss)")
            else:
                print(f"  ✓ CH DEPLETING (loss {total_loss/total_production:.1f}× production)")

        return {
            'CH': densities['CH'],
            'total_production': total_production,
            'total_loss': total_loss,
            'productions': productions,
            'losses': losses,
        }

    except Exception as e:
        print(f"✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("ANALYZE WHY CH IS 6355% OF TARGET (63× TOO HIGH)")
    print("="*80)

    ne_user = 2.3e9

    # Analyze with corrected ne (no new pathways)
    result = analyze_CH_chemistry(ne_user, use_new_pathways=False)

    if result:
        print(f"\n{'='*80}")
        print("CONCLUSION:")
        print(f"{'='*80}")

        if result['productions']:
            dom_prod = max(result['productions'], key=lambda x: x[1])
            print(f"\nCH is {result['CH']/targets['CH']:.1f}× too high because:")
            print(f"  1. {dom_prod[0]} is producing it at {dom_prod[1]:.2e} cm⁻³/s")
            print(f"  2. This is {dom_prod[1]/result['total_loss']:.1f}× faster than total loss")
            print(f"\nTo reduce CH, we need to either:")
            print(f"  - Reduce the rate of {dom_prod[0]}")
            print(f"  - Increase CH loss mechanisms")
            print(f"  - Accept that this is the equilibrium chemistry")
