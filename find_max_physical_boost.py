"""
Find maximum CH3 boost that keeps ALL rates within literature ranges

Question: What's the highest multiplier we can use while staying physical?
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

# Literature ranges
lit_ranges = {
    'e_CH4_CH3_HMinus_cm3_8_1': (1e-13, 1e-8, 'e-impact'),
    'ArStar_CH4_CH3_H_cm3_3_1': (1e-11, 1e-9, 'metastable'),
    'e_CH4Plus_CH3_H_cm3_6_4': (1e-8, 1e-6, 'dissoc. recomb.'),
}

def find_max_multiplier():
    """Find max multiplier that keeps all rates in lit ranges"""
    P = 500.0
    n_total = pressure_to_density(P)
    ne_frac = baseline['Ne'] / pressure_to_density(500.0)
    ne = ne_frac * n_total

    params = {
        'P': P,
        'Te': baseline['Te'],
        'ne': ne,
        'E_field': baseline['E_field'],
        'L_discharge': 0.45,
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
        },
        'species': [],
    }

    k = define_rates(params)

    # Apply baseline tuning
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    print("="*80)
    print("FINDING MAXIMUM PHYSICAL CH3 BOOST")
    print("="*80)

    print("\nBaseline (tuned) rate constants:")
    for rate_key, (min_lit, max_lit, rxn_type) in lit_ranges.items():
        tuned_k = k[rate_key]
        print(f"  {rate_key}")
        print(f"    Current: {tuned_k:.2e} cm³/s")
        print(f"    Lit range: {min_lit:.0e} to {max_lit:.0e} cm³/s ({rxn_type})")
        max_possible_mult = max_lit / tuned_k
        print(f"    Max multiplier before exceeding: {max_possible_mult:.1f}×")
        print()

    # Find the limiting reaction (smallest max multiplier)
    max_multipliers = {}
    for rate_key, (min_lit, max_lit, rxn_type) in lit_ranges.items():
        tuned_k = k[rate_key]
        max_multipliers[rate_key] = max_lit / tuned_k

    limiting_reaction = min(max_multipliers, key=max_multipliers.get)
    max_physical_mult = max_multipliers[limiting_reaction]

    print("="*80)
    print("LIMITING REACTION:")
    print("="*80)
    print(f"\nReaction: {limiting_reaction}")
    print(f"Max physical multiplier: {max_physical_mult:.1f}×")
    print(f"This is the bottleneck keeping us from higher boosts")

    return max_physical_mult

def test_physical_boost(multiplier):
    """Test C2 result with physically constrained boost"""
    P = 500.0
    n_total = pressure_to_density(P)
    ne_frac = baseline['Ne'] / pressure_to_density(500.0)
    ne = ne_frac * n_total

    params = {
        'P': P,
        'Te': baseline['Te'] * 1.1,  # Optimal Te from earlier
        'ne': ne,
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
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    # Apply physically realistic boost
    rate_mults = {
        'e_CH4_CH3_HMinus_cm3_8_1': multiplier,
        'ArStar_CH4_CH3_H_cm3_3_1': multiplier,
        'e_CH4Plus_CH3_H_cm3_6_4': multiplier,
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
    }

    for rate_name, mult in rate_mults.items():
        if rate_name in k:
            k[rate_name] *= mult

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    species = params['species']
    y0 = np.ones(len(species)) * 1e3
    y0[species.index('Ar')] = n_total * 0.85
    y0[species.index('CH4')] = n_total * 0.15
    y0[species.index('e')] = ne

    try:
        ode_func = PlasmaODE_Optimized(params)
        sol = solve_ivp(ode_func, (0, 500), y0, method='BDF', rtol=1e-7, atol=1e-9, max_step=1.0)

        if not sol.success:
            return None

        y_final = sol.y[:, -1]
        H_final = y_final[species.index('H')]
        CH_final = y_final[species.index('CH')]
        C2_final = y_final[species.index('C2')]

        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus', 'CHPlus', 'H2Plus']
        n_i_total = sum(y_final[species.index(ion)] for ion in ions)
        ne_final = y_final[species.index('e')]
        Ni_over_Ne = n_i_total / ne_final if ne_final > 0 else 0

        C2_improvement = C2_final / baseline['target_densities']['C2']

        return {
            'C2': C2_final,
            'H': H_final,
            'CH': CH_final,
            'Ni_over_Ne': Ni_over_Ne,
            'improvement': C2_improvement,
            'stable': Ni_over_Ne < 10
        }

    except Exception as e:
        return None

if __name__ == '__main__':
    # Find maximum physical multiplier
    max_mult = find_max_multiplier()

    print("\n" + "="*80)
    print("TEST WITH PHYSICALLY REALISTIC BOOST")
    print("="*80)

    # Test at max physical multiplier
    print(f"\nTesting at {max_mult:.1f}× CH3 boost (max physical):")
    result_physical = test_physical_boost(max_mult)

    if result_physical:
        c2_pct = result_physical['C2'] / targets['C2'] * 100
        print(f"  C2: {result_physical['C2']:.2e} ({c2_pct:.2f}%)")
        print(f"  H:  {result_physical['H']:.2e} ({result_physical['H']/targets['H']*100:.1f}%)")
        print(f"  CH: {result_physical['CH']:.2e} ({result_physical['CH']/targets['CH']*100:.1f}%)")
        print(f"  Ni/Ne: {result_physical['Ni_over_Ne']:.2f}")
        print(f"  Improvement: {result_physical['improvement']:.1f}×")
        print(f"  Status: {'✓ STABLE' if result_physical['stable'] else '✗ RUNAWAY'}")

    # Compare to 200× (non-physical)
    print(f"\nTesting at 200× CH3 boost (NON-PHYSICAL, for comparison):")
    result_200x = test_physical_boost(200.0)

    if result_200x:
        c2_pct = result_200x['C2'] / targets['C2'] * 100
        print(f"  C2: {result_200x['C2']:.2e} ({c2_pct:.2f}%)")
        print(f"  Improvement: {result_200x['improvement']:.1f}×")
        print(f"  Status: {'✓ STABLE' if result_200x['stable'] else '✗ RUNAWAY'}")

    print("\n" + "="*80)
    print("RECOMMENDATION:")
    print("="*80)
    print(f"\nUse {max_mult:.1f}× boost (physically realistic)")
    print(f"NOT 200× boost (exceeds literature ranges)")
    print("\nIf you need more CH3, should instead:")
    print("  1. Add missing production pathways explicitly")
    print("  2. Increase ne or ArStar densities (physical parameters)")
    print("  3. Include surface chemistry")
    print("  4. Consider experimental uncertainties in baseline rates")
    print("="*80)
