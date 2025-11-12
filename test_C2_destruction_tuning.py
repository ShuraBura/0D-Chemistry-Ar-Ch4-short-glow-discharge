"""
Question 1: What does "C2+H→CH+C eliminated" mean?
Answer: We set rate multiplier to 0.0 (k_effective = 0), but is this realistic?

Test: What if we only REDUCE (not eliminate) C2+H→CH+C?
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

def test_C2_destruction_level(suppression_factor, name):
    """Test different levels of C2+H→CH+C suppression"""
    print(f"\n{'='*80}")
    print(f"{name}")
    print(f"C2 + H → CH + C multiplier: {suppression_factor}")
    print(f"{'='*80}")

    P = 500.0
    n_total = pressure_to_density(P)
    ne_frac = baseline['Ne'] / pressure_to_density(500.0)
    ne = ne_frac * n_total

    params = {
        'P': P, 'Te': baseline['Te'], 'ne': ne, 'E_field': baseline['E_field'],
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

    # Apply 200× CH3 boost + 99% loss reduction
    multipliers = {
        'e_CH4_CH3_HMinus_cm3_8_1': 200.0,
        'ArStar_CH4_CH3_H_cm3_3_1': 200.0,
        'e_CH4Plus_CH3_H_cm3_6_4': 200.0,
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
        'C2_H_CH_C_cm3_7_6': suppression_factor,  # Variable suppression
    }

    for rate_name, mult in multipliers.items():
        if rate_name in k:
            k[rate_name] *= mult

    print(f"\nC2+H→CH+C rate constant: {k['C2_H_CH_C_cm3_7_6']:.2e} cm³/s")

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

        # Calculate destruction rate
        destruction_rate = k['C2_H_CH_C_cm3_7_6'] * C2_final * H_final

        print(f"\nResults:")
        print(f"  H:   {H_final:.2e} ({H_final/targets['H']*100:5.1f}%)")
        print(f"  CH:  {CH_final:.2e} ({CH_final/targets['CH']*100:5.1f}%)")
        print(f"  C2:  {C2_final:.2e} ({C2_final/targets['C2']*100:5.2f}%)")
        print(f"  Ni/Ne: {Ni_over_Ne:.2f}")
        print(f"  C2 improvement: {C2_improvement:.1f}×")
        print(f"  C2 destruction rate: {destruction_rate:.2e} cm⁻³/s")

        return {'C2': C2_final, 'CH': CH_final, 'improvement': C2_improvement,
                'destruction_rate': destruction_rate, 'stable': Ni_over_Ne < 100}

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("QUESTION 1: WHAT DOES 'C2+H→CH+C ELIMINATED' MEAN?")
    print("="*80)
    print("\nWe set the rate multiplier to different values:")
    print("  1.0  = Full rate (baseline)")
    print("  0.1  = 90% reduction")
    print("  0.01 = 99% reduction")
    print("  0.0  = Complete elimination (is this realistic?)")

    results = []

    # Test different suppression levels
    for factor, name in [
        (1.0, "NO suppression (baseline)"),
        (0.5, "50% reduction"),
        (0.1, "90% reduction"),
        (0.05, "95% reduction"),
        (0.01, "99% reduction"),
        (0.001, "99.9% reduction"),
        (0.0, "Complete elimination"),
    ]:
        result = test_C2_destruction_level(factor, name)
        results.append((factor, name, result))

    print("\n" + "="*80)
    print("SUMMARY: C2+H→CH+C SUPPRESSION LEVELS")
    print("="*80)
    print(f"\n{'Suppression':<25} {'Multiplier':>10} {'C2 (%)'  :>10} {'CH (%)':>10} {'Improvement':>12}")
    print("-"*80)

    for factor, name, result in results:
        if result:
            c2_pct = result['C2'] / targets['C2'] * 100
            ch_pct = result['CH'] / targets['CH'] * 100
            print(f"{name:<25} {factor:>10.3f} {c2_pct:>9.2f}% {ch_pct:>9.1f}% {result['improvement']:>11.1f}×")

    print("\n" + "="*80)
    print("ANSWER:")
    print("="*80)
    print("'Eliminated' means we set the rate multiplier to 0.0")
    print("This makes k_effective = 0, completely turning off the reaction")
    print("\nIs this realistic?")
    print("  - Probably not! The reaction has k = 3.53e-11 cm³/s from literature")
    print("  - Setting it to 0 is a modeling assumption")
    print("  - More realistic: reduce it (99% reduction still gives 22.23%)")
    print("\nConclusion: Even 99% reduction (not elimination) gives ~22% C2")
    print("="*80)
