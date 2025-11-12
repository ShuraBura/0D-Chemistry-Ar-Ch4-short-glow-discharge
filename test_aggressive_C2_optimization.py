"""
AGGRESSIVE C2 OPTIMIZATION AT 500 mTorr

Test two strategies:
1. Aggressive CH3 production boost (5×, 10×, 20×) + reduced losses
2. Boost H + C2H2 → C2 rate constant (2×, 5×, 10×)
3. Combination of both

Goal: Push C2 from current 1.3% toward 10%+ while maintaining stability
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

def pressure_to_density(pressure_mTorr, T_K=300):
    """Convert pressure to total density (cm⁻³)"""
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

# Load baseline with tuned rates
with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

targets = {
    'H': 2.52e14,
    'CH': 1.0e9,
    'C2': 5.6e11,
}

def test_scenario(rate_multipliers, name):
    """Test a specific rate multiplier scenario at 500 mTorr"""
    print(f"\n{'='*80}")
    print(f"TEST: {name}")
    print(f"{'='*80}")

    pressure_mTorr = 500.0
    n_total = pressure_to_density(pressure_mTorr)

    # Use baseline plasma parameters
    ne_frac = baseline['Ne'] / pressure_to_density(500.0)
    ne = ne_frac * n_total
    Te = baseline['Te']
    E_field = baseline['E_field']

    params = {
        'P': pressure_mTorr,
        'Te': Te,
        'ne': ne,
        'E_field': E_field,
        'L_discharge': 0.45,
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
        },
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                    'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH',
                    'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C',
                    'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H',
                    'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus',
                    'H2Plus', 'C2H2Star'],
    }

    k = define_rates(params)

    # Apply baseline's 23 tuned rates first
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    # Apply test multipliers
    print("Rate multipliers applied:")
    for rate_name, mult in rate_multipliers.items():
        if rate_name in k:
            k[rate_name] *= mult
            print(f"  {rate_name:50s}: ×{mult:4.1f} → {k[rate_name]:.3e}")

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    species = params['species']
    y0 = np.ones(len(species)) * 1e3
    y0[species.index('Ar')] = n_total * 0.85
    y0[species.index('CH4')] = n_total * 0.15
    y0[species.index('e')] = ne

    try:
        ode_func = PlasmaODE_Optimized(params)
        sol = solve_ivp(
            ode_func, (0, 500), y0,
            method='BDF', rtol=1e-7, atol=1e-9, max_step=1.0
        )

        if not sol.success:
            print(f"  ✗ Integration failed!")
            return None

        y_final = sol.y[:, -1]

        H_final = y_final[species.index('H')]
        CH_final = y_final[species.index('CH')]
        C2_final = y_final[species.index('C2')]
        C2H2_final = y_final[species.index('C2H2')]
        CH3_final = y_final[species.index('CH3')]

        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus', 'CHPlus', 'H2Plus']
        n_i_total = sum(y_final[species.index(ion)] for ion in ions)
        ne_final = y_final[species.index('e')]
        Ni_over_Ne = n_i_total / ne_final if ne_final > 0 else 0

        # Calculate improvements vs baseline
        C2_improvement = C2_final / baseline['target_densities']['C2']
        C2H2_boost = C2H2_final / baseline['all_densities']['C2H2']
        CH3_boost = CH3_final / baseline['all_densities']['CH3']

        print(f"\nResults:")
        print(f"  H:      {H_final:.2e} ({H_final/targets['H']*100:6.1f}%)")
        print(f"  CH:     {CH_final:.2e} ({CH_final/targets['CH']*100:6.1f}%)")
        print(f"  C2:     {C2_final:.2e} ({C2_final/targets['C2']*100:6.1f}%) [×{C2_improvement:.1f} vs baseline]")
        print(f"  C2H2:   {C2H2_final:.2e} [×{C2H2_boost:.1f} baseline]")
        print(f"  CH3:    {CH3_final:.2e} [×{CH3_boost:.1f} baseline]")
        print(f"  Ni/Ne:  {Ni_over_Ne:.2f}")

        # Check stability
        is_stable = (H_final < 5 * targets['H'] and
                     Ni_over_Ne < 100 and
                     CH3_final < 1e15)

        if not is_stable:
            print(f"  ⚠️  RUNAWAY chemistry!")
            return None
        else:
            status = "✓ STABLE"
            if C2_improvement > 1.0:
                status += f" + C2 IMPROVED {C2_improvement:.1f}×!"
            print(f"  {status}")

        return {
            'name': name,
            'H': H_final, 'CH': CH_final, 'C2': C2_final,
            'C2H2': C2H2_final, 'CH3': CH3_final,
            'Ni_over_Ne': Ni_over_Ne,
            'C2_pct': C2_final/targets['C2']*100,
            'C2_improvement': C2_improvement,
            'stable': is_stable,
        }

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("AGGRESSIVE C2 OPTIMIZATION AT 500 mTorr")
    print("="*80)
    print(f"\nCurrent best (3× CH3 + 50% losses):")
    print(f"  C2: 1.3% (7.43e9 cm⁻³)")
    print(f"\nGoal: Push toward 10% while maintaining stability")
    print(f"Target: 5.6e11 cm⁻³ (100%)")

    results = []

    # ========================================================================
    # STRATEGY 1: AGGRESSIVE CH3 PRODUCTION BOOST
    # ========================================================================
    print("\n" + "="*80)
    print("STRATEGY 1: AGGRESSIVE CH3 PRODUCTION BOOST")
    print("="*80)

    # Test 1.1: 5× CH3, 70% loss reduction
    results.append(test_scenario(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 5.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 5.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 5.0,
            'stick_CH3_9_2': 0.3,
            'stick_C2H2_9_11': 0.3,
            'loss_C2H2_11_19': 0.3,
        },
        "5× CH3 + 70% loss reduction"
    ))

    # Test 1.2: 10× CH3, 80% loss reduction
    results.append(test_scenario(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 10.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 10.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 10.0,
            'stick_CH3_9_2': 0.2,
            'stick_C2H2_9_11': 0.2,
            'loss_C2H2_11_19': 0.2,
        },
        "10× CH3 + 80% loss reduction"
    ))

    # Test 1.3: 20× CH3, 90% loss reduction (very aggressive!)
    results.append(test_scenario(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 20.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 20.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 20.0,
            'stick_CH3_9_2': 0.1,
            'stick_C2H2_9_11': 0.1,
            'loss_C2H2_11_19': 0.1,
        },
        "20× CH3 + 90% loss reduction"
    ))

    # ========================================================================
    # STRATEGY 2: BOOST C2 PRODUCTION RATE
    # ========================================================================
    print("\n" + "="*80)
    print("STRATEGY 2: BOOST H + C2H2 → C2 RATE CONSTANT")
    print("="*80)
    print("Baseline rate: C2H2_H_C2_H2_H_cm3_7_50 = 9.478e-12 cm³/s")

    # Test 2.1: 2× C2 production rate only
    results.append(test_scenario(
        {
            'C2H2_H_C2_H2_H_cm3_7_50': 2.0,
        },
        "2× C2 production rate (baseline chemistry)"
    ))

    # Test 2.2: 5× C2 production rate only
    results.append(test_scenario(
        {
            'C2H2_H_C2_H2_H_cm3_7_50': 5.0,
        },
        "5× C2 production rate (baseline chemistry)"
    ))

    # Test 2.3: 10× C2 production rate only
    results.append(test_scenario(
        {
            'C2H2_H_C2_H2_H_cm3_7_50': 10.0,
        },
        "10× C2 production rate (baseline chemistry)"
    ))

    # ========================================================================
    # STRATEGY 3: COMBINED APPROACH
    # ========================================================================
    print("\n" + "="*80)
    print("STRATEGY 3: COMBINED APPROACH (CH3 BOOST + C2 RATE BOOST)")
    print("="*80)

    # Test 3.1: Moderate both (3× CH3 + 5× C2 rate)
    results.append(test_scenario(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 3.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 3.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 3.0,
            'stick_CH3_9_2': 0.5,
            'stick_C2H2_9_11': 0.5,
            'loss_C2H2_11_19': 0.5,
            'C2H2_H_C2_H2_H_cm3_7_50': 5.0,
        },
        "3× CH3 + 50% losses + 5× C2 rate"
    ))

    # Test 3.2: Aggressive both (5× CH3 + 10× C2 rate)
    results.append(test_scenario(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 5.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 5.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 5.0,
            'stick_CH3_9_2': 0.3,
            'stick_C2H2_9_11': 0.3,
            'loss_C2H2_11_19': 0.3,
            'C2H2_H_C2_H2_H_cm3_7_50': 10.0,
        },
        "5× CH3 + 70% losses + 10× C2 rate"
    ))

    # Test 3.3: Very aggressive (10× CH3 + 20× C2 rate)
    results.append(test_scenario(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 10.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 10.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 10.0,
            'stick_CH3_9_2': 0.2,
            'stick_C2H2_9_11': 0.2,
            'loss_C2H2_11_19': 0.2,
            'C2H2_H_C2_H2_H_cm3_7_50': 20.0,
        },
        "10× CH3 + 80% losses + 20× C2 rate"
    ))

    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("\n" + "="*80)
    print("SUMMARY OF STABLE RESULTS")
    print("="*80)

    stable_results = [r for r in results if r and r['stable']]

    if stable_results:
        print(f"\n{'Test Name':<50s} {'C2 (%)':<10s} {'C2 (cm⁻³)':<12s} {'Improvement':<12s}")
        print("-" * 84)

        for r in sorted(stable_results, key=lambda x: x['C2_pct'], reverse=True):
            print(f"{r['name']:<50s} {r['C2_pct']:>7.2f}%   {r['C2']:.2e}   {r['C2_improvement']:>8.1f}×")

        best = max(stable_results, key=lambda x: x['C2_pct'])
        print("\n" + "="*80)
        print(f"BEST RESULT: {best['name']}")
        print("="*80)
        print(f"  H:  {best['H']:.2e} ({best['H']/targets['H']*100:.1f}%)")
        print(f"  CH: {best['CH']:.2e} ({best['CH']/targets['CH']*100:.1f}%)")
        print(f"  C2: {best['C2']:.2e} ({best['C2_pct']:.2f}%)")
        print(f"  C2 improvement: {best['C2_improvement']:.1f}× vs baseline")
        print(f"  Ni/Ne: {best['Ni_over_Ne']:.2f}")
    else:
        print("\n⚠️  No stable results found! All tests went into runaway.")

    print("\n" + "="*80)
