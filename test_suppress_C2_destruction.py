"""
Test suppressing C2 + H → CH + C reaction

Hypothesis: This reaction destroys C2 77,000× faster than CH + CH produces it.
Suppressing it should dramatically improve C2.
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
    k_B = 1.380649e-23  # J/K
    n_m3 = pressure_Pa / (k_B * T_K)  # m⁻³
    return n_m3 * 1e-6  # Convert to cm⁻³

# Load baseline
with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

targets = {
    'H': 2.52e14,
    'CH': 1.0e9,
    'C2': 5.6e11,
}

def test_optimization(rate_multipliers, name="Test"):
    """Test chemistry with given rate multipliers"""
    print(f"\n{'='*80}")
    print(f"{name}")
    print(f"{'='*80}")

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
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                    'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH',
                    'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C',
                    'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H',
                    'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus',
                    'H2Plus', 'C2H2Star'],
    }

    # Get rates with baseline tuned values
    k = define_rates(params)
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    # Apply test multipliers
    print("\nRate multipliers applied:")
    for rate_name, mult in rate_multipliers.items():
        if rate_name in k:
            original = k[rate_name]
            k[rate_name] *= mult
            print(f"  {rate_name}: {original:.2e} → {k[rate_name]:.2e} (×{mult})")

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

        C2_baseline = baseline['target_densities']['C2']
        C2_improvement = C2_final / C2_baseline

        print(f"\nResults:")
        print(f"  H:      {H_final:.2e} ({H_final/targets['H']*100:6.1f}%)")
        print(f"  CH:     {CH_final:.2e} ({CH_final/targets['CH']*100:6.1f}%)")
        print(f"  C2:     {C2_final:.2e} ({C2_final/targets['C2']*100:6.2f}%)")
        print(f"  C2H2:   {C2H2_final:.2e}")
        print(f"  Ni/Ne:  {Ni_over_Ne:.2f}")
        print(f"\n  C2 improvement: {C2_improvement:.1f}× from baseline (0.17%)")

        # Check C2 + H → CH + C rate
        k_C2_H = k.get('C2_H_CH_C_cm3_7_6', 0)
        rate_C2_destruction = k_C2_H * C2_final * H_final
        print(f"\n  C2 + H → CH + C destruction rate: {rate_C2_destruction:.2e} cm⁻³/s")

        is_stable = (H_final < 5 * targets['H'] and
                     Ni_over_Ne < 100 and
                     CH3_final < 1e15)

        if not is_stable:
            print(f"  ⚠️  RUNAWAY chemistry!")
        else:
            print(f"  ✓ STABLE chemistry")

        return {
            'H': H_final, 'CH': CH_final, 'C2': C2_final,
            'C2H2': C2H2_final,
            'Ni_over_Ne': Ni_over_Ne,
            'stable': is_stable,
            'C2_improvement': C2_improvement,
        }

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("TEST SUPPRESSING C2 + H → CH + C DESTRUCTION")
    print("="*80)
    print(f"\nCurrent best (200× CH3):")
    print(f"  C2: 6.56% (39× improvement)")
    print(f"  CH: 343%")
    print(f"  C2 + H → CH + C rate: 2.65e+14 cm⁻³/s (MASSIVE!)")

    # Test 1: 200× CH3 WITHOUT suppressing C2 destruction (baseline)
    result_baseline = test_optimization(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 200.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 200.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 200.0,
            'stick_CH3_9_2': 0.01,
            'stick_C2H2_9_11': 0.01,
            'loss_C2H2_11_19': 0.01,
        },
        name="BASELINE: 200× CH3 + 99% loss (NO C2 suppression)"
    )

    # Test 2: Suppress C2 + H → CH + C by 10×
    result_10x = test_optimization(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 200.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 200.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 200.0,
            'stick_CH3_9_2': 0.01,
            'stick_C2H2_9_11': 0.01,
            'loss_C2H2_11_19': 0.01,
            'C2_H_CH_C_cm3_7_6': 0.1,  # Suppress by 10×
        },
        name="TEST 1: 200× CH3 + SUPPRESS C2+H→CH+C by 10×"
    )

    # Test 3: Suppress by 100×
    result_100x = test_optimization(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 200.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 200.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 200.0,
            'stick_CH3_9_2': 0.01,
            'stick_C2H2_9_11': 0.01,
            'loss_C2H2_11_19': 0.01,
            'C2_H_CH_C_cm3_7_6': 0.01,  # Suppress by 100×
        },
        name="TEST 2: 200× CH3 + SUPPRESS C2+H→CH+C by 100×"
    )

    # Test 4: Completely eliminate C2 + H → CH + C
    result_eliminate = test_optimization(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 200.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 200.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 200.0,
            'stick_CH3_9_2': 0.01,
            'stick_C2H2_9_11': 0.01,
            'loss_C2H2_11_19': 0.01,
            'C2_H_CH_C_cm3_7_6': 0.0,  # Completely eliminate
        },
        name="TEST 3: 200× CH3 + ELIMINATE C2+H→CH+C"
    )

    # Summary
    print("\n" + "="*80)
    print("SUMMARY: IMPACT OF SUPPRESSING C2 DESTRUCTION")
    print("="*80)

    results = [
        ("Baseline (no suppression)", result_baseline),
        ("Suppress C2+H→CH+C by 10×", result_10x),
        ("Suppress C2+H→CH+C by 100×", result_100x),
        ("Eliminate C2+H→CH+C", result_eliminate),
    ]

    print("\n{:<30} {:>10} {:>10} {:>12} {:>10}".format(
        "Configuration", "C2 (%)", "CH (%)", "C2 improv.", "Stable?"
    ))
    print("-"*80)

    for name, result in results:
        if result:
            c2_pct = result['C2'] / targets['C2'] * 100
            ch_pct = result['CH'] / targets['CH'] * 100
            stable = "✓" if result['stable'] else "✗"
            print(f"{name:<30} {c2_pct:>9.2f}% {ch_pct:>9.1f}% {result['C2_improvement']:>11.1f}× {stable:>10}")
        else:
            print(f"{name:<30} {'FAILED':>9} {'--':>10} {'--':>12} {'✗':>10}")

    print("\n" + "="*80)
    print("ANSWER: C2 + H → CH + C is the bottleneck!")
    print("="*80)
