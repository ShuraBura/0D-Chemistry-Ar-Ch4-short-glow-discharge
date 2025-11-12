"""
Test aggressive C2 optimization WITH three-body electron-ion recombination

Compare to previous best result (28.6× improvement) to see if new physics:
1. Maintains stability
2. Allows even more aggressive multipliers
3. Improves C2 beyond 4.80%
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

# Load baseline with tuned rates
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

    # Get default rates (now includes three-body e-ion recombination!)
    k = define_rates(params)

    # Apply baseline's 23 tuned rates
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    # Apply test multipliers
    print("\nRate multipliers applied:")
    for rate_name, mult in rate_multipliers.items():
        if rate_name in k:
            k[rate_name] *= mult
            print(f"  {rate_name}: {k[rate_name]:.3e} (×{mult})")

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

        # Calculate improvements
        C2_baseline = baseline['target_densities']['C2']
        C2_improvement = C2_final / C2_baseline

        print(f"\nResults:")
        print(f"  H:      {H_final:.2e} ({H_final/targets['H']*100:6.1f}%)  [baseline: 80.3%]")
        print(f"  CH:     {CH_final:.2e} ({CH_final/targets['CH']*100:6.1f}%)  [baseline: 113.0%]")
        print(f"  C2:     {C2_final:.2e} ({C2_final/targets['C2']*100:6.2f}%)  [baseline: 0.17%]")
        print(f"  C2H2:   {C2H2_final:.2e} ({C2H2_final/baseline['all_densities']['C2H2']:.2f}× baseline)")
        print(f"  CH3:    {CH3_final:.2e} ({CH3_final/baseline['all_densities']['CH3']:.2f}× baseline)")
        print(f"  Ni/Ne:  {Ni_over_Ne:.2f}  [baseline: 3.12]")
        print(f"\n  C2 improvement: {C2_improvement:.1f}× from baseline (0.17%)")

        # Check stability
        is_stable = (H_final < 3 * targets['H'] and
                     Ni_over_Ne < 100 and
                     CH3_final < 1e15)

        if not is_stable:
            print(f"  ⚠️  RUNAWAY chemistry!")
        else:
            print(f"  ✓ STABLE chemistry")

        return {
            'H': H_final, 'CH': CH_final, 'C2': C2_final,
            'C2H2': C2H2_final, 'CH3': CH3_final,
            'Ni_over_Ne': Ni_over_Ne,
            'stable': is_stable,
            'C2_improvement': C2_improvement,
        }

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("TEST AGGRESSIVE C2 OPTIMIZATION WITH THREE-BODY E-ION RECOMBINATION")
    print("="*80)
    print(f"\nBaseline (500 mTorr, tuned rates WITHOUT three-body e-ion):")
    print(f"  H:  {baseline['target_densities']['H']:.2e} (80.3%)")
    print(f"  CH: {baseline['target_densities']['CH']:.2e} (113.0%)")
    print(f"  C2: {baseline['target_densities']['C2']:.2e} (0.17%)")
    print(f"  C2H2: {baseline['all_densities']['C2H2']:.2e}")
    print(f"  CH3: {baseline['all_densities']['CH3']:.2e}")
    print(f"  Ni/Ne: {baseline['Ni_over_Ne']:.2f}")

    print(f"\nPrevious best (WITHOUT three-body e-ion recombination):")
    print(f"  20× CH3 + 90% loss reduction → 28.6× C2 improvement (0.17% → 4.80%)")

    # Test 1: Reproduce previous best result WITH new physics
    print("\n" + "="*80)
    print("TEST 1: REPRODUCE PREVIOUS BEST (20× CH3 + 90% LOSS) WITH NEW PHYSICS")
    print("="*80)
    result1 = test_optimization(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 20.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 20.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 20.0,
            'stick_CH3_9_2': 0.1,
            'stick_C2H2_9_11': 0.1,
            'loss_C2H2_11_19': 0.1,
        },
        name="20× CH3 + 90% loss (WITH three-body e-ion)"
    )

    # Test 2: Even more aggressive - 30× CH3
    print("\n" + "="*80)
    print("TEST 2: MORE AGGRESSIVE (30× CH3 + 90% LOSS)")
    print("="*80)
    result2 = test_optimization(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 30.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 30.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 30.0,
            'stick_CH3_9_2': 0.1,
            'stick_C2H2_9_11': 0.1,
            'loss_C2H2_11_19': 0.1,
        },
        name="30× CH3 + 90% loss"
    )

    # Test 3: Extreme - 50× CH3
    print("\n" + "="*80)
    print("TEST 3: EXTREME (50× CH3 + 95% LOSS)")
    print("="*80)
    result3 = test_optimization(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 50.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 50.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 50.0,
            'stick_CH3_9_2': 0.05,
            'stick_C2H2_9_11': 0.05,
            'loss_C2H2_11_19': 0.05,
        },
        name="50× CH3 + 95% loss"
    )

    # Test 4: Also boost C2H2 → C2 reaction
    print("\n" + "="*80)
    print("TEST 4: 30× CH3 + 10× C2 PRODUCTION + 90% LOSS")
    print("="*80)
    result4 = test_optimization(
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 30.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 30.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 30.0,
            'H_C2H2_C2_H2_cm3_7_22': 10.0,  # H + C2H2 → C2 + H2
            'stick_CH3_9_2': 0.1,
            'stick_C2H2_9_11': 0.1,
            'loss_C2H2_11_19': 0.1,
        },
        name="30× CH3 + 10× C2 production + 90% loss"
    )

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    results = [
        ("20× CH3 (previous best)", result1),
        ("30× CH3", result2),
        ("50× CH3", result3),
        ("30× CH3 + 10× C2 production", result4),
    ]

    print("\n{:<35} {:>10} {:>12} {:>10}".format("Configuration", "C2 (%)", "C2 improv.", "Stable?"))
    print("-"*80)
    for name, result in results:
        if result:
            c2_pct = result['C2'] / targets['C2'] * 100
            stable = "✓" if result['stable'] else "✗"
            print(f"{name:<35} {c2_pct:>9.2f}% {result['C2_improvement']:>11.1f}× {stable:>10}")
        else:
            print(f"{name:<35} {'FAILED':>9} {'--':>12} {'✗':>10}")

    print("\n" + "="*80)
    print("Check if three-body e-ion recombination allows:")
    print("  1. More aggressive multipliers without runaway")
    print("  2. Better C2 improvement than previous 28.6×")
    print("="*80)
