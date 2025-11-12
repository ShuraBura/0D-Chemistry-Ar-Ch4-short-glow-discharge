"""
Test boosting C2 production pathways to push beyond 22.78%

Strategies:
1. Boost C2H2 + H → C2 + H2 + H (currently 1e-11, very small!)
2. Boost CH + CH → C2 + H2 (currently producing 3.45e9 cm⁻³/s)
3. Boost C + C + M → C2 (three-body, currently negligible)
4. Boost electron impact: e + C2H2 → C2 + H2
5. Combine with eliminating C2 + H → CH + C
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
        C_final = y_final[species.index('C')]

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
        print(f"  C:      {C_final:.2e}")
        print(f"  Ni/Ne:  {Ni_over_Ne:.2f}")
        print(f"\n  C2 improvement: {C2_improvement:.1f}× from baseline (0.17%)")

        is_stable = (H_final < 5 * targets['H'] and
                     Ni_over_Ne < 100 and
                     CH3_final < 1e15)

        if not is_stable:
            print(f"  ⚠️  RUNAWAY chemistry!")
        else:
            print(f"  ✓ STABLE chemistry")

        return {
            'H': H_final, 'CH': CH_final, 'C2': C2_final,
            'C2H2': C2H2_final, 'C': C_final,
            'Ni_over_Ne': Ni_over_Ne,
            'stable': is_stable,
            'C2_improvement': C2_improvement,
        }

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("BOOST C2 PRODUCTION PATHWAYS")
    print("="*80)
    print(f"\nCurrent best (200× CH3 + eliminate C2+H→CH+C):")
    print(f"  C2: 22.78% (135.5× improvement)")

    # Baseline: 200× CH3 + eliminate C2 destruction
    base_multipliers = {
        'e_CH4_CH3_HMinus_cm3_8_1': 200.0,
        'ArStar_CH4_CH3_H_cm3_3_1': 200.0,
        'e_CH4Plus_CH3_H_cm3_6_4': 200.0,
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
        'C2_H_CH_C_cm3_7_6': 0.0,  # Eliminate C2 destruction
    }

    result_base = test_optimization(base_multipliers, "BASELINE: Current best")

    # Strategy 1: Boost C2H2 + H → C2 + H2 + H
    print("\n" + "="*80)
    print("STRATEGY 1: BOOST C2H2 + H → C2 + H2 + H")
    print("="*80)
    print("Current k = 1e-11 (VERY SMALL!)")

    result_1a = test_optimization(
        {**base_multipliers, 'C2H2_H_C2_H2_H_cm3_7_50': 10.0},
        "1a: Boost C2H2+H→C2 by 10×"
    )

    result_1b = test_optimization(
        {**base_multipliers, 'C2H2_H_C2_H2_H_cm3_7_50': 100.0},
        "1b: Boost C2H2+H→C2 by 100×"
    )

    result_1c = test_optimization(
        {**base_multipliers, 'C2H2_H_C2_H2_H_cm3_7_50': 1000.0},
        "1c: Boost C2H2+H→C2 by 1000×"
    )

    # Strategy 2: Boost CH + CH → C2
    print("\n" + "="*80)
    print("STRATEGY 2: BOOST CH + CH → C2 + H2")
    print("="*80)

    result_2a = test_optimization(
        {**base_multipliers,
         'CH_CH_C2_H2_cm3_5_4': 10.0,
         'CH_CH_C2_H2_cm3_7_44': 10.0},
        "2a: Boost CH+CH→C2 by 10×"
    )

    result_2b = test_optimization(
        {**base_multipliers,
         'CH_CH_C2_H2_cm3_5_4': 100.0,
         'CH_CH_C2_H2_cm3_7_44': 100.0},
        "2b: Boost CH+CH→C2 by 100×"
    )

    # Strategy 3: Boost electron impact e + C2H2 → C2
    print("\n" + "="*80)
    print("STRATEGY 3: BOOST e + C2H2 → C2 + H2")
    print("="*80)

    result_3 = test_optimization(
        {**base_multipliers, 'e_C2H2_C2_H2_cm3_1_16': 100.0},
        "3: Boost e+C2H2→C2 by 100×"
    )

    # Strategy 4: Boost C2H2 + C → C2 + CH2
    print("\n" + "="*80)
    print("STRATEGY 4: BOOST C2H2 + C → C2 + CH2")
    print("="*80)

    result_4 = test_optimization(
        {**base_multipliers, 'C2H2_C_C2_CH2_cm3_7_19': 100.0},
        "4: Boost C2H2+C→C2 by 100×"
    )

    # Strategy 5: COMBINED APPROACH
    print("\n" + "="*80)
    print("STRATEGY 5: COMBINED - BOOST ALL C2 PRODUCTION")
    print("="*80)

    result_5 = test_optimization(
        {**base_multipliers,
         'C2H2_H_C2_H2_H_cm3_7_50': 1000.0,  # C2H2 + H → C2
         'CH_CH_C2_H2_cm3_5_4': 100.0,        # CH + CH → C2
         'CH_CH_C2_H2_cm3_7_44': 100.0,       # CH + CH → C2 (alt)
         'e_C2H2_C2_H2_cm3_1_16': 100.0,      # e + C2H2 → C2
         'C2H2_C_C2_CH2_cm3_7_19': 100.0},    # C2H2 + C → C2
        "5: COMBINED - All C2 production boosted"
    )

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    results = [
        ("Baseline (current best)", result_base),
        ("1a: C2H2+H→C2 ×10", result_1a),
        ("1b: C2H2+H→C2 ×100", result_1b),
        ("1c: C2H2+H→C2 ×1000", result_1c),
        ("2a: CH+CH→C2 ×10", result_2a),
        ("2b: CH+CH→C2 ×100", result_2b),
        ("3: e+C2H2→C2 ×100", result_3),
        ("4: C2H2+C→C2 ×100", result_4),
        ("5: COMBINED boost", result_5),
    ]

    print("\n{:<30} {:>10} {:>10} {:>12} {:>10}".format(
        "Strategy", "C2 (%)", "CH (%)", "C2 improv.", "Stable?"
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
    print("Finding the best strategy to push C2 beyond 22.78%...")
    print("="*80)
