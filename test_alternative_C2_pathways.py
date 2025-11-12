"""
Investigate alternative C2 production pathways that bypass C2H2 bottleneck

Since C2H2 is saturated, test pathways like:
1. Boost C production and use C + C + M → C2
2. Boost C2H species (C2H, C2H3) production
3. Boost reactions that convert C2H → C2
4. Test if we can boost atomic carbon and form C2 directly
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
    if rate_multipliers:
        print("\nRate multipliers applied:")
        for rate_name, mult in rate_multipliers.items():
            if rate_name in k:
                original = k[rate_name]
                k[rate_name] *= mult
                print(f"  {rate_name}: ×{mult} ({original:.2e} → {k[rate_name]:.2e})")

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
        C_final = y_final[species.index('C')]
        C2H_final = y_final[species.index('C2H')]

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
        print(f"  C:      {C_final:.2e}")
        print(f"  C2H:    {C2H_final:.2e}")
        print(f"  C2H2:   {C2H2_final:.2e}")
        print(f"  Ni/Ne:  {Ni_over_Ne:.2f}")
        print(f"\n  C2 improvement: {C2_improvement:.1f}× from baseline (0.17%)")

        # Calculate C + C + M → C2 rate
        k_C_C_M = k.get('C_C_M_C2_M_cm6_7_64', 0)
        rate_C_C = k_C_C_M * C_final * C_final
        print(f"\n  C + C + M → C2 rate: {rate_C_C:.2e} cm⁻³/s")

        # Calculate C2H + H → C2 rate if exists
        k_C2H_H = k.get('C2H_H_C2_H2_cm3_7_47', 0)
        rate_C2H = k_C2H_H * C2H_final * H_final
        print(f"  C2H + H → C2 rate: {rate_C2H:.2e} cm⁻³/s")

        is_stable = (H_final < 5 * targets['H'] and
                     Ni_over_Ne < 100)

        if not is_stable:
            print(f"  ⚠️  RUNAWAY chemistry!")
        else:
            print(f"  ✓ STABLE chemistry")

        return {
            'H': H_final, 'CH': CH_final, 'C2': C2_final,
            'C': C_final, 'C2H': C2H_final, 'C2H2': C2H2_final,
            'Ni_over_Ne': Ni_over_Ne,
            'stable': is_stable,
            'C2_improvement': C2_improvement,
        }

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("TEST ALTERNATIVE C2 PATHWAYS (BYPASS C2H2 BOTTLENECK)")
    print("="*80)

    # Baseline: best so far
    base_multipliers = {
        'e_CH4_CH3_HMinus_cm3_8_1': 200.0,
        'ArStar_CH4_CH3_H_cm3_3_1': 200.0,
        'e_CH4Plus_CH3_H_cm3_6_4': 200.0,
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
        'C2_H_CH_C_cm3_7_6': 0.0,
    }

    result_base = test_optimization(base_multipliers, "BASELINE: Current best (22.78%)")

    # Strategy 1: Boost atomic C production from CH
    print("\n" + "="*80)
    print("STRATEGY 1: BOOST ATOMIC CARBON PRODUCTION")
    print("="*80)

    result_1 = test_optimization(
        {**base_multipliers,
         'e_CH_C_H_e_cm3_1_12': 100.0,  # e + CH → C + H (electron impact)
         'C_C_M_C2_M_cm6_7_64': 1000.0},  # C + C + M → C2 (boost three-body)
        "1: Boost C production from CH + boost C+C→C2"
    )

    # Strategy 2: Boost C2H production and conversion to C2
    print("\n" + "="*80)
    print("STRATEGY 2: BOOST C2H PRODUCTION AND C2H → C2")
    print("="*80)

    result_2 = test_optimization(
        {**base_multipliers,
         'C2H_H_C2_H2_cm3_7_47': 1000.0},  # C2H + H → C2 + H2
        "2: Boost C2H + H → C2"
    )

    # Strategy 3: Boost CH + C → C2 pathways
    print("\n" + "="*80)
    print("STRATEGY 3: BOOST CH + C → C2 PATHWAYS")
    print("="*80)

    result_3a = test_optimization(
        {**base_multipliers,
         'CH_C_C2_H_cm3_7_9': 100.0,   # CH + C → C2 + H
         'C_CH_C2_H_cm3_7_4': 100.0},   # C + CH → C2 + H
        "3a: Boost CH + C → C2 reactions"
    )

    result_3b = test_optimization(
        {**base_multipliers,
         'CH_C_C2_H2_cm3_7_24': 100.0},  # CH + C → C2 + H2
        "3b: Boost CH + C → C2 + H2"
    )

    # Strategy 4: Combined atomic carbon pathway
    print("\n" + "="*80)
    print("STRATEGY 4: COMBINED - ATOMIC CARBON PATHWAY")
    print("="*80)

    result_4 = test_optimization(
        {**base_multipliers,
         'e_CH_C_H_e_cm3_1_12': 100.0,      # Boost C production
         'C_C_M_C2_M_cm6_7_64': 10000.0,    # Boost C + C → C2
         'CH_C_C2_H_cm3_7_9': 100.0,        # CH + C → C2
         'C_CH_C2_H_cm3_7_4': 100.0,        # C + CH → C2
         'CH_C_C2_H2_cm3_7_24': 100.0},     # CH + C → C2 + H2
        "4: COMBINED atomic carbon pathway"
    )

    # Strategy 5: Suppress C consumption reactions
    print("\n" + "="*80)
    print("STRATEGY 5: SUPPRESS C CONSUMPTION (BUILD UP C)")
    print("="*80)

    result_5 = test_optimization(
        {**base_multipliers,
         'CH2_C_C2H2_cm3_7_17': 0.1,        # Reduce CH2 + C → C2H2
         'C_C2H2_C2_CH2_cm3_7_16': 0.1,     # Reduce C + C2H2 → C2 + CH2
         'C_C_M_C2_M_cm6_7_64': 10000.0},   # Boost C + C → C2
        "5: Suppress C losses + boost C+C→C2"
    )

    # Summary
    print("\n" + "="*80)
    print("SUMMARY: ALTERNATIVE C2 PATHWAYS")
    print("="*80)

    results = [
        ("Baseline (22.78%)", result_base),
        ("1: Boost C prod + C+C→C2", result_1),
        ("2: Boost C2H→C2", result_2),
        ("3a: Boost CH+C→C2", result_3a),
        ("3b: Boost CH+C→C2+H2", result_3b),
        ("4: Combined C pathway", result_4),
        ("5: Suppress C loss", result_5),
    ]

    print("\n{:<30} {:>10} {:>10} {:>10} {:>12} {:>10}".format(
        "Strategy", "C2 (%)", "C (cm⁻³)", "C2H", "C2 improv.", "Stable?"
    ))
    print("-"*90)

    for name, result in results:
        if result:
            c2_pct = result['C2'] / targets['C2'] * 100
            stable = "✓" if result['stable'] else "✗"
            print(f"{name:<30} {c2_pct:>9.2f}% {result['C']:.2e} {result['C2H']:.2e} {result['C2_improvement']:>11.1f}× {stable:>10}")
        else:
            print(f"{name:<30} {'FAILED':>9} {'--':>11} {'--':>10} {'--':>12} {'✗':>10}")

    print("\n" + "="*80)
    print("Can we bypass the C2H2 bottleneck using atomic carbon?")
    print("="*80)
