"""
Test impact of:
1. New CH3 production pathways (10 added)
2. Increased ne from 1.22e8 to 2.3e9 cm⁻³ (20× higher!)

User: "Ne should be ~2.3e9"
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

def test_configuration(ne_value, use_new_pathways, Te_mult, name):
    """Test with different ne and pathway configuration"""
    print(f"\n{'='*80}")
    print(f"{name}")
    print(f"{'='*80}")
    print(f"ne: {ne_value:.2e} cm⁻³")
    print(f"Te: {baseline['Te'] * Te_mult:.3f} eV (×{Te_mult})")
    print(f"New CH3 pathways: {'YES' if use_new_pathways else 'NO'}")

    P = 500.0
    n_total = pressure_to_density(P)

    params = {
        'P': P,
        'Te': baseline['Te'] * Te_mult,
        'ne': ne_value,  # User-specified ne
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

    # Apply optimization: modest CH3 boost + loss reduction + correct chemistry
    rate_mults = {
        'e_CH4_CH3_HMinus_cm3_8_1': 10.0,  # Modest boost (within literature range)
        'ArStar_CH4_CH3_H_cm3_3_1': 1.4,   # Max physical (at literature limit)
        'e_CH4Plus_CH3_H_cm3_6_4': 1.5,    # Near literature limit
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
        # C2_H_CH_C_cm3_7_6 kept at literature value
    }

    # If not using new pathways, suppress them
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
                k[key] = 0.0  # Turn off

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
            print(f"  ✗ Integration failed!")
            return None

        y_final = sol.y[:, -1]
        H_final = y_final[species.index('H')]
        CH_final = y_final[species.index('CH')]
        C2_final = y_final[species.index('C2')]
        CH3_final = y_final[species.index('CH3')]

        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus', 'CHPlus', 'H2Plus']
        n_i_total = sum(y_final[species.index(ion)] for ion in ions)
        ne_final = y_final[species.index('e')]
        Ni_over_Ne = n_i_total / ne_final if ne_final > 0 else 0

        C2_improvement = C2_final / baseline['target_densities']['C2']

        print(f"\nResults:")
        print(f"  H:    {H_final:.2e} ({H_final/targets['H']*100:5.1f}%)")
        print(f"  CH:   {CH_final:.2e} ({CH_final/targets['CH']*100:5.1f}%)")
        print(f"  C2:   {C2_final:.2e} ({C2_final/targets['C2']*100:5.2f}%)")
        print(f"  CH3:  {CH3_final:.2e}")
        print(f"  Ni/Ne: {Ni_over_Ne:.2f}")
        print(f"  C2 improvement: {C2_improvement:.1f}×")

        is_stable = (Ni_over_Ne < 10 and H_final < 5 * targets['H'])

        if not is_stable:
            print(f"  ⚠️  RUNAWAY!")
        else:
            print(f"  ✓ STABLE")

        return {
            'C2': C2_final, 'CH': CH_final, 'H': H_final, 'CH3': CH3_final,
            'Ni_over_Ne': Ni_over_Ne, 'improvement': C2_improvement, 'stable': is_stable
        }

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("TEST NEW CH3 PATHWAYS + INCREASED ne")
    print("="*80)
    print("\nUser correction: ne should be ~2.3e9 cm⁻³ (currently 1.22e8)")
    print("Added 10 new CH3 production pathways")
    print("Using physically realistic rate multipliers:")
    print("  - e + CH4: 10× (within 1e-13 to 1e-8 range)")
    print("  - Ar* + CH4: 1.4× (at literature limit)")
    print("  - e + CH4+: 1.5× (near literature limit)")

    # Current baseline ne
    ne_baseline = 1.22e8
    # User-specified ne
    ne_user = 2.3e9

    results = []

    # Test 1: Baseline (old ne, no new pathways)
    r1 = test_configuration(ne_baseline, False, 1.0, "BASELINE: Old ne, no new pathways")
    results.append(("Baseline", r1))

    # Test 2: Add new pathways (old ne)
    r2 = test_configuration(ne_baseline, True, 1.0, "TEST 1: Old ne, WITH new pathways")
    results.append(("+ New pathways", r2))

    # Test 3: Increase ne (no new pathways)
    r3 = test_configuration(ne_user, False, 1.0, "TEST 2: New ne (2.3e9), no new pathways")
    results.append(("+ Higher ne", r3))

    # Test 4: New pathways + new ne
    r4 = test_configuration(ne_user, True, 1.0, "TEST 3: New ne + new pathways")
    results.append(("+ Both", r4))

    # Test 5: New pathways + new ne + optimal Te
    r5 = test_configuration(ne_user, True, 1.1, "TEST 4: New ne + new pathways + Te×1.1")
    results.append(("+ Optimal Te", r5))

    # Summary
    print("\n" + "="*80)
    print("SUMMARY: IMPACT OF NEW PATHWAYS AND INCREASED ne")
    print("="*80)

    print(f"\n{'Configuration':<25} {'C2 (%)':<12} {'CH3 (cm⁻³)':<15} {'Improvement':<12} {'Stable?':<10}")
    print("-"*80)

    for name, result in results:
        if result and result['stable']:
            c2_pct = result['C2'] / targets['C2'] * 100
            print(f"{name:<25} {c2_pct:>10.2f}% {result['CH3']:>13.2e} {result['improvement']:>10.1f}× {'✓':<10}")
        elif result:
            c2_pct = result['C2'] / targets['C2'] * 100
            print(f"{name:<25} {c2_pct:>10.2f}% {result['CH3']:>13.2e} {result['improvement']:>10.1f}× {'✗ RUNAWAY':<10}")
        else:
            print(f"{name:<25} {'FAILED':<12} {'--':<15} {'--':<12} {'✗':<10}")

    print("\n" + "="*80)
    print("Can we achieve better C2 with new pathways + correct ne?")
    print("="*80)
