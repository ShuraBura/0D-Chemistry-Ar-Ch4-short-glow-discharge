"""
Systematic test: Boost CH3 production to increase C2H2

Starting from baseline (which achieves H=79.9%, CH=100.7%, C2=16.8%),
test specific multipliers on CH3 production rates.

Current: CH3 = 6.85e11, C2H2 = 4.81e9
Target:  CH3 ~ 2e12 (3× increase), C2H2 ~ 5e10 (10× increase)

Since C2H2 ∝ [CH3]², a 3× CH3 increase should give 9× C2H2 increase.
"""

import numpy as np
import json
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

def pressure_to_density(pressure_mTorr, T_K=300):
    """Convert pressure to total density"""
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    return pressure_Pa / (k_B * T_K)

# Load baseline
with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

# Targets
target_densities = {
    'H': 2.52e14,
    'CH': 1.0e9,
    'C2': 5.6e11,
}

def test_scenario(scenario_name, rate_multipliers, Te_mult=1.0, ne_mult=1.0):
    """
    Test a specific scenario with given rate multipliers

    rate_multipliers: dict of {rate_name: multiplier}
    """
    print(f"\n{'='*80}")
    print(f"SCENARIO: {scenario_name}")
    print(f"{'='*80}")

    Te = baseline['Te'] * Te_mult
    ne = baseline['Ne'] * ne_mult
    E_field = baseline['E_field']
    pressure_mTorr = 500.0
    n_total = pressure_to_density(pressure_mTorr)

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

    # Apply multipliers
    for rate_name, mult in rate_multipliers.items():
        if rate_name in k:
            k[rate_name] *= mult
            print(f"  {rate_name}: ×{mult}")

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
            print("  Integration failed!")
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

        print(f"\nResults:")
        print(f"  H:      {H_final:.2e} ({H_final/target_densities['H']*100:6.1f}%) [target: 79.9%]")
        print(f"  CH:     {CH_final:.2e} ({CH_final/target_densities['CH']*100:6.1f}%) [target: 100.7%]")
        print(f"  C2:     {C2_final:.2e} ({C2_final/target_densities['C2']*100:6.1f}%) [target: >16.8%]")
        print(f"  C2H2:   {C2H2_final:.2e} ({C2H2_final/baseline['all_densities']['C2H2']:.1f}× baseline)")
        print(f"  CH3:    {CH3_final:.2e} ({CH3_final/baseline['all_densities']['CH3']:.1f}× baseline)")
        print(f"  Ni/Ne:  {Ni_over_Ne:.2f} [baseline: 3.12]")

        # Assess
        H_ok = 70 <= H_final/target_densities['H']*100 <= 90
        CH_ok = 90 <= CH_final/target_densities['CH']*100 <= 110
        C2_improved = C2_final > baseline['target_densities']['C2']
        charge_ok = 2 <= Ni_over_Ne <= 7

        status = []
        if H_ok: status.append("H✓")
        if CH_ok: status.append("CH✓")
        if C2_improved: status.append("C2↑")
        if charge_ok: status.append("Ni/Ne✓")

        print(f"  Status: {' '.join(status)}")

        return {
            'H': H_final, 'CH': CH_final, 'C2': C2_final,
            'C2H2': C2H2_final, 'CH3': CH3_final, 'Ni_over_Ne': Ni_over_Ne,
            'all_ok': all([H_ok, CH_ok, C2_improved, charge_ok])
        }

    except Exception as e:
        print(f"  Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("SYSTEMATIC C2H2 BOOST TEST")
    print("="*80)
    print(f"\nBaseline:")
    print(f"  H:  {baseline['target_densities']['H']:.2e} (79.9%)")
    print(f"  CH: {baseline['target_densities']['CH']:.2e} (100.7%)")
    print(f"  C2: {baseline['target_densities']['C2']:.2e} (16.8%)")
    print(f"  C2H2: {baseline['all_densities']['C2H2']:.2e}")
    print(f"  CH3:  {baseline['all_densities']['CH3']:.2e}")
    print(f"  Ni/Ne: {baseline['Ni_over_Ne']:.2f}")

    # Test 1: Boost e + CH4 → CH3 + H⁻ (main CH3 source)
    test_scenario(
        "Boost e + CH4 → CH3 by 2×",
        {'e_CH4_CH3_HMinus_cm3_8_1': 2.0}
    )

    # Test 2: Boost Ar* + CH4 → CH3 + H
    test_scenario(
        "Boost Ar* + CH4 → CH3 by 2×",
        {'ArStar_CH4_CH3_H_cm3_3_1': 2.0}
    )

    # Test 3: Boost both CH3 production pathways
    test_scenario(
        "Boost all CH3 production by 2×",
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 2.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 2.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 2.0,
        }
    )

    # Test 4: Reduce CH3 losses
    test_scenario(
        "Reduce CH3 wall losses by 50%",
        {'stick_CH3_9_2': 0.5}
    )

    # Test 5: Combined approach
    test_scenario(
        "Boost CH3 production 2×, reduce CH3 losses 50%",
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 2.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 2.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 2.0,
            'stick_CH3_9_2': 0.5,
        }
    )

    # Test 6: More aggressive
    test_scenario(
        "Boost CH3 production 3×, reduce CH3 losses 70%",
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 3.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 3.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 3.0,
            'stick_CH3_9_2': 0.3,
        }
    )

    # Test 7: Also reduce C2H2 losses
    test_scenario(
        "Boost CH3 3×, reduce losses for CH3 and C2H2",
        {
            'e_CH4_CH3_HMinus_cm3_8_1': 3.0,
            'ArStar_CH4_CH3_H_cm3_3_1': 3.0,
            'e_CH4Plus_CH3_H_cm3_6_4': 3.0,
            'stick_CH3_9_2': 0.3,
            'stick_C2H2_9_11': 0.5,
            'loss_C2H2_11_19': 0.5,
        }
    )

    print("\n" + "="*80)
    print("SYSTEMATIC TEST COMPLETE")
    print("="*80)
