"""
Test baseline TUNED rates at different pressures

Use the 23 tuned rates from baseline (best_f70.3.json) at different pressures.
Test if lower pressure allows us to boost C2H2 production without instability.
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

# Load baseline with TUNED rates
with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

targets = {
    'H': 2.52e14,
    'CH': 1.0e9,
    'C2': 5.6e11,
}

def test_with_tuned_rates(pressure_mTorr, rate_multipliers=None, name="Baseline"):
    """Test chemistry at given pressure using baseline's tuned rates"""
    print(f"\n{'='*80}")
    print(f"{name} @ {pressure_mTorr} mTorr (WITH TUNED RATES)")
    print(f"{'='*80}")

    n_total = pressure_to_density(pressure_mTorr)

    # Scale Ne proportionally with pressure
    ne_frac = baseline['Ne'] / pressure_to_density(500.0)
    ne = ne_frac * n_total

    # Use same Te and E as baseline
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

    # Get default rates first
    k = define_rates(params)

    # Apply baseline's 23 TUNED rates
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    # Apply multipliers if given
    if rate_multipliers:
        print("Rate multipliers:")
        for rate_name, mult in rate_multipliers.items():
            if rate_name in k:
                k[rate_name] *= mult
                print(f"  {rate_name}: {k[rate_name]:.3e} (×{mult} from tuned value)")

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

        print(f"\nResults:")
        print(f"  H:      {H_final:.2e} ({H_final/targets['H']*100:6.1f}%) [baseline: 79.9%]")
        print(f"  CH:     {CH_final:.2e} ({CH_final/targets['CH']*100:6.1f}%) [baseline: 100.7%]")
        print(f"  C2:     {C2_final:.2e} ({C2_final/targets['C2']*100:6.1f}%) [baseline: 16.8%]")
        print(f"  C2H2:   {C2H2_final:.2e} ({C2H2_final/baseline['all_densities']['C2H2']:.2f}× baseline)")
        print(f"  CH3:    {CH3_final:.2e} ({CH3_final/baseline['all_densities']['CH3']:.2f}× baseline)")
        print(f"  Ni/Ne:  {Ni_over_Ne:.2f} [baseline: 3.12]")

        # Check stability
        is_stable = (H_final < 3 * targets['H'] and
                     Ni_over_Ne < 100 and
                     CH3_final < 1e15)

        # Check improvement
        c2_improved = C2_final > baseline['target_densities']['C2']

        if not is_stable:
            print(f"  ⚠️  RUNAWAY chemistry!")
        elif c2_improved:
            print(f"  ✓ Stable AND C2 IMPROVED!")
        else:
            print(f"  ✓ Stable (but C2 not improved)")

        return {
            'H': H_final, 'CH': CH_final, 'C2': C2_final,
            'C2H2': C2H2_final, 'CH3': CH3_final,
            'Ni_over_Ne': Ni_over_Ne,
            'stable': is_stable,
            'c2_improved': c2_improved,
        }

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("PRESSURE TEST WITH BASELINE'S TUNED RATES")
    print("="*80)
    print(f"\nBaseline @ 500 mTorr:")
    print(f"  H:  {baseline['target_densities']['H']:.2e} (79.9%)")
    print(f"  CH: {baseline['target_densities']['CH']:.2e} (100.7%)")
    print(f"  C2: {baseline['target_densities']['C2']:.2e} (16.8%)")
    print(f"  C2H2: {baseline['all_densities']['C2H2']:.2e}")
    print(f"  CH3: {baseline['all_densities']['CH3']:.2e}")
    print(f"  Ni/Ne: {baseline['Ni_over_Ne']:.2f}")

    # Test 1: Reproduce baseline at 500 mTorr
    print("\n" + "="*80)
    print("TEST 1: REPRODUCE BASELINE @ 500 mTorr")
    print("="*80)
    test_with_tuned_rates(500)

    # Test 2: Baseline tuned rates at lower pressures
    print("\n" + "="*80)
    print("TEST 2: BASELINE TUNED RATES AT LOWER PRESSURES")
    print("="*80)
    for P in [300, 400]:
        test_with_tuned_rates(P)

    # Test 3: Small CH3 boost at different pressures
    print("\n" + "="*80)
    print("TEST 3: SMALL CH3 BOOST (1.5×) AT DIFFERENT PRESSURES")
    print("="*80)
    for P in [300, 400, 500]:
        test_with_tuned_rates(
            P,
            {
                'e_CH4_CH3_HMinus_cm3_8_1': 1.5,
                'ArStar_CH4_CH3_H_cm3_3_1': 1.5,
            },
            name="1.5× CH3 boost"
        )

    # Test 4: Moderate CH3 boost
    print("\n" + "="*80)
    print("TEST 4: MODERATE CH3 BOOST (2×) AT DIFFERENT PRESSURES")
    print("="*80)
    for P in [300, 400, 500]:
        test_with_tuned_rates(
            P,
            {
                'e_CH4_CH3_HMinus_cm3_8_1': 2.0,
                'ArStar_CH4_CH3_H_cm3_3_1': 2.0,
                'e_CH4Plus_CH3_H_cm3_6_4': 2.0,
            },
            name="2× CH3 boost"
        )

    # Test 5: Aggressive approach - boost CH3 AND reduce losses
    print("\n" + "="*80)
    print("TEST 5: BOOST CH3 3× + REDUCE CH3/C2H2 LOSSES 50%")
    print("="*80)
    for P in [300, 400, 500]:
        test_with_tuned_rates(
            P,
            {
                'e_CH4_CH3_HMinus_cm3_8_1': 3.0,
                'ArStar_CH4_CH3_H_cm3_3_1': 3.0,
                'e_CH4Plus_CH3_H_cm3_6_4': 3.0,
                'stick_CH3_9_2': 0.5,
                'stick_C2H2_9_11': 0.5,
                'loss_C2H2_11_19': 0.5,
            },
            name="3× CH3, 50% losses"
        )

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print("Check which pressure + multipliers give:")
    print("  ✓ Stable chemistry (no runaway)")
    print("  ✓ C2 improved over baseline's 16.8%")
    print("="*80)
