"""
Test plasma chemistry at different pressures

Hypothesis: Lower pressure → more stable chemistry

Test pressures: 300, 400, 500 mTorr
For each:
1. Run baseline conditions (Te, Ne, E from 500 mTorr baseline scaled appropriately)
2. Test small perturbations (2× CH3 production)
3. Check if stability improves at lower pressure
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

# Load baseline (at 500 mTorr)
with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline_500 = json.load(f)

targets = {
    'H': 2.52e14,
    'CH': 1.0e9,
    'C2': 5.6e11,
}

def test_pressure(pressure_mTorr, rate_multipliers=None, name="Baseline"):
    """Test chemistry at given pressure"""
    print(f"\n{'='*80}")
    print(f"{name} @ {pressure_mTorr} mTorr")
    print(f"{'='*80}")

    n_total = pressure_to_density(pressure_mTorr)

    # Scale Ne proportionally with pressure (keep ne/n_total constant)
    ne_frac = baseline_500['Ne'] / pressure_to_density(500.0)
    ne = ne_frac * n_total

    # Use same Te and E as baseline
    Te = baseline_500['Te']
    E_field = baseline_500['E_field']

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

    # Apply rate multipliers if given
    if rate_multipliers:
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

        print(f"\nPlasma conditions:")
        print(f"  n_total: {n_total:.2e} cm⁻³")
        print(f"  Te:      {Te:.2f} eV")
        print(f"  Ne:      {ne:.2e} cm⁻³ ({ne/n_total*1e6:.2f} ppm)")
        print(f"  E:       {E_field:.1f} V/cm")

        print(f"\nResults:")
        print(f"  H:      {H_final:.2e} ({H_final/targets['H']*100:6.1f}%)")
        print(f"  CH:     {CH_final:.2e} ({CH_final/targets['CH']*100:6.1f}%)")
        print(f"  C2:     {C2_final:.2e} ({C2_final/targets['C2']*100:6.1f}%)")
        print(f"  C2H2:   {C2H2_final:.2e}")
        print(f"  CH3:    {CH3_final:.2e}")
        print(f"  Ni/Ne:  {Ni_over_Ne:.2f}")

        # Check for runaway (unphysical)
        is_runaway = (H_final > 10 * targets['H'] or
                      Ni_over_Ne > 100 or
                      CH3_final > 1e15)

        if is_runaway:
            print(f"  ⚠️  RUNAWAY chemistry detected!")
        else:
            print(f"  ✓ Stable chemistry")

        return {
            'H': H_final, 'CH': CH_final, 'C2': C2_final,
            'C2H2': C2H2_final, 'CH3': CH3_final,
            'Ni_over_Ne': Ni_over_Ne,
            'runaway': is_runaway,
            'n_total': n_total,
        }

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("PRESSURE STABILITY TEST")
    print("="*80)
    print(f"\nBaseline @ 500 mTorr:")
    print(f"  H:  {baseline_500['target_densities']['H']:.2e} (79.9%)")
    print(f"  CH: {baseline_500['target_densities']['CH']:.2e} (100.7%)")
    print(f"  C2: {baseline_500['target_densities']['C2']:.2e} (16.8%)")
    print(f"  Ni/Ne: {baseline_500['Ni_over_Ne']:.2f}")

    results = {}

    # Test 1: Baseline at each pressure
    print("\n" + "="*80)
    print("TEST 1: BASELINE CONDITIONS AT DIFFERENT PRESSURES")
    print("="*80)

    for P in [300, 400, 500]:
        results[f'baseline_{P}'] = test_pressure(P)

    # Test 2: Small perturbation (1.5× CH3 production) - more conservative
    print("\n" + "="*80)
    print("TEST 2: MODEST CH3 BOOST (1.5×) AT DIFFERENT PRESSURES")
    print("="*80)

    for P in [300, 400, 500]:
        results[f'boost15_{P}'] = test_pressure(
            P,
            {
                'e_CH4_CH3_HMinus_cm3_8_1': 1.5,
                'ArStar_CH4_CH3_H_cm3_3_1': 1.5,
            },
            name="1.5× CH3 boost"
        )

    # Test 3: Moderate perturbation (2× CH3 production)
    print("\n" + "="*80)
    print("TEST 3: 2× CH3 BOOST AT DIFFERENT PRESSURES")
    print("="*80)

    for P in [300, 400, 500]:
        results[f'boost20_{P}'] = test_pressure(
            P,
            {
                'e_CH4_CH3_HMinus_cm3_8_1': 2.0,
                'ArStar_CH4_CH3_H_cm3_3_1': 2.0,
            },
            name="2× CH3 boost"
        )

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    for test_name in ['baseline', 'boost15', 'boost20']:
        print(f"\n{test_name.upper()}:")
        for P in [300, 400, 500]:
            key = f'{test_name}_{P}'
            if results.get(key) and not results[key].get('runaway'):
                print(f"  {P} mTorr: ✓ Stable")
            elif results.get(key):
                print(f"  {P} mTorr: ✗ Runaway")
            else:
                print(f"  {P} mTorr: ✗ Failed")

    print("\n" + "="*80)
    print("CONCLUSION")
    print("="*80)
    print("If lower pressure shows stable chemistry with boosts → use lower pressure!")
    print("="*80)
