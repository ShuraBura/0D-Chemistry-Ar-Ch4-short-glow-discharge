"""
Questions 2 & 3:
2) Can we tune CH3 + CH3 + M → C2H6 → C2H2?
3) Test higher Te and LOWER pressure

Strategy:
- Higher Te: boosts electron-impact reactions (e + C2H6 → C2H2, e + CH4 → CH3, etc.)
- LOWER pressure: reduces three-body rate (CH3+CH3+M), but also reduces losses
- Boost C2H6 → C2H2 conversion reactions
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

def test_conditions(P, Te_multiplier, boost_C2H6_conversion, rate_mults, name):
    """Test different pressure, Te, and C2H6→C2H2 conversion"""
    print(f"\n{'='*80}")
    print(f"{name}")
    print(f"{'='*80}")
    print(f"Pressure: {P} mTorr")
    print(f"Te: {baseline['Te'] * Te_multiplier:.3f} eV (×{Te_multiplier})")
    print(f"Boost C2H6→C2H2: {boost_C2H6_conversion}×")

    n_total = pressure_to_density(P)
    ne_frac = baseline['Ne'] / pressure_to_density(500.0)
    ne = ne_frac * n_total

    params = {
        'P': P,
        'Te': baseline['Te'] * Te_multiplier,  # Modify Te
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

    # Boost C2H6 → C2H2 conversion
    if boost_C2H6_conversion > 1:
        conversion_reactions = [
            'e_C2H6_C2H2_2H2_cm3_1_18',      # e + C2H6 → C2H2 + 2H2
            'CH_C2H6_C2H2_CH3_H_cm3_7_57',   # CH + C2H6 → C2H2 + CH3 + H
            'e_C2H6_C2H4_H2_cm3_1_8',        # e + C2H6 → C2H4 + H2
            'e_C2H4_C2H2_H2_cm3_1_6',        # e + C2H4 → C2H2 + H2
        ]
        print(f"\nBoosting C2H6→C2H2 conversion reactions by {boost_C2H6_conversion}×")
        for rxn in conversion_reactions:
            if rxn in k:
                k[rxn] *= boost_C2H6_conversion

    # Apply test multipliers
    if rate_mults:
        print("\nRate multipliers:")
        for rate_name, mult in rate_mults.items():
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
        sol = solve_ivp(ode_func, (0, 500), y0, method='BDF', rtol=1e-7, atol=1e-9, max_step=1.0)

        if not sol.success:
            print(f"  ✗ Integration failed!")
            return None

        y_final = sol.y[:, -1]
        H_final = y_final[species.index('H')]
        CH_final = y_final[species.index('CH')]
        C2_final = y_final[species.index('C2')]
        C2H2_final = y_final[species.index('C2H2')]
        C2H6_final = y_final[species.index('C2H6')]
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
        print(f"  C2H2: {C2H2_final:.2e}")
        print(f"  C2H6: {C2H6_final:.2e}")
        print(f"  CH3:  {CH3_final:.2e}")
        print(f"  Ni/Ne: {Ni_over_Ne:.2f}")
        print(f"  C2 improvement: {C2_improvement:.1f}×")

        is_stable = (H_final < 5 * targets['H'] and Ni_over_Ne < 100)

        if not is_stable:
            print(f"  ⚠️  RUNAWAY!")
        else:
            print(f"  ✓ STABLE")

        return {
            'C2': C2_final, 'CH': CH_final, 'H': H_final,
            'C2H2': C2H2_final, 'C2H6': C2H6_final,
            'improvement': C2_improvement,
            'Ni_over_Ne': Ni_over_Ne,
            'stable': is_stable
        }

    except Exception as e:
        print(f"  ✗ Error: {e}")
        return None

if __name__ == '__main__':
    print("="*80)
    print("QUESTIONS 2 & 3: TUNE C2H6→C2H2 + HIGHER Te + LOWER PRESSURE")
    print("="*80)

    # Base multipliers: 200× CH3 + 99% loss + KEEP C2+H→CH+C at full rate (user correction!)
    base_mults = {
        'e_CH4_CH3_HMinus_cm3_8_1': 200.0,
        'ArStar_CH4_CH3_H_cm3_3_1': 200.0,
        'e_CH4Plus_CH3_H_cm3_6_4': 200.0,
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
        # C2_H_CH_C_cm3_7_6 NOT modified - keep at literature value!
    }

    results = []

    # Baseline: 500 mTorr, Te = 1.475 eV
    r1 = test_conditions(500, 1.0, 1, base_mults, "BASELINE: 500 mTorr, Te=1.475 eV")
    results.append(("Baseline (500 mTorr)", r1))

    # Test 1: LOWER pressure (300 mTorr, same Te)
    print("\n" + "="*80)
    print("TEST 1: LOWER PRESSURE")
    print("="*80)
    r2 = test_conditions(300, 1.0, 1, base_mults, "300 mTorr, Te=1.475 eV")
    results.append(("Lower P (300 mTorr)", r2))

    # Test 2: HIGHER Te (500 mTorr, Te×1.5)
    print("\n" + "="*80)
    print("TEST 2: HIGHER Te")
    print("="*80)
    r3 = test_conditions(500, 1.5, 1, base_mults, "500 mTorr, Te=2.21 eV")
    results.append(("Higher Te (×1.5)", r3))

    # Test 3: Boost C2H6→C2H2 conversion (500 mTorr, same Te)
    print("\n" + "="*80)
    print("TEST 3: BOOST C2H6→C2H2 CONVERSION")
    print("="*80)
    r4 = test_conditions(500, 1.0, 10, base_mults, "500 mTorr, Te=1.475 eV, 10× C2H6→C2H2")
    results.append(("Boost C2H6→C2H2 (10×)", r4))

    # Test 4: COMBINED - Lower P + Higher Te
    print("\n" + "="*80)
    print("TEST 4: COMBINED - LOWER P + HIGHER Te")
    print("="*80)
    r5 = test_conditions(300, 1.5, 1, base_mults, "300 mTorr, Te=2.21 eV")
    results.append(("Lower P + Higher Te", r5))

    # Test 5: COMBINED - Lower P + Higher Te + Boost C2H6→C2H2
    print("\n" + "="*80)
    print("TEST 5: COMBINED - ALL THREE")
    print("="*80)
    r6 = test_conditions(300, 1.5, 10, base_mults, "300 mTorr, Te=2.21 eV, 10× C2H6→C2H2")
    results.append(("ALL: Lower P + Higher Te + Boost", r6))

    # Test 6: Very high Te
    print("\n" + "="*80)
    print("TEST 6: VERY HIGH Te")
    print("="*80)
    r7 = test_conditions(300, 2.0, 10, base_mults, "300 mTorr, Te=2.95 eV, 10× C2H6→C2H2")
    results.append(("Very high Te (×2.0)", r7))

    # Summary
    print("\n" + "="*80)
    print("SUMMARY: HIGHER Te + LOWER PRESSURE + C2H6→C2H2 TUNING")
    print("="*80)

    print(f"\n{'Configuration':<30} {'C2 (%)':<10} {'C2H2':<12} {'Improvement':<12} {'Stable?':<10}")
    print("-"*80)

    for name, result in results:
        if result and result['stable']:
            c2_pct = result['C2'] / targets['C2'] * 100
            print(f"{name:<30} {c2_pct:>8.2f}% {result['C2H2']:>10.2e} {result['improvement']:>10.1f}× {'✓':<10}")
        elif result:
            c2_pct = result['C2'] / targets['C2'] * 100
            print(f"{name:<30} {c2_pct:>8.2f}% {result['C2H2']:>10.2e} {result['improvement']:>10.1f}× {'✗ RUNAWAY':<10}")
        else:
            print(f"{name:<30} {'FAILED':<10} {'--':<12} {'--':<12} {'✗':<10}")

    print("\n" + "="*80)
    print("Can lower pressure + higher Te push C2 beyond 22.78%?")
    print("="*80)
