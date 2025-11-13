"""
Find optimal Te that maximizes C2 without causing runaway

Higher Te boosts electron-impact reactions BUT causes ionization runaway.
Need to find the sweet spot: max C2 with Ni/Ne < 10 (stable)
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

def test_Te(Te_multiplier, P=500):
    """Test specific Te multiplier"""
    n_total = pressure_to_density(P)
    ne_frac = baseline['Ne'] / pressure_to_density(500.0)
    ne = ne_frac * n_total

    params = {
        'P': P,
        'Te': baseline['Te'] * Te_multiplier,
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
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    # Apply 200× CH3 boost + 99% loss + KEEP C2+H→CH+C (correct chemistry!)
    multipliers = {
        'e_CH4_CH3_HMinus_cm3_8_1': 200.0,
        'ArStar_CH4_CH3_H_cm3_3_1': 200.0,
        'e_CH4Plus_CH3_H_cm3_6_4': 200.0,
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
    }

    for rate_name, mult in multipliers.items():
        if rate_name in k:
            k[rate_name] *= mult

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
            return None

        y_final = sol.y[:, -1]
        H_final = y_final[species.index('H')]
        CH_final = y_final[species.index('CH')]
        C2_final = y_final[species.index('C2')]

        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus', 'CHPlus', 'H2Plus']
        n_i_total = sum(y_final[species.index(ion)] for ion in ions)
        ne_final = y_final[species.index('e')]
        Ni_over_Ne = n_i_total / ne_final if ne_final > 0 else 0

        C2_improvement = C2_final / baseline['target_densities']['C2']
        is_stable = (Ni_over_Ne < 10)  # Strict stability criterion

        return {
            'Te': baseline['Te'] * Te_multiplier,
            'C2': C2_final,
            'CH': CH_final,
            'H': H_final,
            'Ni_over_Ne': Ni_over_Ne,
            'improvement': C2_improvement,
            'stable': is_stable
        }

    except Exception as e:
        return None

if __name__ == '__main__':
    print("="*80)
    print("OPTIMIZE Te TO MAXIMIZE C2 WITHOUT RUNAWAY")
    print("="*80)
    print("\nGoal: Find max C2 with Ni/Ne < 10 (stable)")
    print("Constraint: Keep C2 + H → CH + C at literature value (correct chemistry)\n")

    # Sweep Te from 1.0× to 1.5× in small steps
    Te_multipliers = np.linspace(1.0, 1.5, 21)

    results = []

    print(f"{'Te mult':<10} {'Te (eV)':<10} {'C2 (%)':<12} {'Improvement':<12} {'Ni/Ne':<10} {'Status':<10}")
    print("-"*80)

    for Te_mult in Te_multipliers:
        result = test_Te(Te_mult)
        if result:
            c2_pct = result['C2'] / targets['C2'] * 100
            status = "✓ STABLE" if result['stable'] else "✗ RUNAWAY"
            print(f"{Te_mult:<10.3f} {result['Te']:<10.3f} {c2_pct:<12.2f} {result['improvement']:<12.1f} {result['Ni_over_Ne']:<10.2f} {status:<10}")
            results.append((Te_mult, result))
        else:
            print(f"{Te_mult:<10.3f} {'--':<10} {'FAILED':<12} {'--':<12} {'--':<10} {'✗':<10}")

    # Find best stable result
    stable_results = [(mult, r) for mult, r in results if r['stable']]

    if stable_results:
        best_mult, best_result = max(stable_results, key=lambda x: x[1]['C2'])

        print("\n" + "="*80)
        print("OPTIMAL Te FOUND!")
        print("="*80)
        print(f"\nBest stable configuration:")
        print(f"  Te multiplier: {best_mult:.3f}×")
        print(f"  Te: {best_result['Te']:.3f} eV (baseline: {baseline['Te']:.3f} eV)")
        print(f"  C2: {best_result['C2']:.2e} cm⁻³ ({best_result['C2']/targets['C2']*100:.2f}%)")
        print(f"  H:  {best_result['H']:.2e} cm⁻³ ({best_result['H']/targets['H']*100:.1f}%)")
        print(f"  CH: {best_result['CH']:.2e} cm⁻³ ({best_result['CH']/targets['CH']*100:.1f}%)")
        print(f"  Ni/Ne: {best_result['Ni_over_Ne']:.2f}")
        print(f"  C2 improvement: {best_result['improvement']:.1f}× from baseline")
        print(f"\n  Status: ✓ STABLE (Ni/Ne < 10)")

        # Compare to baseline
        baseline_c2 = 3.67e10  # from earlier test at Te=1.311 eV
        improvement_vs_baseline = best_result['C2'] / baseline_c2
        print(f"\n  Improvement vs baseline (Te={baseline['Te']:.3f} eV): {improvement_vs_baseline:.2f}×")
    else:
        print("\n✗ No stable configuration found! All tested Te values cause runaway.")

    print("\n" + "="*80)
