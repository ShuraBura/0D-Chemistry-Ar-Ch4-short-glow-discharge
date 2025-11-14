#!/usr/bin/env python3
"""
Explore three hypotheses to match all three species (H, CH, C2):
1. Higher Te (2.5-4 eV) - changes electron impact rates
2. Boost C2 production: e + C2H2 → C2 + H2 (maybe too slow)
3. Suppress C2+H→CH+C (maybe too fast for vibrationally excited C2)
"""

import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))

from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized
from define_rates import define_rates
from scipy.integrate import solve_ivp

# CORRECTED Target densities
TARGET_H = 2.3e14    # cm^-3
TARGET_CH = 1.34e9   # cm^-3
TARGET_C2 = 5.60e11  # cm^-3

# Parameters
PRESSURE_MTORR = 400
TGAS_K = 570
NE_NOMINAL = 2.3e9

def pressure_to_density(pressure_mTorr, T_K):
    """Convert pressure to number density"""
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

def run_simulation(ne, Te_eV, P_mTorr, E_field,
                   c2h_suppression=1.0,  # Factor to reduce C2+H→CH+C
                   c2_production_boost=1.0,  # Factor to boost e+C2H2→C2
                   H_drift_gain_factor=1.0):
    """Run simulation with rate modifications"""

    n_total = pressure_to_density(P_mTorr, TGAS_K)

    params = {
        'P': P_mTorr,
        'Tgas': TGAS_K,
        'Te': Te_eV,
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
                    'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4',
                    'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3',
                    'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3',
                    'C3H5', 'C4H', 'C3H6', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                    'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
    }

    k = define_rates(params)

    # Apply modifications
    # 1. C2+H→CH+C suppression (Hypothesis 3)
    if 'C2_H_CH_C_cm3_7_6' in k:
        k['C2_H_CH_C_cm3_7_6'] *= c2h_suppression

    # 2. C2 production boost (Hypothesis 2)
    # Look for e + C2H2 reactions
    c2_production_keys = [
        'e_C2H2_C2_H2_cm3',
        'e_C2H2_C2_H_H_cm3_8_18',
        'e_C2H2_C2H_H_cm3_8_16',
    ]
    for key in c2_production_keys:
        if key in k:
            k[key] *= c2_production_boost

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    species = params['species']
    y0 = np.ones(len(species)) * 1e3
    y0[species.index('Ar')] = n_total * 0.97
    y0[species.index('CH4')] = n_total * 0.03
    y0[species.index('e')] = ne

    try:
        ode_func = PlasmaODE_Optimized(params)
        ode_func.H_drift_gain = 5.7e16 * H_drift_gain_factor

        sol = solve_ivp(ode_func, (0, 500), y0, method='BDF',
                       rtol=1e-7, atol=1e-9, max_step=1.0)

        if not sol.success:
            return None

        y_final = sol.y[:, -1]
        return {
            'H': y_final[species.index('H')],
            'CH': y_final[species.index('CH')],
            'C2': y_final[species.index('C2')]
        }

    except Exception as e:
        return None

def calculate_rms_error(H, CH, C2):
    """Calculate RMS percentage error"""
    h_err = abs(H / TARGET_H - 1.0) * 100
    ch_err = abs(CH / TARGET_CH - 1.0) * 100
    c2_err = abs(C2 / TARGET_C2 - 1.0) * 100
    return np.sqrt((h_err**2 + ch_err**2 + c2_err**2) / 3)

def main():
    print("="*80)
    print("EXPLORING THREE HYPOTHESES TO MATCH ALL THREE SPECIES")
    print("="*80)
    print(f"\nTargets: H={TARGET_H:.2e}, CH={TARGET_CH:.2e}, C2={TARGET_C2:.2e}")
    print()

    # =========================================================================
    # HYPOTHESIS 1: Higher Te (2.5-4 eV)
    # =========================================================================
    print("\n" + "="*80)
    print("HYPOTHESIS 1: Higher Te changes electron impact rates")
    print("="*80)

    print("\nTesting Te range 1.2-4.0 eV at ne=2.3e9, P=400, E=150...")

    Te_values = np.array([1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
    h1_results = []

    for Te in Te_values:
        res = run_simulation(NE_NOMINAL, Te, 400, 150)
        if res:
            rms = calculate_rms_error(res['H'], res['CH'], res['C2'])
            h1_results.append({
                'Te': Te,
                'H': res['H'], 'CH': res['CH'], 'C2': res['C2'],
                'H_pct': res['H']/TARGET_H*100,
                'CH_pct': res['CH']/TARGET_CH*100,
                'C2_pct': res['C2']/TARGET_C2*100,
                'rms': rms
            })
            print(f"  Te={Te:.1f} eV: H={res['H']/TARGET_H*100:6.1f}%, "
                  f"CH={res['CH']/TARGET_CH*100:7.0f}%, "
                  f"C2={res['C2']/TARGET_C2*100:6.1f}%, RMS={rms:6.1f}%")

    if h1_results:
        best_h1 = min(h1_results, key=lambda x: x['rms'])
        print(f"\n>>> Best at Te={best_h1['Te']:.1f} eV: RMS={best_h1['rms']:.1f}%")

    # =========================================================================
    # HYPOTHESIS 2: Boost C2 production (e + C2H2 → C2)
    # =========================================================================
    print("\n" + "="*80)
    print("HYPOTHESIS 2: C2 production rate might be too slow")
    print("="*80)

    print("\nTesting C2 production boost factors at ne=2.3e9, Te=1.5eV, P=400, E=150...")

    boost_factors = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    h2_results = []

    for boost in boost_factors:
        res = run_simulation(NE_NOMINAL, 1.5, 400, 150, c2_production_boost=boost)
        if res:
            rms = calculate_rms_error(res['H'], res['CH'], res['C2'])
            h2_results.append({
                'boost': boost,
                'H': res['H'], 'CH': res['CH'], 'C2': res['C2'],
                'H_pct': res['H']/TARGET_H*100,
                'CH_pct': res['CH']/TARGET_CH*100,
                'C2_pct': res['C2']/TARGET_C2*100,
                'rms': rms
            })
            print(f"  Boost={boost:5.0f}x: H={res['H']/TARGET_H*100:6.1f}%, "
                  f"CH={res['CH']/TARGET_CH*100:7.0f}%, "
                  f"C2={res['C2']/TARGET_C2*100:6.1f}%, RMS={rms:6.1f}%")

    if h2_results:
        best_h2 = min(h2_results, key=lambda x: x['rms'])
        print(f"\n>>> Best at boost={best_h2['boost']:.0f}x: RMS={best_h2['rms']:.1f}%")

    # =========================================================================
    # HYPOTHESIS 3: Suppress C2+H→CH+C (wrong for vibrationally excited C2)
    # =========================================================================
    print("\n" + "="*80)
    print("HYPOTHESIS 3: C2+H→CH+C might be too fast for vibrationally excited C2")
    print("="*80)

    print("\nTesting C2+H suppression factors at ne=2.3e9, Te=1.5eV, P=400, E=150...")

    suppress_factors = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
    h3_results = []

    for suppress in suppress_factors:
        res = run_simulation(NE_NOMINAL, 1.5, 400, 150, c2h_suppression=suppress)
        if res:
            rms = calculate_rms_error(res['H'], res['CH'], res['C2'])
            h3_results.append({
                'suppress': suppress,
                'H': res['H'], 'CH': res['CH'], 'C2': res['C2'],
                'H_pct': res['H']/TARGET_H*100,
                'CH_pct': res['CH']/TARGET_CH*100,
                'C2_pct': res['C2']/TARGET_C2*100,
                'rms': rms
            })
            print(f"  Suppress={suppress:5.2f}x: H={res['H']/TARGET_H*100:6.1f}%, "
                  f"CH={res['CH']/TARGET_CH*100:7.0f}%, "
                  f"C2={res['C2']/TARGET_C2*100:6.1f}%, RMS={rms:6.1f}%")

    if h3_results:
        best_h3 = min(h3_results, key=lambda x: x['rms'])
        print(f"\n>>> Best at suppress={best_h3['suppress']:.2f}x: RMS={best_h3['rms']:.1f}%")

    # =========================================================================
    # COMBINED: Test best combinations
    # =========================================================================
    print("\n" + "="*80)
    print("COMBINED: Testing combinations of all three hypotheses")
    print("="*80)

    print("\nTesting combined factors...")

    # Use insights from individual tests to guide combined search
    combined_tests = []

    # Test grid around promising regions
    for Te in [1.5, 2.0, 2.5, 3.0]:
        for ne_factor in [0.8, 1.0, 1.2, 1.5]:
            ne = NE_NOMINAL * ne_factor
            for c2h_suppress in [0.1, 0.05, 0.02, 0.01]:
                for c2_boost in [5.0, 10.0, 20.0, 50.0]:
                    res = run_simulation(ne, Te, 400, 150,
                                       c2h_suppression=c2h_suppress,
                                       c2_production_boost=c2_boost)
                    if res:
                        rms = calculate_rms_error(res['H'], res['CH'], res['C2'])
                        combined_tests.append({
                            'ne': ne, 'Te': Te,
                            'c2h_suppress': c2h_suppress,
                            'c2_boost': c2_boost,
                            'H': res['H'], 'CH': res['CH'], 'C2': res['C2'],
                            'H_pct': res['H']/TARGET_H*100,
                            'CH_pct': res['CH']/TARGET_CH*100,
                            'C2_pct': res['C2']/TARGET_C2*100,
                            'rms': rms
                        })

    print(f"\n  Tested {len(combined_tests)} combinations")

    if combined_tests:
        combined_tests.sort(key=lambda x: x['rms'])

        print("\n" + "="*80)
        print("TOP 10 COMBINED RESULTS")
        print("="*80)

        for i, r in enumerate(combined_tests[:10]):
            print(f"\n--- Rank {i+1} (RMS={r['rms']:.1f}%) ---")
            print(f"Parameters:")
            print(f"  ne={r['ne']:.2e}, Te={r['Te']:.1f} eV")
            print(f"  C2+H suppression={r['c2h_suppress']:.2f}x")
            print(f"  C2 production boost={r['c2_boost']:.0f}x")
            print(f"Results:")
            print(f"  H:  {r['H']:.2e} ({r['H_pct']:6.1f}%)")
            print(f"  CH: {r['CH']:.2e} ({r['CH_pct']:6.0f}%)")
            print(f"  C2: {r['C2']:.2e} ({r['C2_pct']:6.1f}%)")

        print("\n" + "="*80)
        print("BEST OVERALL RESULT")
        print("="*80)

        best = combined_tests[0]
        print(f"\nParameters:")
        print(f"  ne = {best['ne']:.2e} cm^-3")
        print(f"  Te = {best['Te']:.1f} eV")
        print(f"  P = 400 mTorr")
        print(f"  E = 150 V/cm")
        print(f"  C2+H→CH+C suppression = {best['c2h_suppress']:.2f}x")
        print(f"  C2 production boost = {best['c2_boost']:.0f}x")
        print(f"\nResults:")
        print(f"  H:  {best['H']:.2e} cm^-3 ({best['H_pct']:6.1f}% of target)")
        print(f"  CH: {best['CH']:.2e} cm^-3 ({best['CH_pct']:6.0f}% of target)")
        print(f"  C2: {best['C2']:.2e} cm^-3 ({best['C2_pct']:6.1f}% of target)")
        print(f"  RMS error: {best['rms']:.1f}%")

        print(f"\n>>> INTERPRETATION <<<")
        if best['c2h_suppress'] < 0.2:
            print(f"  C2+H→CH+C rate needs {1/best['c2h_suppress']:.0f}x REDUCTION")
            print(f"  This suggests vibrationally excited C2 is LESS reactive")
        if best['c2_boost'] > 5:
            print(f"  C2 production needs {best['c2_boost']:.0f}x INCREASE")
            print(f"  This suggests e+C2H2→C2 rate is too slow in model")
        if abs(best['ne']/NE_NOMINAL - 1.0) > 0.2:
            print(f"  ne should be {best['ne']/NE_NOMINAL:.2f}x measured value")
        if abs(best['Te'] - 1.5) > 0.5:
            print(f"  Te = {best['Te']:.1f} eV is significantly different from assumed 1.3-1.5 eV")

if __name__ == '__main__':
    main()
