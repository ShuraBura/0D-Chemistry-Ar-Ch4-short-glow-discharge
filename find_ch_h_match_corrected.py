#!/usr/bin/env python3
"""
Find parameter configuration where CH and H match well with CORRECTED parameters.
Then use C2 discrepancy to back-calculate Tv.
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
TARGET_H = 2.3e14    # cm^-3 (CORRECTED from 2.52e14)
TARGET_CH = 1.34e9   # cm^-3 (CORRECTED from 1.0e9)
TARGET_C2 = 5.60e11  # cm^-3

# CORRECTED experimental parameters with ranges
PRESSURE_MTORR = 400  # mTorr (CORRECTED from 500!)
P_RANGE = (380, 440)  # Can vary +40/-20 mTorr
NE_NOMINAL = 2.3e9    # cm^-3
NE_RANGE = (1.15e9, 4.6e9)  # Can vary ~100%
TE_EV_NOMINAL = 1.3   # eV, but can vary a lot (especially higher)
TE_RANGE = (1.0, 3.0)  # Explore higher Te values
E_FIELD_RANGE = (50, 300)  # V/cm
TGAS_K = 570          # K

def pressure_to_density(pressure_mTorr, T_K):
    """Convert pressure to number density"""
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6  # Convert to cm^-3

def c2_vibrational_fraction(Tv, v=1):
    """Calculate fraction of C2 in vibrational level v at temperature Tv"""
    omega_e = 1641  # cm^-1 for C2
    delta_E_K = 1.4388 * omega_e  # Convert to Kelvin (2361 K)
    Z = sum(np.exp(-vi * delta_E_K / Tv) for vi in range(10))
    return np.exp(-v * delta_E_K / Tv) / Z

def calculate_required_Tv(measured_C2, predicted_C2, Tv_assumed=1500):
    """Back-calculate what Tv would make measured C2 equal predicted C2"""
    f_v1_assumed = c2_vibrational_fraction(Tv_assumed, v=1)
    target_f_v1 = f_v1_assumed * measured_C2 / predicted_C2

    for Tv_test in np.linspace(500, 5000, 1000):
        f_v1_test = c2_vibrational_fraction(Tv_test, v=1)
        if abs(f_v1_test - target_f_v1) / max(target_f_v1, 1e-10) < 0.01:
            return Tv_test
    return None

def run_simulation(ne, Te_eV, P_mTorr, E_field, c2h_factor=1.0, chh_factor=1.0,
                   dust_density=0, H_drift_gain_factor=1.0):
    """Run simulation with given parameters"""

    n_total = pressure_to_density(P_mTorr, TGAS_K)

    params = {
        'P': P_mTorr,
        'Tgas': TGAS_K,
        'Te': Te_eV,
        'ne': ne,
        'E_field': E_field,
        'L_discharge': 0.45,  # cm
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

    if dust_density > 0:
        params['enable_dust_loss'] = True
        params['dust_density'] = dust_density
        params['dust_radius'] = 50e-7
        params['dust_sticking'] = 0.5

    k = define_rates(params)

    # Apply rate modifications
    if 'C2_H_CH_C_cm3_7_6' in k:
        k['C2_H_CH_C_cm3_7_6'] *= c2h_factor
    if 'CH_H_C_H2_cm3_7_3' in k:
        k['CH_H_C_H2_cm3_7_3'] *= chh_factor

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

def main():
    print("="*80)
    print("FINDING CONFIGURATION WHERE CH AND H MATCH WELL")
    print("="*80)
    print(f"\nCORRECTED Target densities:")
    print(f"  H:  {TARGET_H:.2e} cm^-3")
    print(f"  CH: {TARGET_CH:.2e} cm^-3")
    print(f"  C2: {TARGET_C2:.2e} cm^-3")
    print(f"\nVariable parameter ranges:")
    print(f"  P:  {P_RANGE[0]}-{P_RANGE[1]} mTorr")
    print(f"  ne: {NE_RANGE[0]:.2e}-{NE_RANGE[1]:.2e} cm^-3")
    print(f"  Te: {TE_RANGE[0]}-{TE_RANGE[1]} eV")
    print(f"  E:  {E_FIELD_RANGE[0]}-{E_FIELD_RANGE[1]} V/cm")
    print(f"  Tgas: {TGAS_K} K (fixed)")
    print()

    # Parameter sweep - focus on key variables
    ne_values = np.array([1.5e9, 2.0e9, 2.3e9, 3.0e9, 4.0e9])
    Te_values = np.array([1.2, 1.5, 2.0, 2.5])
    P_values = np.array([380, 400, 420, 440])
    E_values = np.array([50, 100, 150, 200, 250])

    # Rate modification factors
    c2h_factors = [0.5, 1.0]
    chh_factors = [1.0, 2.0]
    h_drift_factors = [0.5, 1.0]

    results = []

    # Start with a coarse sweep on Te and ne (most important)
    print("Phase 1: Sweeping Te and ne (E=150, P=400)...")
    for Te in Te_values:
        for ne in ne_values:
            for h_drift in h_drift_factors:
                res = run_simulation(ne, Te, 400, 150,
                                    c2h_factor=1.0, chh_factor=1.0,
                                    H_drift_gain_factor=h_drift)
                if res and res['H'] > 1e13 and res['CH'] > 1e8:
                    H_ratio = res['H'] / TARGET_H
                    CH_ratio = res['CH'] / TARGET_CH
                    C2_ratio = res['C2'] / TARGET_C2

                    ch_h_match = abs(np.log10(H_ratio) - np.log10(CH_ratio))
                    avg_match = np.sqrt(H_ratio * CH_ratio)
                    score = ch_h_match + abs(np.log10(avg_match))

                    results.append({
                        'ne': ne, 'Te': Te, 'P': 400, 'E': 150,
                        'c2h': 1.0, 'chh': 1.0, 'h_drift': h_drift,
                        'H': res['H'], 'CH': res['CH'], 'C2': res['C2'],
                        'H_ratio': H_ratio, 'CH_ratio': CH_ratio, 'C2_ratio': C2_ratio,
                        'score': score, 'ch_h_match': ch_h_match
                    })

    print(f"  Found {len(results)} valid configurations")

    if len(results) == 0:
        print("ERROR: No successful simulations in Phase 1!")
        return

    # Sort and take top 5
    results.sort(key=lambda x: x['score'])
    top5 = results[:5]

    print(f"\nPhase 2: Refining top configurations with E and P variations...")
    refined_results = list(top5)  # Keep original top 5

    for config in top5[:3]:  # Refine top 3
        for E in E_values:
            for P in P_values:
                res = run_simulation(config['ne'], config['Te'], P, E,
                                    H_drift_gain_factor=config['h_drift'])
                if res and res['H'] > 1e13 and res['CH'] > 1e8:
                    H_ratio = res['H'] / TARGET_H
                    CH_ratio = res['CH'] / TARGET_CH
                    C2_ratio = res['C2'] / TARGET_C2

                    ch_h_match = abs(np.log10(H_ratio) - np.log10(CH_ratio))
                    avg_match = np.sqrt(H_ratio * CH_ratio)
                    score = ch_h_match + abs(np.log10(avg_match))

                    refined_results.append({
                        'ne': config['ne'], 'Te': config['Te'], 'P': P, 'E': E,
                        'c2h': 1.0, 'chh': 1.0, 'h_drift': config['h_drift'],
                        'H': res['H'], 'CH': res['CH'], 'C2': res['C2'],
                        'H_ratio': H_ratio, 'CH_ratio': CH_ratio, 'C2_ratio': C2_ratio,
                        'score': score, 'ch_h_match': ch_h_match
                    })

    refined_results.sort(key=lambda x: x['score'])

    print(f"\n{'='*80}")
    print("TOP 10 CONFIGURATIONS WHERE CH AND H MATCH BEST")
    print(f"{'='*80}\n")

    for i, r in enumerate(refined_results[:10]):
        print(f"--- Rank {i+1} ---")
        print(f"Parameters:")
        print(f"  ne = {r['ne']:.2e} cm^-3,  Te = {r['Te']:.2f} eV")
        print(f"  P = {r['P']:.0f} mTorr,  E = {r['E']:.0f} V/cm")
        print(f"  H_drift = {r['h_drift']:.2f} x 5.7e16")
        print(f"Results:")
        print(f"  H:  {r['H']:.2e} ({r['H_ratio']*100:5.1f}% of target)")
        print(f"  CH: {r['CH']:.2e} ({r['CH_ratio']*100:5.1f}% of target)")
        print(f"  C2: {r['C2']:.2e} ({r['C2_ratio']*100:5.1f}% of target)")
        print(f"Match score: {r['ch_h_match']:.3f}\n")

        # Calculate required Tv
        if r['C2'] > 1e10:
            Tv_req = calculate_required_Tv(TARGET_C2, r['C2'], 1500)
            if Tv_req and 800 < Tv_req < 4500:
                print(f">>> Tv correction: {Tv_req:.0f} K (assumed 1500K in LIF)")
                print()

    print(f"{'='*80}")
    print("BEST MATCH SUMMARY")
    print(f"{'='*80}\n")
    best = refined_results[0]
    print(f"Parameters:")
    print(f"  ne = {best['ne']:.2e} cm^-3")
    print(f"  Te = {best['Te']:.2f} eV")
    print(f"  P = {best['P']:.0f} mTorr")
    print(f"  E = {best['E']:.0f} V/cm")
    print(f"  H_drift_gain = {best['h_drift']:.2f} x 5.7e16 cm^-3/s")
    print(f"\nResults:")
    print(f"  H:  {best['H']:.2e} cm^-3 ({best['H_ratio']*100:.1f}% of {TARGET_H:.2e})")
    print(f"  CH: {best['CH']:.2e} cm^-3 ({best['CH_ratio']*100:.1f}% of {TARGET_CH:.2e})")
    print(f"  C2: {best['C2']:.2e} cm^-3 ({best['C2_ratio']*100:.1f}% of {TARGET_C2:.2e})")

    if best['C2'] > 1e10:
        Tv_req = calculate_required_Tv(TARGET_C2, best['C2'], 1500)
        if Tv_req and 800 < Tv_req < 4500:
            print(f"\n>>> VIBRATIONAL TEMPERATURE CORRECTION <<<")
            print(f"Required Tv = {Tv_req:.0f} K")
            print(f"\nIf you re-analyze your C2 LIF data with Tv={Tv_req:.0f}K")
            print(f"instead of 1500K, all three species should match!")

if __name__ == '__main__':
    main()
