#!/usr/bin/env python3
"""
Find parameter configuration where CH and H match well (even if both low/high by similar factors).
Then use C2 discrepancy to back-calculate what Tv should be.
"""

import numpy as np
from pathlib import Path
import sys

# Add the directory containing the modules to the Python path
sys.path.insert(0, str(Path(__file__).parent))

from build_reactions import build_reactions
from odefun_optimized import create_optimized_odefun
from define_rates import define_rates
from scipy.integrate import solve_ivp

# Target densities (experimental, from sheath edge)
TARGET_H = 2.52e14   # cm^-3
TARGET_CH = 1.00e9   # cm^-3
TARGET_C2 = 5.60e11  # cm^-3

# C2 vibrational temperature correction
# If measured with Tv_assumed but actual Tv is different:
# True_C2 = Measured_C2 * [f_v1(Tv_assumed) / f_v1(Tv_actual)]
def c2_vibrational_fraction(Tv, v=1):
    """Calculate fraction of C2 in vibrational level v at temperature Tv"""
    omega_e = 1641  # cm^-1 for C2
    delta_E_K = 1.4388 * omega_e  # Convert to Kelvin (2361 K)

    # Boltzmann distribution
    Z = 0.0
    for vi in range(10):  # Sum over first 10 levels
        Z += np.exp(-vi * delta_E_K / Tv)

    return np.exp(-v * delta_E_K / Tv) / Z

def calculate_required_Tv(measured_C2, predicted_C2, Tv_assumed=1500):
    """
    Back-calculate what Tv would make measured C2 equal predicted C2.

    measured_C2: What was measured (5.60e11) assuming Tv_assumed
    predicted_C2: What the model predicts
    Tv_assumed: Temperature assumed in LIF analysis (1500K)

    Returns: Actual Tv that would reconcile measurement with model
    """
    # True_C2 = Measured_C2 * [f_v1(Tv_assumed) / f_v1(Tv_actual)]
    # We want: True_C2 = Predicted_C2
    # So: Predicted_C2 = Measured_C2 * [f_v1(Tv_assumed) / f_v1(Tv_actual)]
    # => f_v1(Tv_actual) = f_v1(Tv_assumed) * Measured_C2 / Predicted_C2

    f_v1_assumed = c2_vibrational_fraction(Tv_assumed, v=1)
    target_f_v1 = f_v1_assumed * measured_C2 / predicted_C2

    # Search for Tv that gives this f_v1
    for Tv_test in np.linspace(500, 5000, 1000):
        f_v1_test = c2_vibrational_fraction(Tv_test, v=1)
        if abs(f_v1_test - target_f_v1) / target_f_v1 < 0.01:  # Within 1%
            return Tv_test

    return None  # Not found

def run_simulation(ne, c2h_factor, dust_density, chh_factor, H_drift_gain=5.7e16):
    """Run simulation with given parameters"""

    # Parameters
    params = {
        'P': 500e-3,  # Torr
        'Tgas': 570,  # K - CORRECTED from 400K
        'Te_eV': 1.3,
        'ne': ne,
        'diffusion_length': 0.6,  # cm
        'wall_loss_fraction': 0.1,
        'electrode_area': 177,  # cm^2
        'AR': 1.0,
        'CH4': 0.03
    }

    # Build reactions
    sto, k_values, k_names, species, species_charges = build_reactions(params)

    # Apply rate modifications
    modified_ks = {}
    for i, name in enumerate(k_names):
        if name == 'C2_H_CH_C_cm3':
            modified_ks[name] = k_values[i] * c2h_factor
        elif name == 'CH_H_C_H2_cm3':
            modified_ks[name] = k_values[i] * chh_factor
        else:
            modified_ks[name] = k_values[i]

    # Create ODE function
    odefun = create_optimized_odefun(
        sto=sto,
        k=modified_ks,
        k_names=k_names,
        species=species,
        species_charges=species_charges,
        params=params,
        dust_density=dust_density,
        H_drift_gain=H_drift_gain
    )

    # Initial conditions
    n0 = np.zeros(len(species))
    n0[species.index('AR')] = params['P'] * 760 * params['AR'] / (1.38e-16 * params['Tgas'])
    n0[species.index('CH4')] = params['P'] * 760 * params['CH4'] / (1.38e-16 * params['Tgas'])
    n0[species.index('E')] = params['ne']
    n0[species.index('AR+')] = params['ne']

    # Integrate
    t_span = (0, 1e-3)  # 1 ms
    t_eval = np.logspace(-9, -3, 100)

    try:
        sol = solve_ivp(odefun, t_span, n0, method='BDF', t_eval=t_eval,
                       rtol=1e-6, atol=1e-3)

        if not sol.success:
            return None, None, None

        # Extract final densities
        n_final = sol.y[:, -1]

        H_pred = n_final[species.index('H')] if 'H' in species else 0
        CH_pred = n_final[species.index('CH')] if 'CH' in species else 0
        C2_pred = n_final[species.index('C2')] if 'C2' in species else 0

        return H_pred, CH_pred, C2_pred

    except Exception as e:
        return None, None, None

def main():
    print("="*80)
    print("FINDING CONFIGURATION WHERE CH AND H MATCH WELL")
    print("="*80)
    print()
    print(f"Target H:  {TARGET_H:.2e} cm^-3")
    print(f"Target CH: {TARGET_CH:.2e} cm^-3")
    print(f"Target C2: {TARGET_C2:.2e} cm^-3")
    print()
    print("Strategy: Find where |H_error| ≈ |CH_error| (both match by similar factor)")
    print()

    # Parameter sweep - focus on lower ne range where earlier results were better
    ne_values = [3e8, 4e8, 5e8, 6e8, 7e8, 8e8]
    c2h_factors = [0.3, 0.4, 0.5, 0.6]
    dust_values = [5e7, 1e8, 1.5e8, 2e8]
    chh_factors = [1.5, 2.0, 2.5, 3.0]

    results = []
    total_runs = len(ne_values) * len(c2h_factors) * len(dust_values) * len(chh_factors)
    run_count = 0

    print(f"Running {total_runs} simulations...")
    print()

    for ne in ne_values:
        for c2h_factor in c2h_factors:
            for dust in dust_values:
                for chh_factor in chh_factors:
                    run_count += 1
                    if run_count % 20 == 0:
                        print(f"Progress: {run_count}/{total_runs} ({100*run_count/total_runs:.1f}%)")

                    H_pred, CH_pred, C2_pred = run_simulation(
                        ne, c2h_factor, dust, chh_factor
                    )

                    if H_pred is None:
                        continue

                    # Calculate errors
                    H_ratio = H_pred / TARGET_H
                    CH_ratio = CH_pred / TARGET_CH
                    C2_ratio = C2_pred / TARGET_C2

                    # Metric: How well do H and CH match (by similar factors)?
                    # We want log(H_ratio) ≈ log(CH_ratio)
                    ch_h_match_score = abs(np.log10(H_ratio) - np.log10(CH_ratio))

                    # Also want them not too far from 1.0
                    average_match = np.sqrt(H_ratio * CH_ratio)

                    # Combined score: prefer when H and CH match each other,
                    # and both are reasonably close to target
                    combined_score = ch_h_match_score + abs(np.log10(average_match))

                    results.append({
                        'ne': ne,
                        'c2h_factor': c2h_factor,
                        'dust': dust,
                        'chh_factor': chh_factor,
                        'H_pred': H_pred,
                        'CH_pred': CH_pred,
                        'C2_pred': C2_pred,
                        'H_ratio': H_ratio,
                        'CH_ratio': CH_ratio,
                        'C2_ratio': C2_ratio,
                        'ch_h_match_score': ch_h_match_score,
                        'combined_score': combined_score
                    })

    print()
    print(f"Completed {len(results)} successful runs")
    print()

    # Sort by combined score (best match of CH and H)
    results.sort(key=lambda x: x['combined_score'])

    print("="*80)
    print("TOP 10 CONFIGURATIONS WHERE CH AND H MATCH BEST")
    print("="*80)
    print()

    for i, res in enumerate(results[:10]):
        print(f"\n--- Rank {i+1} ---")
        print(f"Parameters:")
        print(f"  ne = {res['ne']:.2e} cm^-3")
        print(f"  C2+H factor = {res['c2h_factor']:.2f}")
        print(f"  dust = {res['dust']:.2e} cm^-3")
        print(f"  CH+H factor = {res['chh_factor']:.2f}")
        print(f"\nPredicted densities:")
        print(f"  H:  {res['H_pred']:.2e} ({res['H_ratio']*100:5.1f}% of target)")
        print(f"  CH: {res['CH_pred']:.2e} ({res['CH_ratio']*100:5.1f}% of target)")
        print(f"  C2: {res['C2_pred']:.2e} ({res['C2_ratio']*100:5.1f}% of target)")
        print(f"\nMatch quality:")
        print(f"  CH-H match score: {res['ch_h_match_score']:.3f} (lower is better)")
        print(f"  Combined score: {res['combined_score']:.3f}")

        # Calculate required Tv for C2
        if res['C2_pred'] > 0:
            Tv_required = calculate_required_Tv(TARGET_C2, res['C2_pred'], Tv_assumed=1500)
            if Tv_required:
                print(f"\n>>> Tv correction <<<")
                print(f"  If C2 model predicts {res['C2_pred']:.2e}")
                print(f"  And LIF measures {TARGET_C2:.2e} assuming Tv=1500K")
                print(f"  Then ACTUAL Tv must be: {Tv_required:.0f} K")

                # Verify
                f_v1_1500 = c2_vibrational_fraction(1500, v=1)
                f_v1_actual = c2_vibrational_fraction(Tv_required, v=1)
                true_C2 = TARGET_C2 * (f_v1_1500 / f_v1_actual)
                print(f"  Verification: True C2 = {true_C2:.2e} (should match prediction)")

    print("\n" + "="*80)
    print("BEST CONFIGURATION SUMMARY")
    print("="*80)

    best = results[0]
    print(f"\nne = {best['ne']:.2e} cm^-3")
    print(f"C2+H suppression factor = {best['c2h_factor']:.2f}")
    print(f"dust_density = {best['dust']:.2e} cm^-3")
    print(f"CH+H multiplier = {best['chh_factor']:.2f}")
    print(f"\nH:  {best['H_pred']:.2e} ({best['H_ratio']*100:5.1f}% of target)")
    print(f"CH: {best['CH_pred']:.2e} ({best['CH_ratio']*100:5.1f}% of target)")
    print(f"C2: {best['C2_pred']:.2e} ({best['C2_ratio']*100:5.1f}% of target)")

    if best['C2_pred'] > 0:
        Tv_required = calculate_required_Tv(TARGET_C2, best['C2_pred'], Tv_assumed=1500)
        if Tv_required:
            print(f"\n>>> REQUIRED Tv = {Tv_required:.0f} K <<<")
            print(f"\nThis means:")
            print(f"  - If the chemistry model is correct")
            print(f"  - And C2 is truly at {best['C2_pred']:.2e} cm^-3")
            print(f"  - Then your LIF measurement analysis should use Tv={Tv_required:.0f}K")
            print(f"  - Not the currently assumed Tv=1500K")

if __name__ == '__main__':
    main()
