"""
CRITICAL: All three species (H, CH, C2) measured in SAME setup
Must match ALL THREE simultaneously!

Targets (all from same experiment):
- H:  2.52×10¹⁴ cm⁻³
- CH: 1.00×10⁹ cm⁻³
- C2: 5.60×10¹¹ cm⁻³

Current status:
- Baseline (no dust):    H=104%, CH=6810%, C2=107%
- Sweet spot dust:       H=11%,  CH=1575%, C2=100%

Neither works! Need new approach.

Question: What chemistry change allows ALL THREE to match simultaneously?
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

def run_simulation(params_override, verbose=False):
    """Run simulation with parameter overrides"""

    P = 500.0
    n_total = pressure_to_density(P)

    params = {
        'P': P,
        'Te': baseline['Te'],
        'ne': 2.3e9,
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

    params.update(params_override)

    k = define_rates(params)

    # Apply baseline rates
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    # Apply physically realistic multipliers
    rate_mults = {
        'e_CH4_CH3_HMinus_cm3_8_1': 10.0,
        'ArStar_CH4_CH3_H_cm3_3_1': 1.4,
        'e_CH4Plus_CH3_H_cm3_6_4': 1.5,
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
    }

    # Apply any rate overrides from params_override
    if 'rate_overrides' in params_override:
        rate_mults.update(params_override['rate_overrides'])

    for rate_name, mult in rate_mults.items():
        if rate_name in k:
            k[rate_name] *= mult

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    species = params['species']
    y0 = np.ones(len(species)) * 1e3
    y0[species.index('Ar')] = n_total * 0.85
    y0[species.index('CH4')] = n_total * 0.15
    y0[species.index('e')] = params['ne']

    try:
        ode_func = PlasmaODE_Optimized(params)
        sol = solve_ivp(ode_func, (0, 500), y0, method='BDF', rtol=1e-7, atol=1e-9, max_step=1.0)

        if not sol.success:
            return None

        y_final = sol.y[:, -1]

        H = y_final[species.index('H')]
        CH = y_final[species.index('CH')]
        C2 = y_final[species.index('C2')]

        # Calculate Ni/Ne ratio
        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus', 'CHPlus', 'H2Plus']
        Ni_total = sum(y_final[species.index(ion)] for ion in ions if ion in species)
        Ne = y_final[species.index('e')]

        if Ni_total / Ne > 10:
            return None  # Runaway

        return {
            'H': H,
            'CH': CH,
            'C2': C2,
            'Ni_Ne': Ni_total / Ne,
            'H_pct': H / targets['H'] * 100,
            'CH_pct': CH / targets['CH'] * 100,
            'C2_pct': C2 / targets['C2'] * 100,
        }

    except Exception as e:
        if verbose:
            print(f"Error: {e}")
        return None


if __name__ == '__main__':
    print("="*80)
    print("CRITICAL CONSTRAINT: MUST MATCH ALL THREE SPECIES SIMULTANEOUSLY")
    print("="*80)
    print("\nAll measured in same experimental setup:")
    print(f"  H:  {targets['H']:.2e} cm⁻³")
    print(f"  CH: {targets['CH']:.2e} cm⁻³")
    print(f"  C2: {targets['C2']:.2e} cm⁻³")

    print("\n" + "="*80)
    print("HYPOTHESIS TESTING")
    print("="*80)

    print("\n" + "="*80)
    print("HYPOTHESIS 1: Reduce C2 + H → CH + C rate constant")
    print("="*80)
    print("\nC2 + H → CH + C is creating massive CH (5.6×10¹⁵ cm⁻³/s)")
    print("Literature value: k = 3.53×10⁻¹¹ cm³/s")
    print("Testing: What if this rate is overestimated?")
    print("-"*80)

    # Sweep C2+H→CH+C suppression
    suppressions = [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]

    print(f"{'Suppression':<12} {'H %':<10} {'CH %':<10} {'C2 %':<10} {'Score':<12} {'Status'}")
    print("-"*80)

    best_score = 1e10
    best_config = None

    for supp in suppressions:
        params_override = {
            'rate_overrides': {
                'C2_H_CH_C_cm3_7_6': supp,  # Suppress C2 + H → CH + C
            }
        }

        result = run_simulation(params_override)

        if result is None:
            print(f"{supp:<12.3f} {'FAILED':<10} {'FAILED':<10} {'FAILED':<10} {'---':<12} ✗")
            continue

        # Score: sum of squared errors from 100%
        score = (result['H_pct'] - 100)**2 + (result['CH_pct'] - 100)**2 + (result['C2_pct'] - 100)**2
        score = np.sqrt(score / 3)  # RMS error

        status = "✓" if result['Ni_Ne'] < 10 else "✗"

        marker = ""
        if score < 50:  # All within 50% of target
            marker = " ← Good!"
        if score < 20:  # All within 20% of target
            marker = " ← EXCELLENT!"

        print(f"{supp:<12.3f} {result['H_pct']:<10.1f} {result['CH_pct']:<10.1f} {result['C2_pct']:<10.1f} {score:<12.1f} {status}{marker}")

        if score < best_score:
            best_score = score
            best_config = (supp, result)

    print("\n" + "="*80)
    print("HYPOTHESIS 2: Adjust dust + C2+H suppression")
    print("="*80)
    print("\nCombine moderate dust with C2+H rate reduction")
    print("-"*80)

    # Try combinations
    dust_levels = [0, 1e8, 3e8, 5e8]
    c2h_suppressions = [0.01, 0.02, 0.05, 0.1, 0.2]

    print(f"{'Dust (cm⁻³)':<12} {'C2+H supp':<12} {'H %':<10} {'CH %':<10} {'C2 %':<10} {'Score':<12}")
    print("-"*80)

    for dust in dust_levels:
        for supp in c2h_suppressions:
            if dust > 0:
                params_override = {
                    'enable_dust_loss': True,
                    'dust_density': dust,
                    'dust_radius': 50e-7,
                    'dust_sticking': 0.5,
                    'rate_overrides': {
                        'C2_H_CH_C_cm3_7_6': supp,
                    }
                }
            else:
                params_override = {
                    'rate_overrides': {
                        'C2_H_CH_C_cm3_7_6': supp,
                    }
                }

            result = run_simulation(params_override)

            if result is None:
                continue

            score = np.sqrt(((result['H_pct'] - 100)**2 + (result['CH_pct'] - 100)**2 + (result['C2_pct'] - 100)**2) / 3)

            marker = ""
            if score < 30:
                marker = " ← Good!"
            if score < 15:
                marker = " ← EXCELLENT!"
            if score < 10:
                marker = " ← PERFECT!"

            print(f"{dust:<12.1e} {supp:<12.3f} {result['H_pct']:<10.1f} {result['CH_pct']:<10.1f} {result['C2_pct']:<10.1f} {score:<12.1f}{marker}")

            if score < best_score:
                best_score = score
                best_config = (params_override, result)

    print("\n" + "="*80)
    print("BEST OVERALL CONFIGURATION:")
    print("="*80)

    if best_config:
        config, result = best_config

        print(f"\nConfiguration:")
        if isinstance(config, dict) and 'enable_dust_loss' in config:
            print(f"  Dust enabled: YES")
            print(f"  n_dust = {config.get('dust_density', 0):.2e} cm⁻³")
            print(f"  r_dust = {config.get('dust_radius', 0)*1e7:.0f} nm")
            print(f"  α_dust = {config.get('dust_sticking', 0):.2f}")
        else:
            print(f"  Dust enabled: NO")

        if 'rate_overrides' in config:
            for key, val in config['rate_overrides'].items():
                print(f"  {key}: ×{val:.3f}")

        print(f"\nFinal results:")
        print(f"  H:  {result['H']:.2e} cm⁻³ ({result['H_pct']:6.2f}% of target)")
        print(f"  CH: {result['CH']:.2e} cm⁻³ ({result['CH_pct']:6.2f}% of target)")
        print(f"  C2: {result['C2']:.2e} cm⁻³ ({result['C2_pct']:6.2f}% of target)")
        print(f"  Ni/Ne: {result['Ni_Ne']:.2f}")

        print(f"\nRMS error from targets: {best_score:.1f}%")

        print(f"\n{'='*80}")
        print("CONCLUSION:")
        print(f"{'='*80}\n")

        if best_score < 20:
            print("✓ SUCCESS! Found configuration matching ALL THREE species!")
            print(f"  All within {best_score:.0f}% RMS error of targets")
        elif best_score < 50:
            print("⚠️  PARTIAL SUCCESS: Reasonably close to all three targets")
            print(f"  RMS error: {best_score:.0f}%")
        else:
            print("✗ Could not find configuration matching all three simultaneously")
            print(f"  Best RMS error: {best_score:.0f}%")
            print("\nPossible issues:")
            print("  1. C2 + H → CH + C rate may need even more suppression")
            print("  2. Missing major CH loss mechanism beyond dust")
            print("  3. Plasma conditions (Te, ne) may be incorrect")
            print("  4. Other rate constants may need adjustment")
