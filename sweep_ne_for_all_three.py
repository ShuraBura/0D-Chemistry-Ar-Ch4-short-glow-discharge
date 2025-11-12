"""
Test if adjusting ne (electron density) can help match all three species

Current ne = 2.3e9 cm⁻³ was set to achieve C2 target
But maybe a different ne + dust combination works better for ALL THREE?
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

def run_simulation(ne_value, dust_density=0, c2h_suppression=1.0):
    """Run simulation with given ne, dust, and C2+H suppression"""

    P = 500.0
    n_total = pressure_to_density(P)

    params = {
        'P': P,
        'Te': baseline['Te'],
        'ne': ne_value,
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

    if dust_density > 0:
        params['enable_dust_loss'] = True
        params['dust_density'] = dust_density
        params['dust_radius'] = 50e-7
        params['dust_sticking'] = 0.5

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
        'C2_H_CH_C_cm3_7_6': c2h_suppression,
    }

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

    except Exception:
        return None


if __name__ == '__main__':
    print("="*80)
    print("SWEEP ne TO FIND CONFIGURATION MATCHING ALL THREE SPECIES")
    print("="*80)

    print("\n" + "="*80)
    print("STRATEGY: Sweep ne with different dust levels")
    print("="*80 + "\n")

    ne_values = np.logspace(8, 10, 15)  # 1e8 to 1e10 cm⁻³
    dust_levels = [0, 1e8, 3e8, 5e8]

    best_overall_score = 1e10
    best_overall_config = None

    for dust in dust_levels:
        print(f"\n{'='*80}")
        print(f"DUST: {dust:.1e} cm⁻³, C2+H suppression: 0.01")
        print(f"{'='*80}")
        print(f"{'ne (cm⁻³)':<15} {'H %':<10} {'CH %':<10} {'C2 %':<10} {'RMS err':<12} {'Status'}")
        print("-"*80)

        best_score = 1e10
        best_result = None
        best_ne = None

        for ne in ne_values:
            result = run_simulation(ne, dust_density=dust, c2h_suppression=0.01)

            if result is None:
                continue

            # RMS error
            score = np.sqrt(((result['H_pct'] - 100)**2 + (result['CH_pct'] - 100)**2 + (result['C2_pct'] - 100)**2) / 3)

            marker = ""
            if score < 100:
                marker = " ← Good!"
            if score < 50:
                marker = " ← Excellent!"
            if score < 30:
                marker = " ← BEST!"

            status = "✓" if result['Ni_Ne'] < 10 else "✗"

            print(f"{ne:<15.2e} {result['H_pct']:<10.1f} {result['CH_pct']:<10.1f} {result['C2_pct']:<10.1f} {score:<12.1f} {status}{marker}")

            if score < best_score:
                best_score = score
                best_result = result
                best_ne = ne

            if score < best_overall_score:
                best_overall_score = score
                best_overall_config = (ne, dust, 0.01, result)

        if best_result:
            print(f"\nBest for this dust level:")
            print(f"  ne = {best_ne:.2e} cm⁻³")
            print(f"  RMS error: {best_score:.1f}%")

    print("\n" + "="*80)
    print("FINAL BEST CONFIGURATION:")
    print("="*80)

    if best_overall_config:
        ne, dust, c2h_supp, result = best_overall_config

        print(f"\nOptimal parameters:")
        print(f"  ne = {ne:.3e} cm⁻³")
        print(f"  Dust density = {dust:.2e} cm⁻³")
        print(f"  Dust radius = 50 nm")
        print(f"  Dust sticking = 0.5")
        print(f"  C2 + H → CH + C suppression = {c2h_supp:.3f}")

        print(f"\nFinal results:")
        print(f"  H:  {result['H']:.2e} cm⁻³ ({result['H_pct']:6.2f}% of target)")
        print(f"  CH: {result['CH']:.2e} cm⁻³ ({result['CH_pct']:6.2f}% of target)")
        print(f"  C2: {result['C2']:.2e} cm⁻³ ({result['C2_pct']:6.2f}% of target)")
        print(f"  Ni/Ne: {result['Ni_Ne']:.2f}")

        print(f"\nRMS error from targets: {best_overall_score:.1f}%")

        print(f"\n{'='*80}")
        print("ASSESSMENT:")
        print(f"{'='*80}\n")

        if best_overall_score < 30:
            print("✓ SUCCESS! All three species within 30% RMS")
        elif best_overall_score < 50:
            print("✓ GOOD! All three species within 50% RMS")
        elif best_overall_score < 100:
            print("⚠️  MODERATE: All three within 100% RMS")
        else:
            print("✗ POOR: Cannot match all three simultaneously")
            print("\nThis suggests:")
            print("  • Fundamental chemistry issue")
            print("  • Missing physics (e.g., transport, spatial gradients)")
            print("  • Experimental targets may not be self-consistent")
