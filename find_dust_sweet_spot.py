"""
Find the "sweet spot" for dust parameters

Goal: Maximize C2 (keep it near 100%) while reducing CH as much as possible
Accept that we may not get CH to 100%, but find the best balance
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

def run_simulation(dust_params, ne_value=2.3e9):
    """Run simulation with given dust parameters"""

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
        'enable_dust_loss': True,
    }

    params.update(dust_params)

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
    print("FIND DUST 'SWEET SPOT' - BALANCE CH REDUCTION WITH C2 PRESERVATION")
    print("="*80)

    print("\n" + "="*80)
    print("STRATEGY 1: Optimize for C2 ≈ 100%, minimize CH")
    print("="*80 + "\n")

    print("Fine sweep around n_dust = 5.46e8 cm⁻³ (where C2 was 99.6%)")
    print("-"*80)

    # Focus around the region where C2 ≈ 100%
    dust_densities = np.linspace(3e8, 8e8, 30)
    dust_radius = 50e-7
    dust_sticking = 0.5

    print(f"{'n_dust (cm⁻³)':<15} {'CH %':<10} {'C2 %':<10} {'H %':<10} {'CH/C2 balance':<15}")
    print("-"*80)

    best_c2_config = None
    best_c2_error = 1e10

    for n_dust in dust_densities:
        dust_params = {
            'dust_density': n_dust,
            'dust_radius': dust_radius,
            'dust_sticking': dust_sticking,
        }

        result = run_simulation(dust_params)

        if result is None:
            continue

        # Prioritize C2 being close to 100%
        c2_error = abs(result['C2_pct'] - 100)

        marker = ""
        if c2_error < 5 and result['C2_pct'] > 95:  # C2 within 5% of target
            marker = " ← C2 ≈ 100%!"

        print(f"{n_dust:<15.2e} {result['CH_pct']:<10.1f} {result['C2_pct']:<10.2f} {result['H_pct']:<10.1f} {result['CH_pct']/result['C2_pct']:<15.1f}{marker}")

        if c2_error < best_c2_error:
            best_c2_error = c2_error
            best_c2_config = (dust_params.copy(), result.copy())

    print("\n" + "="*80)
    print("BEST CONFIGURATION FOR C2 ≈ 100%:")
    print("="*80)

    if best_c2_config:
        params, result = best_c2_config
        print(f"\nDust parameters:")
        print(f"  n_dust = {params['dust_density']:.3e} cm⁻³")
        print(f"  r_dust = {params['dust_radius']*1e7:.0f} nm")
        print(f"  α_dust = {params['dust_sticking']:.2f}")

        print(f"\nResults:")
        print(f"  C2: {result['C2']:.2e} cm⁻³ ({result['C2_pct']:6.2f}%) ← PRIMARY GOAL")
        print(f"  CH: {result['CH']:.2e} cm⁻³ ({result['CH_pct']:6.1f}%)")
        print(f"  H:  {result['H']:.2e} cm⁻³ ({result['H_pct']:6.1f}%)")
        print(f"  Ni/Ne: {result['Ni_Ne']:.2f}")

        print(f"\nImprovements vs baseline:")
        print(f"  CH: 6810% → {result['CH_pct']:.1f}% ({6810/result['CH_pct']:.1f}× reduction)")
        print(f"  C2:  107% → {result['C2_pct']:.1f}% ({result['C2_pct']/107*100:.0f}% of baseline)")

    print("\n" + "="*80)
    print("STRATEGY 2: Vary radius and sticking at optimal density")
    print("="*80 + "\n")

    if best_c2_config:
        best_n = best_c2_config[0]['dust_density']
        print(f"Using n_dust = {best_n:.2e} cm⁻³")
        print("-"*80)

        # Try different combinations
        radii = [30, 40, 50, 60, 70]
        stickings = [0.3, 0.4, 0.5, 0.6, 0.7]

        print(f"{'r (nm)':<8} {'α':<6} {'CH %':<10} {'C2 %':<10} {'H %':<10} {'Note'}")
        print("-"*80)

        best_balance = None
        best_balance_score = 1e10

        for r in radii:
            for alpha in stickings:
                dust_params = {
                    'dust_density': best_n,
                    'dust_radius': r * 1e-7,
                    'dust_sticking': alpha,
                }

                result = run_simulation(dust_params)

                if result is None:
                    continue

                # Balance score: prefer C2 > 90%, minimize CH
                if result['C2_pct'] > 90:
                    balance_score = result['CH_pct'] + 2 * abs(result['C2_pct'] - 100)
                else:
                    balance_score = 1e10  # Penalize C2 < 90%

                note = ""
                if result['C2_pct'] > 95 and result['C2_pct'] < 105:
                    note = " C2 optimal"
                if result['CH_pct'] < 2000:
                    note += " CH<2000%"

                print(f"{r:<8} {alpha:<6.1f} {result['CH_pct']:<10.1f} {result['C2_pct']:<10.2f} {result['H_pct']:<10.1f} {note}")

                if balance_score < best_balance_score:
                    best_balance_score = balance_score
                    best_balance = (dust_params.copy(), result.copy())

    print("\n" + "="*80)
    print("RECOMMENDED 'SWEET SPOT' CONFIGURATION:")
    print("="*80)

    if best_balance:
        params, result = best_balance
        print(f"\nOptimal dust parameters (best CH reduction while keeping C2 ≈ 100%):")
        print(f"  n_dust = {params['dust_density']:.3e} cm⁻³")
        print(f"  r_dust = {params['dust_radius']*1e7:.1f} nm")
        print(f"  α_dust = {params['dust_sticking']:.2f}")

        print(f"\nFinal balanced results:")
        print(f"  C2: {result['C2']:.2e} cm⁻³ ({result['C2_pct']:6.2f}% of target) ✓")
        print(f"  CH: {result['CH']:.2e} cm⁻³ ({result['CH_pct']:6.1f}% of target)")
        print(f"  H:  {result['H']:.2e} cm⁻³ ({result['H_pct']:6.1f}% of target)")
        print(f"  Ni/Ne: {result['Ni_Ne']:.2f} ✓")

        print(f"\nKey metrics:")
        print(f"  CH reduction: 6810% → {result['CH_pct']:.0f}% ({6810/result['CH_pct']:.1f}× improvement)")
        print(f"  C2 status: {result['C2_pct']:.1f}% of target (baseline was 107%)")
        print(f"  H status: {result['H_pct']:.1f}% of target (baseline was 104%)")

        print(f"\n{'='*80}")
        print("CONCLUSION:")
        print(f"{'='*80}\n")

        if result['C2_pct'] > 90 and result['CH_pct'] < 2000:
            print("✓ SUCCESS! Found sweet spot that:")
            print(f"  • Keeps C2 near target ({result['C2_pct']:.1f}%)")
            print(f"  • Reduces CH by {6810/result['CH_pct']:.1f}× (to {result['CH_pct']:.0f}%)")
            print(f"  • Maintains stable plasma (Ni/Ne = {result['Ni_Ne']:.2f})")
            print(f"\nThis is a REALISTIC dust scenario for Ar/CH4 plasmas!")
            print(f"CH is still {result['CH_pct']/100:.1f}× target, but this may reflect")
            print("experimental reality with nanoparticle formation.")
        else:
            print("⚠️  Trade-off remains: reducing CH further will destroy C2")
