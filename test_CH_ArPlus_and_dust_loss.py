"""
Test impact of adding:
1. CH + Ar⁺ → CHPlus + Ar (ion-neutral reaction)
2. Dust/nanoparticle loss mechanisms

Can these reduce CH from 6355% to reasonable levels?
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

def run_test(ne_value, enable_dust=False, dust_params=None, description=""):
    """Run simulation with new CH loss mechanisms"""
    print(f"\n{'='*80}")
    print(f"TEST: {description}")
    print(f"{'='*80}")
    print(f"ne: {ne_value:.2e} cm⁻³")
    print(f"Dust loss: {'ENABLED' if enable_dust else 'DISABLED'}")
    if enable_dust and dust_params:
        print(f"  n_dust = {dust_params.get('dust_density', 1e8):.1e} cm⁻³")
        print(f"  r_dust = {dust_params.get('dust_radius', 50e-7)*1e7:.0f} nm")
        print(f"  α_dust = {dust_params.get('dust_sticking', 0.5):.2f}")

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
        'enable_dust_loss': enable_dust,
    }

    # Add dust parameters if specified
    if enable_dust and dust_params:
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

    print(f"  Number of reactions: {len(params['R'])}")
    dust_reactions = [tag for tag in params['tags'] if 'dust_loss' in tag]
    print(f"  Dust reactions: {len(dust_reactions)}")
    if dust_reactions:
        ch_dust_idx = params['tags'].index('dust_loss_CH_12') if 'dust_loss_CH_12' in params['tags'] else -1
        if ch_dust_idx >= 0:
            print(f"    dust_loss_CH_12 rate: {params['R'][ch_dust_idx].rate:.2e} s⁻¹")

    species = params['species']
    y0 = np.ones(len(species)) * 1e3
    y0[species.index('Ar')] = n_total * 0.85
    y0[species.index('CH4')] = n_total * 0.15
    y0[species.index('e')] = ne_value

    try:
        ode_func = PlasmaODE_Optimized(params)

        # Debug: Check if dust reactions are in rate_constants
        if enable_dust:
            dust_idx = params['tags'].index('dust_loss_CH_12') if 'dust_loss_CH_12' in params['tags'] else -1
            if dust_idx >= 0:
                print(f"  DEBUG: dust_loss_CH_12 in rate_constants[{dust_idx}] = {ode_func.rate_constants[dust_idx]:.2e}")

        sol = solve_ivp(ode_func, (0, 500), y0, method='BDF', rtol=1e-7, atol=1e-9, max_step=1.0)

        if not sol.success:
            print("✗ Integration failed!")
            return None

        y_final = sol.y[:, -1]

        # Get final densities
        H = y_final[species.index('H')]
        CH = y_final[species.index('CH')]
        C2 = y_final[species.index('C2')]
        ArPlus = y_final[species.index('ArPlus')]

        # Calculate Ni/Ne ratio
        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus', 'CHPlus', 'H2Plus']
        Ni_total = sum(y_final[species.index(ion)] for ion in ions if ion in species)
        Ne = y_final[species.index('e')]

        print(f"\nResults:")
        print(f"  H:    {H:.2e} ({H/targets['H']*100:6.2f}%)")
        print(f"  CH:   {CH:.2e} ({CH/targets['CH']*100:6.2f}%)")
        print(f"  C2:   {C2:.2e} ({C2/targets['C2']*100:6.2f}%)")
        print(f"  Ar⁺:  {ArPlus:.2e} cm⁻³")
        print(f"  Ni/Ne: {Ni_total/Ne:.2f}")

        # Calculate CH + Ar⁺ loss rate if reaction exists
        if 'ArPlus_CH_CHPlus_Ar_cm3_5_16' in k:
            rate_ArPlus_CH = k['ArPlus_CH_CHPlus_Ar_cm3_5_16'] * CH * ArPlus
            print(f"\nCH + Ar⁺ loss rate: {rate_ArPlus_CH:.2e} cm⁻³/s")

        # Calculate dust loss rate if enabled
        if enable_dust and 'dust_loss_CH_12' in k:
            rate_dust_CH = k['dust_loss_CH_12'] * CH
            print(f"CH dust loss rate:  {rate_dust_CH:.2e} cm⁻³/s")
            print(f"  (k_dust = {k['dust_loss_CH_12']:.2e} s⁻¹)")

        if Ni_total / Ne > 10:
            print("  ⚠️  RUNAWAY (Ni/Ne > 10)")
            return None
        else:
            print("  ✓ STABLE")

        return {
            'H': H,
            'CH': CH,
            'C2': C2,
            'Ni_Ne': Ni_total / Ne,
            'ArPlus': ArPlus,
        }

    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == '__main__':
    print("="*80)
    print("TEST CH LOSS MECHANISMS: CH + Ar⁺ AND DUST")
    print("="*80)

    ne_user = 2.3e9

    print("\n" + "="*80)
    print("BASELINE (no new mechanisms)")
    print("="*80)
    baseline_result = run_test(ne_user, enable_dust=False, description="Baseline (no CH+Ar⁺, no dust)")

    print("\n" + "="*80)
    print("TEST 1: Add CH + Ar⁺ only")
    print("="*80)
    result1 = run_test(ne_user, enable_dust=False, description="With CH + Ar⁺ → CHPlus + Ar")

    print("\n" + "="*80)
    print("TEST 2: Add moderate dust only (no CH+Ar⁺)")
    print("="*80)
    dust_moderate = {
        'dust_density': 1e8,
        'dust_radius': 50e-7,  # 50 nm
        'dust_sticking': 0.5,
    }
    result2 = run_test(ne_user, enable_dust=True, dust_params=dust_moderate,
                      description="Moderate dust (n=1e8, r=50nm, α=0.5)")

    print("\n" + "="*80)
    print("TEST 3: Add high dust only")
    print("="*80)
    dust_high = {
        'dust_density': 1e10,
        'dust_radius': 100e-7,  # 100 nm
        'dust_sticking': 0.8,
    }
    result3 = run_test(ne_user, enable_dust=True, dust_params=dust_high,
                      description="High dust (n=1e10, r=100nm, α=0.8)")

    print("\n" + "="*80)
    print("TEST 4: Add BOTH CH+Ar⁺ AND moderate dust")
    print("="*80)
    result4 = run_test(ne_user, enable_dust=True, dust_params=dust_moderate,
                      description="Both CH+Ar⁺ AND moderate dust")

    print("\n" + "="*80)
    print("TEST 5: Optimized dust scenario (target CH reduction)")
    print("="*80)
    # Calculate dust parameters needed to bring CH to ~100% of target
    # Current CH: 6355%, need to reduce by factor of ~60
    # Current CH loss: ~6e3 s⁻¹, need additional ~4.9e4 s⁻¹
    # k_dust = n_dust × π × r² × v_th × α
    # For n=5e9, r=70nm, α=0.6: k ≈ 4.9e4 s⁻¹
    dust_optimized = {
        'dust_density': 5e9,
        'dust_radius': 70e-7,  # 70 nm
        'dust_sticking': 0.6,
    }
    result5 = run_test(ne_user, enable_dust=True, dust_params=dust_optimized,
                      description="Optimized dust (n=5e9, r=70nm, α=0.6)")

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    results = [
        ("Baseline", baseline_result),
        ("CH+Ar⁺ only", result1),
        ("Moderate dust", result2),
        ("High dust", result3),
        ("CH+Ar⁺ + dust", result4),
        ("Optimized dust", result5),
    ]

    print(f"\n{'Configuration':<20} {'CH (cm⁻³)':<15} {'CH %':<10} {'C2 %':<10} {'H %':<10} {'Status'}")
    print("-"*80)
    for name, res in results:
        if res:
            ch_pct = res['CH'] / targets['CH'] * 100
            c2_pct = res['C2'] / targets['C2'] * 100
            h_pct = res['H'] / targets['H'] * 100
            status = "✓" if res['Ni_Ne'] < 10 else "✗"
            print(f"{name:<20} {res['CH']:<15.2e} {ch_pct:<10.1f} {c2_pct:<10.1f} {h_pct:<10.1f} {status}")
        else:
            print(f"{name:<20} {'FAILED':<15} {'':<10} {'':<10} {'':<10} ✗")

    print("\n" + "="*80)
    print("CONCLUSION:")
    print("="*80)

    if result5 and result5['CH'] / targets['CH'] < 2.0:
        print(f"\n✓ SUCCESS! Dust/nanoparticle loss CAN reduce CH to near-target levels!")
        print(f"  Optimized dust scenario gives CH = {result5['CH']/targets['CH']*100:.1f}% of target")
        print(f"  This suggests nanoparticle formation is CRITICAL in real plasmas")
    elif result3:
        print(f"\n⚠️  High dust scenario gives CH = {result3['CH']/targets['CH']*100:.1f}% of target")
        print(f"  Dust helps significantly but may not fully explain the discrepancy")
    else:
        print("\n✗ Dust loss alone insufficient to reduce CH to target levels")
