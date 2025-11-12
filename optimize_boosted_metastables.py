#!/usr/bin/env python3
"""
Optimizer with FIXED ionization rates but BOOSTED Ar* metastable production.

Key insight: Ar* is the bottleneck! Ar* + CH4 → CH3 drives C2H2 production.

Strategy:
1. Keep ionization rates FIXED (e + Ar → Ar+ + 2e, e + CH4 → ions)
2. BOOST Ar* excitation rate (e + Ar → Ar*) to literature maximum
3. REDUCE Ar* wall loss to allow accumulation
4. This creates more CH3 via Ar* + CH4 → CH3 + H
5. More CH3 → More C2H2 (∝ [CH3]²) WITHOUT breaking charge balance!
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import os

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

os.makedirs('optimization_results_boosted_metastables', exist_ok=True)

def pressure_to_density(pressure_mTorr, T_K=400):
    kB = 1.38064852e-23
    Torr_to_Pa = 133.322
    P_Pa = pressure_mTorr * 1e-3 * Torr_to_Pa
    n_m3 = P_Pa / (kB * T_K)
    return n_m3 * 1e-6

target_densities = {
    'H': 2.52e14,
    'CH': 1.0e9,
    'C2': 5.6e11,
}

eval_count = 0
best_f = np.inf

def objective(x):
    global eval_count, best_f
    eval_count += 1

    Te = x[0]
    ne = 10**x[1]
    E_field = x[2]

    pressure_mTorr = 500.0
    n_total = pressure_to_density(pressure_mTorr)

    params = {
        'P': pressure_mTorr, 'Te': Te, 'ne': ne, 'E_field': E_field,
        'L_discharge': 0.45,
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
        },
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                    'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar',
                    'C2H4', 'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C',
                    'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H',
                    'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus',
                    'H2Plus', 'C2H2Star'],
    }

    k = define_rates(params)

    # Tunable reactions - NOW INCLUDING Ar* production/loss!
    critical_reactions = [
        # === Ar* PRODUCTION/LOSS (NEW: TUNABLE!) ===
        'e_Ar_ArStar_cm3_1_7',           # e + Ar → Ar* (BOOST THIS!) - FIXED NAME!
        'stick_ArStar_9_5',              # Ar* wall loss (REDUCE THIS!)

        # === Ar* PENNING REACTIONS (keep tunable) ===
        'ArStar_CH4_CH3_H_cm3_3_1',      # Ar* + CH4 → CH3 + H (KEY!)

        # === C2 CHEMISTRY ===
        'CH_CH_C2_H2_cm3_5_4',
        'e_C2H2_C2_H2_cm3_1_16',
        'C2HPlus_e_C2_H_cm3_6_18',
        'C_CH_C2_H_cm3_7_4',
        'CH_C_C2_H_cm3_7_9',
        'C2H_H_C2_H2_cm3_7_47',
        'C2H2_H_C2_H2_H_cm3_7_50',
        'C2_H_CH_C_cm3_7_6',

        # === C2 LOSSES ===
        'stick_C2_9_9',
        'loss_C2_11_3',

        # === C2H2 CHEMISTRY ===
        'stick_C2H2_9_11',
        'loss_C2H2_11_19',
        'CH3_CH3_C2H2_H2_H2_cm3_7_49',   # 2CH3 → C2H2 (CRITICAL!)
        'C_CH3_C2_H2_H_cm3_7_8',
        'CH3_CH_C2H2_H2_cm3_7_16',

        # === CH/CH3 CONTROL ===
        'loss_CH_11_9',
        'stick_CH_9_3',
        'stick_CH3_9_2',
        'loss_CH3_11_2',

        # === RECOMBINATION (narrow bounds) ===
        'e_CH4Plus_CH3_H_cm3_6_4',
        'e_CH3Plus_CH2_H_cm3_6_6',
        'e_ArPlus_Ar_cm3_6_1',
    ]

    # NOTE: Ionization rates (e + Ar → Ar+ + 2e, e + CH4 → ions) remain FIXED!
    # Only Ar* excitation is tunable now!

    for i, rate_name in enumerate(critical_reactions):
        if rate_name in k:
            k[rate_name] = x[3 + i]

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
            return 1e6

        y_final = sol.y[:, -1]

        H_final = y_final[species.index('H')]
        CH_final = y_final[species.index('CH')]
        C2_final = y_final[species.index('C2')]
        C2H2_final = y_final[species.index('C2H2')]
        ArStar_final = y_final[species.index('ArStar')]
        CH3_final = y_final[species.index('CH3')]

        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus', 'CHPlus', 'H2Plus']
        n_i_total = sum(y_final[species.index(ion)] for ion in ions)
        ne_final = y_final[species.index('e')]
        Ni_over_Ne = n_i_total / ne_final if ne_final > 0 else 0

        weights = {'H': 20.0, 'CH': 20.0, 'C2': 20.0}

        err_H = ((H_final - target_densities['H']) / target_densities['H'])**2
        err_CH = ((CH_final - target_densities['CH']) / target_densities['CH'])**2
        err_C2 = ((C2_final - target_densities['C2']) / target_densities['C2'])**2

        f = (weights['H'] * err_H + weights['CH'] * err_CH + weights['C2'] * err_C2)

        # Charge balance constraint
        if Ni_over_Ne < 2.0 or Ni_over_Ne > 7.0:
            charge_penalty = 200.0 * min((2.0 - Ni_over_Ne)**2, (Ni_over_Ne - 7.0)**2)
            f += charge_penalty

        if f < best_f:
            best_f = f
            result = {
                'pressure_mTorr': pressure_mTorr,
                'n_total': n_total,
                'Te': Te,
                'Ne': ne_final,
                'E_field': E_field,
                'n_i_total': n_i_total,
                'Ni_over_Ne': Ni_over_Ne,
                'charge_imbalance_pct': abs(ne_final - n_i_total) / ne_final * 100,
                'rate_values': {rate_name: k[rate_name] for rate_name in critical_reactions if rate_name in k},
                'target_densities': {'H': H_final, 'CH': CH_final, 'C2': C2_final},
                'all_densities': {sp: y_final[i] for i, sp in enumerate(species)}
            }

            filename = f'optimization_results_boosted_metastables/best_f{f:.1f}.json'
            with open(filename, 'w') as f_out:
                json.dump(result, f_out, indent=2)

            print(f"\nEval {eval_count}: f(x) = {f:.2f}")
            print(f"  H:    {H_final:.2e} ({H_final/target_densities['H']*100:5.1f}%)")
            print(f"  CH:   {CH_final:.2e} ({CH_final/target_densities['CH']*100:5.1f}%)")
            print(f"  C2:   {C2_final:.2e} ({C2_final/target_densities['C2']*100:5.1f}%)")
            print(f"  C2H2: {C2H2_final:.2e}")
            print(f"  Ar*:  {ArStar_final:.2e} (metastable!)")
            print(f"  CH3:  {CH3_final:.2e}")
            print(f"  Ni/Ne: {Ni_over_Ne:.2f}")

        return f

    except Exception as e:
        print(f"Error in eval {eval_count}: {e}")
        return 1e6

if __name__ == '__main__':
    print("=" * 80)
    print("OPTIMIZER: BOOST Ar* METASTABLES TO LITERATURE MAX")
    print("=" * 80)
    print()
    print("Strategy: Ar* is the bottleneck!")
    print("  - e + Ar → Ar* (excitation): BOOST to literature max")
    print("  - Ar* wall loss: REDUCE to allow accumulation")
    print("  - Ar* + CH4 → CH3: Tune for optimal CH3 production")
    print()
    print("This increases CH3 → increases C2H2 (∝ [CH3]²)")
    print("WITHOUT creating ions (so charge balance stays good!)")
    print()

    bounds = [
        (0.5, 3.0),      # Te
        (8.0, 10.5),     # log10(ne)
        (10.0, 300.0),   # E_field
    ]

    rate_bounds = {
        # === Ar* PRODUCTION/LOSS (CRITICAL!) ===
        'e_Ar_ArStar_cm3_1_7': (1e-09, 1e-07),        # BOOST! (literature max ~1e-7) - FIXED NAME!
        'stick_ArStar_9_5': (10, 500),                # REDUCE wall loss!

        # === Ar* PENNING ===
        'ArStar_CH4_CH3_H_cm3_3_1': (1e-11, 5e-10),   # Ar* + CH4 → CH3

        # C2 chemistry
        'CH_CH_C2_H2_cm3_5_4': (5e-11, 2e-10),
        'e_C2H2_C2_H2_cm3_1_16': (3e-11, 1e-10),
        'C2HPlus_e_C2_H_cm3_6_18': (1e-07, 3e-07),
        'C_CH_C2_H_cm3_7_4': (5e-11, 2e-10),
        'CH_C_C2_H_cm3_7_9': (5e-11, 2e-10),
        'C2H_H_C2_H2_cm3_7_47': (5e-11, 2e-10),
        'C2H2_H_C2_H2_H_cm3_7_50': (5e-12, 2e-11),
        'C2_H_CH_C_cm3_7_6': (1e-12, 5e-11),

        # C2 losses
        'stick_C2_9_9': (1000, 8000),
        'loss_C2_11_3': (500, 3000),

        # C2H2 chemistry
        'stick_C2H2_9_11': (100, 2000),
        'loss_C2H2_11_19': (100, 2000),
        'CH3_CH3_C2H2_H2_H2_cm3_7_49': (1e-12, 5e-11),  # 2CH3 → C2H2
        'C_CH3_C2_H2_H_cm3_7_8': (5e-11, 3e-10),
        'CH3_CH_C2H2_H2_cm3_7_16': (5e-11, 3e-10),

        # CH/CH3 control
        'loss_CH_11_9': (500, 5000),
        'stick_CH_9_3': (500, 5000),
        'stick_CH3_9_2': (100, 2000),
        'loss_CH3_11_2': (100, 2000),

        # Recombination (narrow)
        'e_CH4Plus_CH3_H_cm3_6_4': (5e-07, 1.5e-06),
        'e_CH3Plus_CH2_H_cm3_6_6': (5e-07, 1.5e-06),
        'e_ArPlus_Ar_cm3_6_1': (1e-07, 5e-07),
    }

    critical_reactions = list(rate_bounds.keys())
    bounds.extend(rate_bounds.values())

    print(f"Optimizing {len(bounds)} parameters...")
    print(f"  - 3 plasma parameters")
    print(f"  - {len(critical_reactions)} reaction rates")
    print()
    print("KEY: e + Ar → Ar* can now reach 1e-7 cm³/s (literature max)")
    print("     Ar* wall loss can be as low as 10 s⁻¹")
    print()

    result = differential_evolution(
        objective,
        bounds,
        maxiter=100,
        popsize=20,
        workers=1,
        updating='deferred',
        polish=False,
        seed=47
    )

    print("\n" + "=" * 80)
    print("OPTIMIZATION COMPLETE")
    print("=" * 80)
    print(f"Best f(x) = {result.fun:.2f}")
    print(f"Total evaluations: {eval_count}")
