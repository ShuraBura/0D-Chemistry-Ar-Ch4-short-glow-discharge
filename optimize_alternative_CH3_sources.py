#!/usr/bin/env python3
"""
Optimizer focusing on ALTERNATIVE CH3 sources (not Ar*).

Strategy:
1. Keep Ar* at realistic levels (~1e8)
2. Boost e + CH4 → CH3 + H (electron-impact dissociation)
3. Boost ArPlus + CH4 → ArHPlus + CH3 (ion-molecule)
4. Boost CH + CH4 → CH2 + CH3
5. Keep ionization fixed at literature values
6. This should increase CH3 → C2H2 WITHOUT high Ar*!
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import os

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

os.makedirs('optimization_results_alternative_CH3', exist_ok=True)

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

    # Focus on ALTERNATIVE CH3 sources (not Ar*!)
    critical_reactions = [
        # === ALTERNATIVE CH3 SOURCES (KEY!) ===
        'e_CH4_CH3_H_cm3_1_1',           # e + CH4 → CH3 + H (12% currently!)
        'ArPlus_CH4_ArHPlus_CH3_cm3_5_9', # ArPlus + CH4 → ArHPlus + CH3 (25.5%!)
        'CH_CH4_CH2_CH3_cm3_7_39',       # CH + CH4 → CH2 + CH3 (3%)
        'H_CH4_CH3_H2_cm3_7_25',         # H + CH4 → CH3 + H2 (0.3%)

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
        'CH3_CH3_C2H2_H2_H2_cm3_7_49',   # 2CH3 → C2H2 (98%!)
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

    # NOTE: Ionization rates AND Ar* rates are FIXED!
    # We accept Ar* ~ 1e8 as realistic and focus on other CH3 sources

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
        CH3_final = y_final[species.index('CH3')]
        ArStar_final = y_final[species.index('ArStar')]

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

            filename = f'optimization_results_alternative_CH3/best_f{f:.1f}.json'
            with open(filename, 'w') as f_out:
                json.dump(result, f_out, indent=2)

            print(f"\nEval {eval_count}: f(x) = {f:.2f}")
            print(f"  H:    {H_final:.2e} ({H_final/target_densities['H']*100:5.1f}%)")
            print(f"  CH:   {CH_final:.2e} ({CH_final/target_densities['CH']*100:5.1f}%)")
            print(f"  C2:   {C2_final:.2e} ({C2_final/target_densities['C2']*100:5.1f}%)")
            print(f"  C2H2: {C2H2_final:.2e}")
            print(f"  CH3:  {CH3_final:.2e}")
            print(f"  Ar*:  {ArStar_final:.2e} (should stay ~1e8!)")
            print(f"  Ni/Ne: {Ni_over_Ne:.2f}")

        return f

    except Exception as e:
        print(f"Error in eval {eval_count}: {e}")
        return 1e6

if __name__ == '__main__':
    print("=" * 80)
    print("OPTIMIZER: BOOST ALTERNATIVE CH3 SOURCES (NOT Ar*!)")
    print("=" * 80)
    print()
    print("Strategy: Accept Ar* ~ 1e8 as realistic.")
    print("Instead, boost:")
    print("  1. e + CH4 → CH3 + H (electron-impact, currently 12%)")
    print("  2. ArPlus + CH4 → ArHPlus + CH3 (ion-molecule, currently 25.5%)")
    print("  3. CH + CH4 → CH2 + CH3 (neutral-neutral, currently 3%)")
    print()
    print("This increases CH3 → C2H2 without requiring high Ar*!")
    print()

    bounds = [
        (0.5, 3.0),      # Te
        (8.0, 10.5),     # log10(ne)
        (10.0, 300.0),   # E_field
    ]

    rate_bounds = {
        # === ALTERNATIVE CH3 SOURCES (BOOST THESE!) ===
        'e_CH4_CH3_H_cm3_1_1': (1e-10, 1e-08),           # e + CH4 → CH3 (wide range!)
        'ArPlus_CH4_ArHPlus_CH3_cm3_5_9': (5e-10, 5e-09), # ArPlus + CH4 → CH3
        'CH_CH4_CH2_CH3_cm3_7_39': (1e-12, 1e-10),       # CH + CH4 → CH3
        'H_CH4_CH3_H2_cm3_7_25': (1e-19, 1e-16),         # H + CH4 → CH3

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
        'CH3_CH3_C2H2_H2_H2_cm3_7_49': (1e-12, 5e-11),
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
    print("KEY: e + CH4 → CH3 can vary 100× (1e-10 to 1e-08 cm³/s)")
    print("     ArPlus + CH4 → CH3 can vary 10× (5e-10 to 5e-09 cm³/s)")
    print()

    result = differential_evolution(
        objective,
        bounds,
        maxiter=100,
        popsize=20,
        workers=1,
        updating='deferred',
        polish=False,
        seed=48
    )

    print("\n" + "=" * 80)
    print("OPTIMIZATION COMPLETE")
    print("=" * 80)
    print(f"Best f(x) = {result.fun:.2f}")
    print(f"Total evaluations: {eval_count}")
