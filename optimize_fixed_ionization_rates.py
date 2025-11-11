#!/usr/bin/env python3
"""
Optimizer with FIXED ionization rates at literature values.
Only tune C2/C2H2 chemistry, not fundamental plasma processes.

Root cause identified: Previous optimizers boosted ionization rates 10×
to achieve high C2H2, creating unphysical charge imbalance (Ni/Ne=215).

This optimizer uses CORRECT PHYSICS: ionization rates are fixed, we accept
whatever C2H2 level is achievable with proper charge balance.
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import os

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

os.makedirs('optimization_results_fixed_ionization', exist_ok=True)

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

    # ONLY TUNE C2/C2H2 CHEMISTRY - NOT ionization!
    critical_reactions = [
        # === C2 CHEMISTRY (tunable) ===
        'CH_CH_C2_H2_cm3_5_4',           # CH + CH → C2 + H2
        'e_C2H2_C2_H2_cm3_1_16',         # e + C2H2 → C2 + H2
        'C2HPlus_e_C2_H_cm3_6_18',       # C2H+ + e → C2 + H
        'C_CH_C2_H_cm3_7_4',             # C + CH → C2 + H
        'CH_C_C2_H_cm3_7_9',             # CH + C → C2 + H
        'C2H_H_C2_H2_cm3_7_47',          # C2H + H → C2 + H2
        'C2H2_H_C2_H2_H_cm3_7_50',       # H + C2H2 → C2 + H2 + H (MAIN!)
        'C2_H_CH_C_cm3_7_6',             # H + C2 → CH + C (DESTRUCTION!)

        # === C2 LOSSES (tunable) ===
        'stick_C2_9_9',
        'loss_C2_11_3',

        # === C2H2 CHEMISTRY (tunable) ===
        'stick_C2H2_9_11',               # C2H2 wall loss
        'loss_C2H2_11_19',               # C2H2 volumetric loss
        'CH3_CH3_C2H2_H2_H2_cm3_7_49',   # 2CH3 → C2H2 (67% of production!)
        'C_CH3_C2_H2_H_cm3_7_8',         # C + CH3 → C2H2
        'CH3_CH_C2H2_H2_cm3_7_16',       # CH + CH3 → C2H2

        # === CH/CH3 CONTROL (tunable) ===
        'loss_CH_11_9',
        'stick_CH_9_3',
        'stick_CH3_9_2',
        'loss_CH3_11_2',

        # === ELECTRON-ION RECOMBINATION (tunable with narrow bounds) ===
        'e_CH4Plus_CH3_H_cm3_6_4',
        'e_CH3Plus_CH2_H_cm3_6_6',
        'e_ArPlus_Ar_cm3_6_1',

        # === Ar* CONTROL (affects Penning ionization indirectly) ===
        'stick_ArStar_9_5',              # Ar* wall loss (tunable)
    ]

    # NOTE: Ionization rates (e + Ar, e + CH4, Ar* + CH4) are NOT tunable!
    # They remain at default literature values.

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

            filename = f'optimization_results_fixed_ionization/best_f{f:.1f}.json'
            with open(filename, 'w') as f_out:
                json.dump(result, f_out, indent=2)

            print(f"\nEval {eval_count}: f(x) = {f:.2f}")
            print(f"  H:    {H_final:.2e} ({H_final/target_densities['H']*100:5.1f}%)")
            print(f"  CH:   {CH_final:.2e} ({CH_final/target_densities['CH']*100:5.1f}%)")
            print(f"  C2:   {C2_final:.2e} ({C2_final/target_densities['C2']*100:5.1f}%)")
            print(f"  C2H2: {C2H2_final:.2e}")
            print(f"  Ni/Ne: {Ni_over_Ne:.2f}")

        return f

    except Exception as e:
        print(f"Error in eval {eval_count}: {e}")
        return 1e6

if __name__ == '__main__':
    print("=" * 80)
    print("OPTIMIZER WITH FIXED IONIZATION RATES")
    print("=" * 80)
    print()
    print("KEY CHANGE: Ionization rates (e+Ar, e+CH4, Ar*+CH4) are FIXED")
    print("at literature values. Only C2/C2H2 chemistry is tunable.")
    print()
    print("This prevents unphysical high-ionization states that were")
    print("achieving high C2H2 but violating charge balance (Ni/Ne >> 7).")
    print()
    print("We accept whatever C2H2 level is achievable with CORRECT PHYSICS.")
    print()

    bounds = [
        (0.5, 3.0),      # Te
        (8.0, 10.5),     # log10(ne)
        (10.0, 300.0),   # E_field
    ]

    rate_bounds = {
        # C2 chemistry
        'CH_CH_C2_H2_cm3_5_4': (5e-11, 2e-10),
        'e_C2H2_C2_H2_cm3_1_16': (3e-11, 1e-10),
        'C2HPlus_e_C2_H_cm3_6_18': (1e-07, 3e-07),
        'C_CH_C2_H_cm3_7_4': (5e-11, 2e-10),
        'CH_C_C2_H_cm3_7_9': (5e-11, 2e-10),
        'C2H_H_C2_H2_cm3_7_47': (5e-11, 2e-10),
        'C2H2_H_C2_H2_H_cm3_7_50': (5e-12, 2e-11),
        'C2_H_CH_C_cm3_7_6': (1e-12, 5e-11),          # Allow reduction!

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

        # Recombination (narrow bounds - these are well-known!)
        'e_CH4Plus_CH3_H_cm3_6_4': (5e-07, 1.5e-06),  # Narrow!
        'e_CH3Plus_CH2_H_cm3_6_6': (5e-07, 1.5e-06),  # Narrow!
        'e_ArPlus_Ar_cm3_6_1': (1e-07, 5e-07),        # Narrow!

        # Ar* control
        'stick_ArStar_9_5': (100, 2000),
    }

    critical_reactions = list(rate_bounds.keys())
    bounds.extend(rate_bounds.values())

    print(f"Optimizing {len(bounds)} parameters...")
    print(f"  - 3 plasma parameters")
    print(f"  - {len(critical_reactions)} reaction rates (NO ionization tuning!)")
    print()

    result = differential_evolution(
        objective,
        bounds,
        maxiter=100,
        popsize=20,
        workers=1,
        updating='deferred',
        polish=False,
        seed=46
    )

    print("\n" + "=" * 80)
    print("OPTIMIZATION COMPLETE")
    print("=" * 80)
    print(f"Best f(x) = {result.fun:.2f}")
    print(f"Total evaluations: {eval_count}")
