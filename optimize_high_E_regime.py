#!/usr/bin/env python3
"""
User's new approach: Don't stress about Ni/Ne, focus on hitting species targets!

Constraints:
  - Ne ~ 2.3e9 cm⁻³ (fixed, very specific!)
  - Te ~ 1.5 eV
  - E = 50-300 V/cm
  - Target: H, CH, C2 densities
  - Monitor Ni/Ne but don't constrain it

Focus: Hit the species targets without worrying about charge balance.
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import os

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

os.makedirs('optimization_results_high_E_regime', exist_ok=True)

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
    ne = 2.3e9  # FIXED per user
    E_field = x[1]

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

    # Tunable reactions
    critical_reactions = [
        'CH_CH_C2_H2_cm3_5_4',
        'e_C2H2_C2_H2_cm3_1_16',
        'C2HPlus_e_C2_H_cm3_6_18',
        'C_CH_C2_H_cm3_7_4',
        'CH_C_C2_H_cm3_7_9',
        'C2H_H_C2_H2_cm3_7_47',
        'C2H2_H_C2_H2_H_cm3_7_50',
        'C2_H_CH_C_cm3_7_6',
        'stick_C2_9_9',
        'loss_C2_11_3',
        'stick_C2H2_9_11',
        'loss_C2H2_11_19',
        'CH3_CH3_C2H2_H2_H2_cm3_7_49',
        'C_CH3_C2_H2_H_cm3_7_8',
        'CH3_CH_C2H2_H2_cm3_7_16',
        'loss_CH_11_9',
        'stick_CH_9_3',
        'stick_CH3_9_2',
        'e_CH4Plus_CH3_H_cm3_6_4',
    ]

    for i, rate_name in enumerate(critical_reactions):
        if rate_name in k:
            k[rate_name] = x[2 + i]

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

        # ONLY optimize for species targets, NO charge balance penalty!
        weights = {'H': 20.0, 'CH': 20.0, 'C2': 20.0}

        err_H = ((H_final - target_densities['H']) / target_densities['H'])**2
        err_CH = ((CH_final - target_densities['CH']) / target_densities['CH'])**2
        err_C2 = ((C2_final - target_densities['C2']) / target_densities['C2'])**2

        f = (weights['H'] * err_H + weights['CH'] * err_CH + weights['C2'] * err_C2)

        # NO charge balance penalty per user: "don't stress about it"

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

            filename = f'optimization_results_high_E_regime/best_f{f:.1f}.json'
            with open(filename, 'w') as f_out:
                json.dump(result, f_out, indent=2)

            print(f"\nEval {eval_count}: f(x) = {f:.2f}")
            print(f"  Te:     {Te:.2f} eV")
            print(f"  Ne:     {ne_final:.2e} cm⁻³ (fixed at 2.3e9)")
            print(f"  E:      {E_field:.1f} V/cm")
            print(f"  H:      {H_final:.2e} ({H_final/target_densities['H']*100:5.1f}%)")
            print(f"  CH:     {CH_final:.2e} ({CH_final/target_densities['CH']*100:5.1f}%)")
            print(f"  C2:     {C2_final:.2e} ({C2_final/target_densities['C2']*100:5.1f}%)")
            print(f"  C2H2:   {C2H2_final:.2e}")
            print(f"  CH3:    {CH3_final:.2e}")
            print(f"  Ni/Ne:  {Ni_over_Ne:.2f} (monitoring only)")

        return f

    except Exception as e:
        print(f"Error in eval {eval_count}: {e}")
        return 1e6

if __name__ == '__main__':
    print("=" * 80)
    print("HIGH E-FIELD REGIME: Focus on species targets, not charge balance")
    print("=" * 80)
    print()
    print("User's approach:")
    print("  - Ne = 2.3e9 cm⁻³ (FIXED)")
    print("  - Te ~ 1.5 eV")
    print("  - E = 50-300 V/cm")
    print("  - Optimize: H, CH, C2 densities")
    print("  - Monitor Ni/Ne but don't constrain it")
    print()
    print("Focus: Hit species targets without charge balance constraint!")
    print()

    bounds = [
        (1.0, 2.0),      # Te: around 1.5 eV per user
        (50.0, 300.0),   # E_field: 50-300 V/cm
    ]

    rate_bounds = {
        'CH_CH_C2_H2_cm3_5_4': (5e-11, 2e-10),
        'e_C2H2_C2_H2_cm3_1_16': (3e-11, 1e-10),
        'C2HPlus_e_C2_H_cm3_6_18': (1e-07, 3e-07),
        'C_CH_C2_H_cm3_7_4': (5e-11, 2e-10),
        'CH_C_C2_H_cm3_7_9': (5e-11, 2e-10),
        'C2H_H_C2_H2_cm3_7_47': (5e-11, 2e-10),
        'C2H2_H_C2_H2_H_cm3_7_50': (5e-12, 2e-11),
        'C2_H_CH_C_cm3_7_6': (1e-12, 5e-11),
        'stick_C2_9_9': (1000, 8000),
        'loss_C2_11_3': (500, 3000),
        'stick_C2H2_9_11': (100, 2000),
        'loss_C2H2_11_19': (100, 2000),
        'CH3_CH3_C2H2_H2_H2_cm3_7_49': (1e-12, 5e-11),
        'C_CH3_C2_H2_H_cm3_7_8': (5e-11, 3e-10),
        'CH3_CH_C2H2_H2_cm3_7_16': (5e-11, 3e-10),
        'loss_CH_11_9': (500, 5000),
        'stick_CH_9_3': (500, 5000),
        'stick_CH3_9_2': (100, 2000),
        'e_CH4Plus_CH3_H_cm3_6_4': (5e-07, 1.5e-06),
    }

    bounds.extend(rate_bounds.values())

    print(f"Optimizing {len(bounds)} parameters:")
    print(f"  - Te (1.0-2.0 eV)")
    print(f"  - E-field (50-300 V/cm)")
    print(f"  - {len(rate_bounds)} reaction rates")
    print()

    result = differential_evolution(
        objective,
        bounds,
        maxiter=100,
        popsize=20,
        workers=1,
        updating='deferred',
        polish=False,
        seed=50
    )

    print("\n" + "=" * 80)
    print("OPTIMIZATION COMPLETE")
    print("=" * 80)
    print(f"Best f(x) = {result.fun:.2f}")
    print(f"Total evaluations: {eval_count}")
