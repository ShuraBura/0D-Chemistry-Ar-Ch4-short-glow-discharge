#!/usr/bin/env python3
"""
Optimization with CH constrained to lower values to force carbon redistribution.

Key insight: CH is overshooting to 3510% at steady state. By constraining CH
to be closer to target, we force carbon to redistribute into other pathways,
potentially boosting C2H2 and then C2 production.
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import os

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

# Create output directory
os.makedirs('optimization_results_CH_constrained', exist_ok=True)

def pressure_to_density(pressure_mTorr, T_K=400):
    kB = 1.38064852e-23
    Torr_to_Pa = 133.322
    P_Pa = pressure_mTorr * 1e-3 * Torr_to_Pa
    n_m3 = P_Pa / (kB * T_K)
    return n_m3 * 1e-6

# Target densities from experiment
target_densities = {
    'H': 2.52e14,    # cm⁻³
    'CH': 1.0e9,     # cm⁻³
    'C2': 5.6e11,    # cm⁻³
}

# CH constraint: penalize if CH > 5× target (currently at 35× target!)
CH_max_allowed = 5.0 * target_densities['CH']  # 5e+9

# C2H2 target from user's chemical insight
C2H2_min_required = 5e12

# Global counter
eval_count = 0
best_f = np.inf

def objective(x):
    """
    Objective function for optimization.
    x = [Te, log10(ne), E_field, tunable_rate_1, tunable_rate_2, ...]
    """
    global eval_count, best_f
    eval_count += 1

    # Extract parameters
    Te = x[0]
    ne = 10**x[1]
    E_field = x[2]

    # Plasma parameters
    pressure_mTorr = 500.0
    n_total = pressure_to_density(pressure_mTorr)

    params = {
        'P': pressure_mTorr,
        'Te': Te,
        'ne': ne,
        'E_field': E_field,
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

    # Define all rate constants
    k = define_rates(params)

    # Critical reactions to include in optimization (C2 producers + C2H2 pathways)
    critical_reactions = [
        # C2 production reactions
        'CH_CH_C2_H2_cm3_5_4',
        'e_C2H2_C2_H2_cm3_1_16',
        'C2HPlus_e_C2_H_cm3_6_18',
        'C_CH_C2_H_cm3_7_4',
        'CH_C_C2_H_cm3_7_9',
        'C2H_H_C2_H2_cm3_7_47',

        # C2H2 loss rates (CRITICAL for C2H2 accumulation!)
        'stick_C2H2_9_11',
        'loss_C2H2_11_19',

        # CH loss rates (NEW: reduce CH by increasing its loss!)
        'loss_CH_11_9',
        'stick_CH_9_3',

        # Top C2H2 producers (boost these to increase C2H2)
        'CH3_CH3_C2H2_H2_H2_cm3_7_49',  # 2CH3 → 2H2 + C2H2 (67% of C2H2 production!)
        'C_CH3_C2_H2_H_cm3_7_8',         # C + CH3 → H + C2H2 (9%)
        'CH3_CH_C2H2_H2_cm3_7_16',       # CH + CH3 → H2 + C2H2 (7%)
    ]

    # Apply tunable rates
    for i, rate_name in enumerate(critical_reactions):
        if rate_name in k:
            k[rate_name] = x[3 + i]

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    # Initial conditions (start from typical discharge)
    species = params['species']
    y0 = np.ones(len(species)) * 1e3
    y0[species.index('Ar')] = n_total * 0.85
    y0[species.index('CH4')] = n_total * 0.15
    y0[species.index('e')] = ne

    # Integrate to steady state (use t=500s for speed)
    try:
        ode_func = PlasmaODE_Optimized(params)
        sol = solve_ivp(
            ode_func,
            (0, 500),
            y0,
            method='BDF',
            rtol=1e-7,
            atol=1e-9,
            max_step=1.0
        )

        if not sol.success:
            return 1e6

        y_final = sol.y[:, -1]

        # Extract target species
        H_final = y_final[species.index('H')]
        CH_final = y_final[species.index('CH')]
        C2_final = y_final[species.index('C2')]
        C2H2_final = y_final[species.index('C2H2')]

        # Calculate ion density for charge balance
        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus',
                'CHPlus', 'H2Plus']
        n_i_total = sum(y_final[species.index(ion)] for ion in ions)
        ne_final = y_final[species.index('e')]

        Ni_over_Ne = n_i_total / ne_final if ne_final > 0 else 0

        # Objective function components
        weights = {
            'H': 20.0,
            'CH': 20.0,
            'C2': 20.0
        }

        # Normalized errors
        err_H = ((H_final - target_densities['H']) / target_densities['H'])**2
        err_CH = ((CH_final - target_densities['CH']) / target_densities['CH'])**2
        err_C2 = ((C2_final - target_densities['C2']) / target_densities['C2'])**2

        f = (weights['H'] * err_H +
             weights['CH'] * err_CH +
             weights['C2'] * err_C2)

        # Charge balance constraint (Ni/Ne should be 2-7)
        if Ni_over_Ne < 2.0 or Ni_over_Ne > 7.0:
            charge_penalty = 200.0 * min((2.0 - Ni_over_Ne)**2, (Ni_over_Ne - 7.0)**2)
            f += charge_penalty

        # NEW: CH constraint - penalize if CH > 5× target
        if CH_final > CH_max_allowed:
            CH_penalty = 100.0 * ((CH_final - CH_max_allowed) / target_densities['CH'])**2
            f += CH_penalty

        # C2H2 constraint from user's chemical insight
        if C2H2_final < C2H2_min_required:
            C2H2_penalty = 100.0 * ((C2H2_min_required - C2H2_final) / C2H2_min_required)**2
            f += C2H2_penalty

        # Save if best
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
                'target_densities': {
                    'H': H_final,
                    'CH': CH_final,
                    'C2': C2_final
                },
                'all_densities': {sp: y_final[i] for i, sp in enumerate(species)}
            }

            filename = f'optimization_results_CH_constrained/best_f{f:.1f}.json'
            with open(filename, 'w') as f_out:
                json.dump(result, f_out, indent=2)

            print(f"\nEval {eval_count}: f(x) = {f:.2f}")
            print(f"  H:    {H_final:.2e} ({H_final/target_densities['H']*100:5.1f}%)")
            print(f"  CH:   {CH_final:.2e} ({CH_final/target_densities['CH']*100:5.1f}% - max allowed: 500%)")
            print(f"  C2:   {C2_final:.2e} ({C2_final/target_densities['C2']*100:5.1f}%)")
            print(f"  C2H2: {C2H2_final:.2e} ({C2H2_final/C2H2_min_required*100:5.1f}% of 5e+12)")
            print(f"  Ni/Ne: {Ni_over_Ne:.2f}")

        return f

    except Exception as e:
        print(f"Error in eval {eval_count}: {e}")
        return 1e6

if __name__ == '__main__':
    print("=" * 80)
    print("OPTIMIZATION WITH CH CONSTRAINT")
    print("=" * 80)
    print()
    print("Strategy: Constrain CH ≤ 5× target to force carbon redistribution")
    print("  - Increase CH loss rates (stick_CH, loss_CH) - tunable")
    print("  - Boost C2H2 producers (CH3+CH3, C+CH3, CH+CH3) - tunable")
    print("  - Maintain C2H2 ≥ 5e+12 constraint")
    print()

    # Define bounds for optimization
    # [Te, log10(ne), E_field, ...rate constants...]
    bounds = [
        (0.5, 3.0),      # Te (eV)
        (8.0, 10.5),     # log10(ne)
        (10.0, 300.0),   # E_field (V/cm)
    ]

    # Add bounds for critical reactions
    rate_bounds = {
        # C2-producing reactions (standard bounds)
        'CH_CH_C2_H2_cm3_5_4': (5e-11, 2e-10),
        'e_C2H2_C2_H2_cm3_1_16': (3e-11, 1e-10),
        'C2HPlus_e_C2_H_cm3_6_18': (1e-07, 3e-07),
        'C_CH_C2_H_cm3_7_4': (5e-11, 2e-10),
        'CH_C_C2_H_cm3_7_9': (5e-11, 2e-10),
        'C2H_H_C2_H2_cm3_7_47': (5e-11, 2e-10),

        # C2H2 loss rates (allow low values for accumulation)
        'stick_C2H2_9_11': (500, 2000),
        'loss_C2H2_11_19': (1000, 2000),

        # CH loss rates (NEW: allow HIGH values to reduce CH!)
        'loss_CH_11_9': (1000, 5000),    # Increase range to allow more CH loss
        'stick_CH_9_3': (1000, 5000),    # Increase range to allow more CH loss

        # C2H2 producers (allow boost)
        'CH3_CH3_C2H2_H2_H2_cm3_7_49': (5e-12, 5e-11),
        'C_CH3_C2_H2_H_cm3_7_8': (5e-11, 2e-10),
        'CH3_CH_C2H2_H2_cm3_7_16': (5e-11, 2e-10),
    }

    critical_reactions = list(rate_bounds.keys())
    bounds.extend(rate_bounds.values())

    print(f"Optimizing {len(bounds)} parameters...")
    print(f"  - 3 plasma parameters (Te, Ne, E)")
    print(f"  - {len(critical_reactions)} reaction rates")
    print()

    # Run optimization
    result = differential_evolution(
        objective,
        bounds,
        maxiter=50,
        popsize=15,
        workers=1,
        updating='deferred',
        polish=False,
        seed=42
    )

    print("\n" + "=" * 80)
    print("OPTIMIZATION COMPLETE")
    print("=" * 80)
    print(f"Best f(x) = {result.fun:.2f}")
    print(f"Total evaluations: {eval_count}")
