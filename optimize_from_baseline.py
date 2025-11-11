#!/usr/bin/env python3
"""
Optimizer starting from BASELINE SUCCESS CASE.

Baseline (best_f70.3.json):
  H:  2.01e14 (79.9%)  ✓
  CH: 1.01e09 (100.7%) ✓✓✓
  C2: 9.41e08 (16.8%)  ← IMPROVE THIS
  Ni/Ne: 3.12 ✓

  Te: 1.31 eV, Ne: 1.22e8, E: 250 V/cm
  C2H2: 4.81e09 (need 5e12, 1000× higher!)
  CH3: 6.85e11 (need 2e13, 30× higher)

Strategy: Start near these conditions and push C2 higher
by increasing CH3/C2H2 production without breaking H/CH balance.
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import os

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

os.makedirs('optimization_results_from_baseline', exist_ok=True)

# Load baseline
with open('optimization_results_comprehensive_1e12/best_f70.3.json', 'r') as f:
    baseline = json.load(f)

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

# Define critical reactions list at module level (used in both objective and bounds setup)
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
    'C_CH3_C2H2_H_cm3_7_8',
    'CH3_CH_C2H2_H2_cm3_7_16',
    'e_CH4Plus_CH3_H_cm3_6_4',
    'e_CH4Plus_CH2_H2_cm3_6_9',
    'ArStar_CH4_CH3_H_cm3_3_1',
    'stick_ArStar_9_5',
    'e_CH4_CH3_HMinus_cm3_8_1',
    'loss_CH_11_9',
    'stick_CH_9_3',
    'stick_CH3_9_2',
]

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

    # Apply tunable reactions (critical_reactions defined at module level)
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

        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus', 'CHPlus', 'H2Plus']
        n_i_total = sum(y_final[species.index(ion)] for ion in ions)
        ne_final = y_final[species.index('e')]
        Ni_over_Ne = n_i_total / ne_final if ne_final > 0 else 0

        # Balance all three targets - baseline already achieved H and CH
        weights = {'H': 20.0, 'CH': 30.0, 'C2': 50.0}  # C2 weight highest, but keep CH important too

        err_H = ((H_final - target_densities['H']) / target_densities['H'])**2
        err_CH = ((CH_final - target_densities['CH']) / target_densities['CH'])**2
        err_C2 = ((C2_final - target_densities['C2']) / target_densities['C2'])**2

        f = (weights['H'] * err_H + weights['CH'] * err_CH + weights['C2'] * err_C2)

        # Charge balance constraint (user said monitor but don't stress)
        # Use soft constraint: penalty only if really bad
        if Ni_over_Ne < 1.5 or Ni_over_Ne > 10.0:
            charge_penalty = 50.0 * min((1.5 - Ni_over_Ne)**2, (Ni_over_Ne - 10.0)**2)
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

            filename = f'optimization_results_from_baseline/best_f{f:.1f}.json'
            with open(filename, 'w') as f_out:
                json.dump(result, f_out, indent=2)

            print(f"\nEval {eval_count}: f(x) = {f:.2f}")
            print(f"  Te:     {Te:.2f} eV")
            print(f"  Ne:     {ne_final:.2e}")
            print(f"  E:      {E_field:.1f} V/cm")
            print(f"  H:      {H_final:.2e} ({H_final/target_densities['H']*100:5.1f}%)")
            print(f"  CH:     {CH_final:.2e} ({CH_final/target_densities['CH']*100:5.1f}%)")
            print(f"  C2:     {C2_final:.2e} ({C2_final/target_densities['C2']*100:5.1f}%)")
            print(f"  C2H2:   {C2H2_final:.2e}")
            print(f"  CH3:    {CH3_final:.2e}")
            print(f"  Ni/Ne:  {Ni_over_Ne:.2f}")

        return f

    except Exception as e:
        print(f"Error in eval {eval_count}: {e}")
        return 1e6

if __name__ == '__main__':
    print("=" * 80)
    print("OPTIMIZE FROM BASELINE SUCCESS CASE")
    print("=" * 80)
    print()
    print("Baseline:")
    print(f"  H:  {baseline['target_densities']['H']:.2e} (79.9%) ✓")
    print(f"  CH: {baseline['target_densities']['CH']:.2e} (100.7%) ✓✓✓")
    print(f"  C2: {baseline['target_densities']['C2']:.2e} (16.8%) ← IMPROVE")
    print(f"  Te: {baseline['Te']:.2f} eV, E: {baseline['E_field']:.1f} V/cm")
    print(f"  Ni/Ne: {baseline['Ni_over_Ne']:.2f} ✓")
    print()
    print("Strategy: Start near baseline conditions, push C2 higher")
    print("  - C2 weight 4× higher than H/CH in objective")
    print("  - Soft charge balance constraint (Ni/Ne: 1.5-10 okay)")
    print()

    # Start search NEAR baseline conditions
    bounds = [
        (1.0, 2.0),       # Te: baseline was 1.31, allow 1.0-2.0
        (7.5, 9.5),       # log10(ne): baseline was log10(1.22e8)=8.09
        (150.0, 300.0),   # E_field: baseline was 250
    ]

    # Rate bounds - MUST be in same order as critical_reactions list in objective()!
    baseline_rates = baseline['rate_values']

    # Set bounds for each rate in the SAME order as critical_reactions
    rate_bounds_ordered = []

    for rate_name in critical_reactions:
        val = baseline_rates.get(rate_name, None)

        # Chemistry reactions (cm³/s)
        if '_cm3_' in rate_name:
            if val is None:
                val = 1e-10  # Default for chemistry
            rate_bounds_ordered.append((val * 0.2, val * 5.0))

        # Loss/sticking (s⁻¹)
        elif 'stick_' in rate_name or 'loss_' in rate_name:
            if val is None:
                val = 1000  # Default for losses
            rate_bounds_ordered.append((val * 0.1, val * 3.0))

        else:
            # Fallback
            if val is None:
                val = 1e-10
            rate_bounds_ordered.append((val * 0.1, val * 10.0))

    bounds.extend(rate_bounds_ordered)

    print(f"Optimizing {len(bounds)} parameters")
    print(f"  - 3 plasma params (near baseline)")
    print(f"  - 23 reaction rates (near baseline values)")
    print()

    result = differential_evolution(
        objective,
        bounds,
        maxiter=100,
        popsize=20,
        workers=1,
        updating='deferred',
        polish=False,
        seed=51
    )

    print("\n" + "=" * 80)
    print("OPTIMIZATION COMPLETE")
    print("=" * 80)
    print(f"Best f(x) = {result.fun:.2f}")
    print(f"Total evaluations: {eval_count}")
