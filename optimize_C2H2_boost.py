"""
STAGE 1: Boost C2H2 from 4.81e9 to ~5e10 (10× increase)

Key insight from H vs C2 analysis:
- H + C2H2 → C2 is 14× faster than H + C2 → CH
- Higher H HELPS C2, doesn't hurt it!
- At C2H2 = 5e10, C2 reaches 92% of target (baseline H)

Strategy:
1. Start from baseline conditions (which achieve H=79.9%, CH=100.7%)
2. Only tune rates related to C2H2 production:
   - CH3 production: e + CH4, Ar* + CH4, e + CH4+
   - C2H2 formation: CH3 + CH3, C + CH3, CH3 + CH
   - CH3 losses (reduce to keep CH3 high)
3. Lock all other rates at baseline values
4. Target: C2H2 ~ 5e10, maintain H and CH balance

C2H2 formation: C2H2 ∝ [CH3]²
Current CH3: 6.85e11
Target CH3: ~2e12 (3× increase → 9× C2H2 increase)
"""

import sys
sys.path.append('/home/user/0D-Chemistry-Ar-Ch4-short-glow-discharge')

import numpy as np
import json
import os
from scipy.optimize import differential_evolution
from scipy.integrate import solve_ivp

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

os.makedirs('optimization_results_C2H2_boost', exist_ok=True)

def pressure_to_ntotal(pressure_mTorr, T_K=300):
    """Convert pressure to total density"""
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    return pressure_Pa / (k_B * T_K)

# Load baseline
with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

# Targets
target_densities = {
    'H': 2.52e14,
    'CH': 1.0e9,
    'C2': 5.6e11,
    'C2H2': 5e10,  # NEW TARGET!
}

eval_count = 0
best_f = np.inf

# Only tune these rates (CH3 and C2H2 production)
tunable_rates = [
    # CH3 production
    'e_CH4_CH3_HMinus_cm3_8_1',       # e + CH4 → CH3 + H⁻ (Te-dependent)
    'ArStar_CH4_CH3_H_cm3_3_1',       # Ar* + CH4 → CH3 + H
    'e_CH4Plus_CH3_H_cm3_6_4',        # e + CH4+ → CH3 + H

    # C2H2 formation
    'CH3_CH3_C2H2_H2_H2_cm3_7_49',    # CH3 + CH3 → C2H2 + 2H2 (main!)
    'C_CH3_C2H2_H_cm3_7_8',           # C + CH3 → C2H2 + H
    'CH3_CH_C2H2_H2_cm3_7_16',        # CH3 + CH → C2H2 + H2

    # CH3 losses (reduce these to keep CH3 high)
    'stick_CH3_9_2',                  # CH3 wall sticking

    # C2H2 losses (reduce to accumulate C2H2)
    'stick_C2H2_9_11',                # C2H2 wall sticking
    'loss_C2H2_11_19',                # C2H2 volumetric loss
]

def objective(x):
    global eval_count, best_f
    eval_count += 1

    Te = x[0]
    ne = 10**x[1]
    E_field = x[2]

    n_total = pressure_to_ntotal(500.0, 300.0)

    pressure_mTorr = 500.0

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
                    'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH',
                    'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C',
                    'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H',
                    'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus',
                    'H2Plus', 'C2H2Star'],
    }

    k = define_rates(params)

    # Apply tunable rates
    for i, rate_name in enumerate(tunable_rates):
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
            return 1e10

        y_final = sol.y[:, -1]

        H_final = y_final[species.index('H')]
        CH_final = y_final[species.index('CH')]
        C2_final = y_final[species.index('C2')]
        C2H2_final = y_final[species.index('C2H2')]
        CH3_final = y_final[species.index('CH3')]

        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                'H3Plus', 'CHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                'C2H3Plus', 'C2HPlus', 'H2Plus']
        n_i_total = sum(y_final[species.index(ion)] for ion in ions)
        ne_final = y_final[species.index('e')]
        Ni_over_Ne = n_i_total / ne_final if ne_final > 0 else 0

        # NEW OBJECTIVE: Prioritize C2H2 boost!
        weights = {'H': 10.0, 'CH': 20.0, 'C2H2': 50.0, 'C2': 30.0}

        err_H = ((H_final - target_densities['H']) / target_densities['H'])**2
        err_CH = ((CH_final - target_densities['CH']) / target_densities['CH'])**2
        err_C2H2 = ((C2H2_final - target_densities['C2H2']) / target_densities['C2H2'])**2
        err_C2 = ((C2_final - target_densities['C2']) / target_densities['C2'])**2

        f = (weights['H'] * err_H +
             weights['CH'] * err_CH +
             weights['C2H2'] * err_C2H2 +
             weights['C2'] * err_C2)

        # Charge balance constraint
        if Ni_over_Ne < 1.5 or Ni_over_Ne > 10.0:
            penalty = 1000 * (min(1.5 - Ni_over_Ne, 0)**2 + max(Ni_over_Ne - 10.0, 0)**2)
            f += penalty

        if f < best_f:
            best_f = f
            print(f"\nEval {eval_count}: f(x) = {f:.2f}")
            print(f"  Te:     {Te:.2f} eV")
            print(f"  Ne:     {ne:.2e}")
            print(f"  E:      {E_field:.1f} V/cm")
            print(f"  H:      {H_final:.2e} ({H_final/target_densities['H']*100:6.1f}%)")
            print(f"  CH:     {CH_final:.2e} ({CH_final/target_densities['CH']*100:6.1f}%)")
            print(f"  C2H2:   {C2H2_final:.2e} ({C2H2_final/target_densities['C2H2']*100:6.1f}%) ← TARGET")
            print(f"  C2:     {C2_final:.2e} ({C2_final/target_densities['C2']*100:6.1f}%)")
            print(f"  CH3:    {CH3_final:.2e}")
            print(f"  Ni/Ne:  {Ni_over_Ne:.2f}")

            # Save result
            result = {
                'pressure_mTorr': 500.0,
                'n_total': n_total,
                'Te': Te,
                'Ne': ne,
                'E_field': E_field,
                'n_i_total': n_i_total,
                'Ni_over_Ne': Ni_over_Ne,
                'charge_imbalance_pct': abs(Ni_over_Ne - 1.0) * 100,
                'rate_values': {rate: k[rate] for rate in tunable_rates if rate in k},
                'target_densities': {
                    'H': float(H_final),
                    'CH': float(CH_final),
                    'C2': float(C2_final),
                    'C2H2': float(C2H2_final),
                    'CH3': float(CH3_final),
                },
                'all_densities': {sp: float(y_final[species.index(sp)]) for sp in species},
            }

            with open(f'optimization_results_C2H2_boost/best_f{f:.1f}.json', 'w') as outf:
                json.dump(result, outf, indent=2)

        return f

    except Exception as e:
        return 1e10

if __name__ == '__main__':
    print("=" * 80)
    print("STAGE 1: BOOST C2H2 FROM 4.81e9 TO ~5e10")
    print("=" * 80)
    print()
    print("Baseline:")
    print(f"  H:  {baseline['target_densities']['H']:.2e} ({baseline['target_densities']['H']/target_densities['H']*100:.1f}%) ✓")
    print(f"  CH: {baseline['target_densities']['CH']:.2e} ({baseline['target_densities']['CH']/target_densities['CH']*100:.1f}%) ✓✓✓")
    print(f"  C2H2: {baseline['all_densities']['C2H2']:.2e} ← BOOST 10×")
    print(f"  C2: {baseline['target_densities']['C2']:.2e} ({baseline['target_densities']['C2']/target_densities['C2']*100:.1f}%)")
    print(f"  Te: {baseline['Te']:.2f} eV, E: {baseline['E_field']:.1f} V/cm")
    print(f"  Ni/Ne: {baseline['Ni_over_Ne']:.2f} ✓")
    print()
    print("Strategy: Only tune CH3 and C2H2 production pathways")
    print(f"  - {len(tunable_rates)} tunable rates (vs 23 in previous attempt)")
    print("  - All other rates locked at baseline values")
    print("  - Target: C2H2 = 5e10 (92% → C2 target)")
    print()

    # Bounds centered on baseline
    baseline_rates = baseline['rate_values']

    bounds = [
        (1.0, 2.0),       # Te: baseline was 1.31
        (7.5, 9.0),       # log10(ne): baseline was 8.09
        (150.0, 300.0),   # E_field: baseline was 250
    ]

    # Rate bounds - allow significant variation to boost CH3/C2H2
    for rate_name in tunable_rates:
        val = baseline_rates.get(rate_name, None)

        if '_cm3_' in rate_name:
            if val is None:
                val = 1e-10
            # Allow more boost for production rates
            if 'CH3' in rate_name or 'C2H2' in rate_name:
                bounds.append((val * 0.1, val * 10.0))  # Wide range
            else:
                bounds.append((val * 0.5, val * 2.0))

        elif 'stick_' in rate_name or 'loss_' in rate_name:
            if val is None:
                val = 1000
            # Allow reducing losses significantly
            bounds.append((val * 0.01, val * 2.0))  # Can go down to 1%

        else:
            if val is None:
                val = 1e-10
            bounds.append((val * 0.1, val * 10.0))

    print(f"Optimizing {len(bounds)} parameters")
    print(f"  - 3 plasma params")
    print(f"  - {len(tunable_rates)} rates (CH3 & C2H2 production only)")
    print()

    result = differential_evolution(
        objective,
        bounds,
        maxiter=300,
        popsize=15,
        seed=42,
        workers=1,
        updating='deferred',
        atol=1e-3,
        tol=0.01,
    )

    print("\n" + "=" * 80)
    print("OPTIMIZATION COMPLETE")
    print("=" * 80)
    print(f"Best f(x) = {result.fun:.2f}")
    print(f"Total evaluations: {eval_count}")
