#!/usr/bin/env python3
"""
Manual sweep of C2_H_CH_C (H + C2 -> CH + C) reaction rate
This reaction accounts for 28% of C2 consumption but has NEVER been optimized!
Starting from best known parameters (checkpoint_f3407.json)
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from odefun import PlasmaODE
from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from rate_database_complete import get_complete_rate_database

# Load best known parameters
with open('checkpoint_f3407.json', 'r') as f:
    best = json.load(f)

params_base = best['params'].copy()

# Target densities
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

# Get literature range for C2_H_CH_C
db = get_complete_rate_database()
baseline_rate = 9.6e-11  # From define_rates_tunable.py
lit_min = db['C2_H_CH_C_cm3_7_6'].min  # 8.0e-11
lit_max = db['C2_H_CH_C_cm3_7_6'].max  # 1.2e-10
lit_baseline = db['C2_H_CH_C_cm3_7_6'].value  # 9.6e-11

# Test these C2_H_CH_C values (in cm^3/s)
# Try reducing the rate to save C2 from being consumed
c2_h_ch_c_values = [
    8.0e-11,   # Literature minimum
    8.5e-11,
    9.0e-11,
    9.6e-11,   # Baseline/Literature value
    1.0e-10,
    1.1e-10,
    1.2e-10,   # Literature maximum
]

print("="*80)
print("MANUAL SWEEP: C2_H_CH_C (H + C2 -> CH + C) reaction rate")
print("="*80)
print(f"\nThis reaction accounts for 28% of C2 consumption (2nd largest loss)")
print(f"Has NEVER been included in rate_values optimization!")
print()
print(f"Best known result (f(x)={best['objective']:.1f}):")
print(f"  C2_H_CH_C = {baseline_rate:.2e} (baseline, never tuned)")
print(f"  H:  {best['densities']['H']:.2e} ({best['densities']['H']/TARGETS['H']:.2f}× target)")
print(f"  C2: {best['densities']['C2']:.2e} ({best['densities']['C2']/TARGETS['C2']:.2f}× target)")
print(f"  CH: {best['densities']['CH']:.2e} ({best['densities']['CH']/TARGETS['CH']:.2f}× target)")
print()
print(f"Literature range: [{lit_min:.2e}, {lit_max:.2e}] cm³/s")
print(f"Testing values: {[f'{v:.2e}' for v in c2_h_ch_c_values]}")
print("="*80)
print()

results = []

for c2_h_rate in c2_h_ch_c_values:
    print(f"Testing C2_H_CH_C = {c2_h_rate:.2e}...", end='', flush=True)

    # Create params with modified C2_H_CH_C
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()
    params['rate_values']['C2_H_CH_C_cm3_7_6'] = c2_h_rate

    try:
        # Get rate dictionary with modified C2_H_CH_C
        k_dict = define_rates_tunable(params)

        # CRITICAL: Apply rate_values overrides to k_dict
        for name, val in params.get('rate_values', {}).items():
            if name in k_dict and name in db:
                k_dict[name] = np.clip(val, db[name].min, db[name].max)

        params['k'] = k_dict

        # Build reactions (needs k in params)
        reactions, tags = build_reactions(params)

        # Add reactions to params for PlasmaODE
        params['R'] = reactions
        params['tags'] = tags

        # Setup ODE system
        ode = PlasmaODE(params)

        # Initial conditions
        y0 = np.zeros(ode.ns)
        y0[ode.species.index('e')] = params['ne']
        y0[ode.species.index('Ar')] = 1.29e16  # From P=0.4 Torr
        y0[ode.species.index('CH4')] = 1.29e15

        # Solve to steady-state
        sol = solve_ivp(
            ode,
            t_span=[0, 1e-3],
            y0=y0,
            method='BDF',
            rtol=1e-6,
            atol=1e-10,
            max_step=1e-5
        )

        if not sol.success:
            print(f" FAILED: {sol.message}")
            continue

        # Extract final densities
        y_final = sol.y[:, -1]
        H = y_final[ode.species.index('H')]
        CH = y_final[ode.species.index('CH')]
        C2 = y_final[ode.species.index('C2')]
        C2H2 = y_final[ode.species.index('C2H2')]

        # Calculate objective
        errors = []
        for sp, target in TARGETS.items():
            idx = ode.species.index(sp)
            density = y_final[idx]
            ratio = density / target
            if ratio < 0.6:
                err = (0.6 / ratio) ** 2
            elif ratio > 1.4:
                err = (ratio / 1.4) ** 2
            else:
                err = 1.0
            errors.append(err)

        objective = np.prod(errors)

        results.append({
            'c2_h_rate': c2_h_rate,
            'H': H,
            'CH': CH,
            'C2': C2,
            'C2H2': C2H2,
            'objective': objective
        })

        print(f" f(x)={objective:.1f} | H={H:.2e} ({H/TARGETS['H']:.2f}×) | C2={C2:.2e} ({C2/TARGETS['C2']:.2f}×) | CH={CH:.2e} ({CH/TARGETS['CH']:.2f}×)")

    except Exception as e:
        import traceback
        print(f" ERROR: {e}")
        traceback.print_exc()
        continue

print()
print("="*80)
print("SUMMARY")
print("="*80)
print(f"{'C2_H_CH_C':>12} | {'f(x)':>10} | {'H ratio':>8} | {'C2 ratio':>8} | {'CH ratio':>8}")
print("-"*80)

for r in results:
    marker = " ← baseline" if abs(r['c2_h_rate'] - baseline_rate) < 1e-13 else ""
    marker = " ← lit min" if abs(r['c2_h_rate'] - lit_min) < 1e-13 else marker
    marker = " ← lit max" if abs(r['c2_h_rate'] - lit_max) < 1e-13 else marker

    print(f"{r['c2_h_rate']:>12.2e} | {r['objective']:>10.1f} | "
          f"{r['H']/TARGETS['H']:>8.2f} | {r['C2']/TARGETS['C2']:>8.2f} | {r['CH']/TARGETS['CH']:>8.2f}{marker}")

# Find best
if results:
    best_result = min(results, key=lambda x: x['objective'])
    print()
    print(f"BEST: C2_H_CH_C = {best_result['c2_h_rate']:.2e} cm³/s, f(x) = {best_result['objective']:.1f}")
    print(f"  H:  {best_result['H']:.2e} ({best_result['H']/TARGETS['H']:.2f}× target)")
    print(f"  C2: {best_result['C2']:.2e} ({best_result['C2']/TARGETS['C2']:.2f}× target)")
    print(f"  CH: {best_result['CH']:.2e} ({best_result['CH']/TARGETS['CH']:.2f}× target)")
    print()

    # Compare to baseline
    baseline_result = [r for r in results if abs(r['c2_h_rate'] - baseline_rate) < 1e-13]
    if baseline_result:
        baseline_obj = baseline_result[0]['objective']
        improvement = baseline_obj / best_result['objective']
        print(f"Improvement vs baseline: {improvement:.2f}× better")
