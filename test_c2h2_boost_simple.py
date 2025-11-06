#!/usr/bin/env python3
"""
Boost C2H2 by reducing stick_C2H2 wall losses
Strategy: More C2H2 → More C2 production (92% of C2 comes from C2H2 + H)
Starting from best combined parameters (stick_C2 at min, loss_C2 at min, C2_H_CH_C at min)
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from odefun import PlasmaODE
from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from rate_database_complete import get_complete_rate_database

# Load checkpoint
with open('checkpoint_f3407.json', 'r') as f:
    checkpoint = json.load(f)

params_base = checkpoint['params'].copy()

# Target densities
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

# Get literature range for stick_C2H2
db = get_complete_rate_database()
stick_c2h2_checkpoint = params_base['rate_values']['stick_C2H2_9_11']  # 944
lit_min = db['stick_C2H2_9_11'].min  # Need to check
lit_max = db['stick_C2H2_9_11'].max  # Need to check

print(f"stick_C2H2 in checkpoint: {stick_c2h2_checkpoint:.2e}")
print(f"Literature range: [{lit_min:.2e}, {lit_max:.2e}]")
print()

# Test reduction factors for stick_C2H2
# Lower = less C2H2 wall loss = more C2H2 available = more C2 production
stick_c2h2_factors = [
    0.2,   # 80% reduction
    0.3,   # 70% reduction
    0.4,   # 60% reduction
    0.5,   # 50% reduction
    0.7,   # 30% reduction
    1.0,   # Baseline (current best)
]

print("="*80)
print("C2H2 BOOST TEST - Starting from BEST combined parameters")
print("="*80)
print(f"\nStrategy: Reduce stick_C2H2 to boost C2H2 density")
print(f"92% of C2 production comes from C2H2 + H → C2 + H2")
print()
print("Starting point (best combined):")
print(f"  stick_C2:   1.25e+03 (lit min)")
print(f"  loss_C2:    1.00e+02 (near min)")
print(f"  C2_H_CH_C:  8.00e-11 (lit min)")
print(f"  stick_C2H2: {stick_c2h2_checkpoint:.2e} (checkpoint)")
print()
print(f"Testing stick_C2H2 reduction factors: {stick_c2h2_factors}")
print("="*80)
print()

results = []

for factor in stick_c2h2_factors:
    stick_c2h2_val = stick_c2h2_checkpoint * factor
    print(f"Testing stick_C2H2 factor = {factor:.2f} ({stick_c2h2_val:.2e})...", end='', flush=True)

    # Start from best combined parameters
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply BEST parameters from previous tests
    params['rate_values']['stick_C2_9_9'] = 1.25e3      # lit min (was 5.47e3)
    params['rate_values']['loss_C2_11_3'] = 100         # near min (was 1199)
    params['rate_values']['C2_H_CH_C_cm3_7_6'] = 8.0e-11  # lit min (was 9.6e-11)

    # Apply C2H2 reduction
    params['rate_values']['stick_C2H2_9_11'] = stick_c2h2_val

    try:
        # Get rate dictionary
        k_dict = define_rates_tunable(params)

        # Apply rate_values overrides
        for name, val in params.get('rate_values', {}).items():
            if name in k_dict and name in db:
                k_dict[name] = np.clip(val, db[name].min, db[name].max)

        params['k'] = k_dict

        # Build reactions
        reactions, tags = build_reactions(params)
        params['R'] = reactions
        params['tags'] = tags

        # Setup ODE system
        ode = PlasmaODE(params)

        # Initial conditions
        y0 = np.zeros(ode.ns)
        y0[ode.species.index('e')] = params['ne']
        y0[ode.species.index('Ar')] = 1.29e16
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
            'factor': factor,
            'stick_c2h2': stick_c2h2_val,
            'H': H,
            'CH': CH,
            'C2': C2,
            'C2H2': C2H2,
            'objective': objective
        })

        print(f" f(x)={objective:.1f} | H={H:.2e} ({H/TARGETS['H']:.2f}×) | C2={C2:.2e} ({C2/TARGETS['C2']:.2f}×) | C2H2={C2H2:.2e} | CH={CH:.2e} ({CH/TARGETS['CH']:.2f}×)")

    except Exception as e:
        import traceback
        print(f" ERROR: {e}")
        traceback.print_exc()
        continue

print()
print("="*80)
print("SUMMARY")
print("="*80)
print(f"{'Factor':>8} | {'stick_C2H2':>12} | {'f(x)':>10} | {'H ratio':>8} | {'C2 ratio':>8} | {'C2H2':>12} | {'CH ratio':>8}")
print("-"*80)

for r in results:
    marker = " ← baseline" if abs(r['factor'] - 1.0) < 1e-6 else ""
    print(f"{r['factor']:>8.2f} | {r['stick_c2h2']:>12.2e} | {r['objective']:>10.1f} | "
          f"{r['H']/TARGETS['H']:>8.2f} | {r['C2']/TARGETS['C2']:>8.2f} | "
          f"{r['C2H2']:>12.2e} | {r['CH']/TARGETS['CH']:>8.2f}{marker}")

# Find best
if results:
    best_result = min(results, key=lambda x: x['objective'])
    print()
    print(f"BEST: Factor = {best_result['factor']:.2f}, f(x) = {best_result['objective']:.1f}")
    print(f"  H:    {best_result['H']:.2e} ({best_result['H']/TARGETS['H']:.2f}× target)")
    print(f"  C2:   {best_result['C2']:.2e} ({best_result['C2']/TARGETS['C2']:.2f}× target)")
    print(f"  C2H2: {best_result['C2H2']:.2e} (precursor)")
    print(f"  CH:   {best_result['CH']:.2e} ({best_result['CH']/TARGETS['CH']:.2f}× target)")
    print()

    # Compare to baseline
    baseline_result = [r for r in results if abs(r['factor'] - 1.0) < 1e-6]
    if baseline_result:
        baseline_obj = baseline_result[0]['objective']
        baseline_c2 = baseline_result[0]['C2']
        baseline_c2h2 = baseline_result[0]['C2H2']
        improvement = baseline_obj / best_result['objective']
        c2_improvement = best_result['C2'] / baseline_c2
        c2h2_improvement = best_result['C2H2'] / baseline_c2h2
        print(f"Improvement vs baseline (factor=1.0):")
        print(f"  Objective: {improvement:.2f}× better")
        print(f"  C2 density: {c2_improvement:.2f}× higher ({(c2_improvement-1)*100:+.0f}%)")
        print(f"  C2H2 density: {c2h2_improvement:.2f}× higher ({(c2h2_improvement-1)*100:+.0f}%)")
