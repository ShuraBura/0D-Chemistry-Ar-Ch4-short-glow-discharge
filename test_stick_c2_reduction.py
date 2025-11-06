#!/usr/bin/env python3
"""
Test reducing stick_C2 from checkpoint value (5.47e3) to literature minimum (1.25e3)
stick_C2 accounts for 59% of C2 consumption - this could be a MAJOR improvement!
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

# Get literature range for stick_C2
db = get_complete_rate_database()
checkpoint_val = params_base['rate_values']['stick_C2_9_9']  # 5.47e3
lit_min = db['stick_C2_9_9'].min  # 1.25e3
lit_max = db['stick_C2_9_9'].max  # 6.25e3

# Test these stick_C2 values
stick_c2_values = [
    lit_min,          # 1.25e3 - Literature minimum
    2.0e3,
    3.0e3,
    4.0e3,
    checkpoint_val,   # 5.47e3 - Current "optimized" value
    lit_max,          # 6.25e3 - Literature maximum
]

print("="*80)
print("STICK_C2 REDUCTION TEST - THE DOMINANT C2 LOSS MECHANISM")
print("="*80)
print(f"\nCRITICAL FINDING: stick_C2 accounts for 59% of C2 consumption")
print(f"Checkpoint has stick_C2 = {checkpoint_val:.2e} ({checkpoint_val/lit_min:.2f}× literature min!)")
print()
print(f"Best known result (f(x)={best['objective']:.1f}):")
print(f"  H:  {best['densities']['H']:.2e} ({best['densities']['H']/TARGETS['H']:.2f}× target)")
print(f"  C2: {best['densities']['C2']:.2e} ({best['densities']['C2']/TARGETS['C2']:.2f}× target)")
print(f"  CH: {best['densities']['CH']:.2e} ({best['densities']['CH']/TARGETS['CH']:.2f}× target)")
print()
print(f"Literature range: [{lit_min:.2e}, {lit_max:.2e}]")
print(f"Testing values: {[f'{v:.2e}' for v in stick_c2_values]}")
print("="*80)
print()

results = []

for stick_c2 in stick_c2_values:
    print(f"Testing stick_C2 = {stick_c2:.2e}...", end='', flush=True)

    # Create params with modified stick_C2
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()
    params['rate_values']['stick_C2_9_9'] = stick_c2

    try:
        # Get rate dictionary
        k_dict = define_rates_tunable(params)

        # CRITICAL: Apply rate_values overrides
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
            'stick_c2': stick_c2,
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
print(f"{'stick_C2':>12} | {'f(x)':>10} | {'H ratio':>8} | {'C2 ratio':>8} | {'CH ratio':>8}")
print("-"*80)

for r in results:
    marker = ""
    if abs(r['stick_c2'] - checkpoint_val) < 1:
        marker = " ← checkpoint"
    elif abs(r['stick_c2'] - lit_min) < 1:
        marker = " ← LIT MIN"
    elif abs(r['stick_c2'] - lit_max) < 1:
        marker = " ← LIT MAX"

    print(f"{r['stick_c2']:>12.2e} | {r['objective']:>10.1f} | "
          f"{r['H']/TARGETS['H']:>8.2f} | {r['C2']/TARGETS['C2']:>8.2f} | {r['CH']/TARGETS['CH']:>8.2f}{marker}")

# Find best
if results:
    best_result = min(results, key=lambda x: x['objective'])
    print()
    print(f"BEST: stick_C2 = {best_result['stick_c2']:.2e}, f(x) = {best_result['objective']:.1f}")
    print(f"  H:  {best_result['H']:.2e} ({best_result['H']/TARGETS['H']:.2f}× target)")
    print(f"  C2: {best_result['C2']:.2e} ({best_result['C2']/TARGETS['C2']:.2f}× target)")
    print(f"  CH: {best_result['CH']:.2e} ({best_result['CH']/TARGETS['CH']:.2f}× target)")
    print()

    # Compare to checkpoint
    checkpoint_result = [r for r in results if abs(r['stick_c2'] - checkpoint_val) < 1]
    if checkpoint_result:
        checkpoint_obj = checkpoint_result[0]['objective']
        checkpoint_c2 = checkpoint_result[0]['C2']
        improvement = checkpoint_obj / best_result['objective']
        c2_improvement = best_result['C2'] / checkpoint_c2
        print(f"Improvement vs checkpoint:")
        print(f"  Objective: {improvement:.2f}× better")
        print(f"  C2 density: {c2_improvement:.2f}× higher ({c2_improvement - 1:.0%} increase)")
        print()
        print(f"This makes sense! Lower stick_C2 = less C2 wall loss = more C2 retained!")
