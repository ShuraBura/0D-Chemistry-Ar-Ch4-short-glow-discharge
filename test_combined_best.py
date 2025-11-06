#!/usr/bin/env python3
"""
Test combining ALL the best parameters discovered:
1. stick_C2 = 1.25e3 (lit min) - 72% improvement alone!
2. loss_C2 = 100 (near min) - ~14% improvement
3. C2_H_CH_C = 8.0e-11 (lit min) - ~11% improvement
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from odefun import PlasmaODE
from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions
from rate_database_complete import get_complete_rate_database

# Load checkpoint as baseline
with open('checkpoint_f3407.json', 'r') as f:
    checkpoint = json.load(f)

params_base = checkpoint['params'].copy()

# Target densities
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

db = get_complete_rate_database()

print("="*80)
print("COMBINED OPTIMIZATION - Best parameters from all tests")
print("="*80)
print()

# Test configurations
configs = [
    {
        'name': 'Checkpoint (baseline)',
        'stick_C2': checkpoint['params']['rate_values']['stick_C2_9_9'],  # 5.47e3
        'loss_C2': checkpoint['params']['rate_values']['loss_C2_11_3'],    # 1199
        'C2_H_CH_C': 9.6e-11,  # baseline (not in rate_values)
    },
    {
        'name': 'Only stick_C2 optimized',
        'stick_C2': 1.25e3,  # lit min ← 72% improvement!
        'loss_C2': checkpoint['params']['rate_values']['loss_C2_11_3'],
        'C2_H_CH_C': 9.6e-11,
    },
    {
        'name': 'stick_C2 + loss_C2 optimized',
        'stick_C2': 1.25e3,  # lit min
        'loss_C2': 100,      # near min ← ~14% improvement
        'C2_H_CH_C': 9.6e-11,
    },
    {
        'name': 'ALL THREE optimized',
        'stick_C2': 1.25e3,      # lit min
        'loss_C2': 100,          # near min
        'C2_H_CH_C': 8.0e-11,    # lit min ← ~11% improvement
    },
]

results = []

for config in configs:
    print(f"Testing: {config['name']}...")
    print(f"  stick_C2:   {config['stick_C2']:.2e}")
    print(f"  loss_C2:    {config['loss_C2']:.2e}")
    print(f"  C2_H_CH_C:  {config['C2_H_CH_C']:.2e}")

    # Create params
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply optimized parameters
    params['rate_values']['stick_C2_9_9'] = config['stick_C2']
    params['rate_values']['loss_C2_11_3'] = config['loss_C2']
    params['rate_values']['C2_H_CH_C_cm3_7_6'] = config['C2_H_CH_C']

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
            print(f"  FAILED: {sol.message}\n")
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

        result = {
            'name': config['name'],
            'H': H,
            'CH': CH,
            'C2': C2,
            'C2H2': C2H2,
            'objective': objective,
            'H_ratio': H / TARGETS['H'],
            'CH_ratio': CH / TARGETS['CH'],
            'C2_ratio': C2 / TARGETS['C2'],
        }
        results.append(result)

        print(f"  f(x) = {objective:.1f}")
        print(f"  H:    {H:.2e} ({H/TARGETS['H']:.2f}× target)")
        print(f"  CH:   {CH:.2e} ({CH/TARGETS['CH']:.2f}× target)")
        print(f"  C2:   {C2:.2e} ({C2/TARGETS['C2']:.2f}× target)")
        print(f"  C2H2: {C2H2:.2e}")
        print()

    except Exception as e:
        import traceback
        print(f"  ERROR: {e}\n")
        traceback.print_exc()
        continue

print("="*80)
print("SUMMARY - Best Achievable Densities")
print("="*80)
print()

if results:
    best = results[-1]  # Last one should be "ALL THREE optimized"
    baseline = results[0]

    print(f"BASELINE (checkpoint f(x)={baseline['objective']:.1f}):")
    print(f"  H:  {baseline['H']:.2e} ({baseline['H_ratio']:.2f}× target)")
    print(f"  CH: {baseline['CH']:.2e} ({baseline['CH_ratio']:.2f}× target)")
    print(f"  C2: {baseline['C2']:.2e} ({baseline['C2_ratio']:.2f}× target)")
    print()

    print(f"BEST ACHIEVABLE (all optimizations combined, f(x)={best['objective']:.1f}):")
    print(f"  H:  {best['H']:.2e} ({best['H_ratio']:.2f}× target)")
    print(f"  CH: {best['CH']:.2e} ({best['CH_ratio']:.2f}× target)")
    print(f"  C2: {best['C2']:.2e} ({best['C2_ratio']:.2f}× target)")
    print()

    print("IMPROVEMENTS:")
    print(f"  Objective: {baseline['objective']/best['objective']:.2f}× better")
    print(f"  H:  {best['H']/baseline['H']:.2f}× ({(best['H']/baseline['H']-1)*100:+.0f}%)")
    print(f"  CH: {best['CH']/baseline['CH']:.2f}× ({(best['CH']/baseline['CH']-1)*100:+.0f}%)")
    print(f"  C2: {best['C2']/baseline['C2']:.2f}× ({(best['C2']/baseline['C2']-1)*100:+.0f}%)")
    print()

    print("DISTANCE TO TARGET:")
    if best['H_ratio'] < 0.6:
        print(f"  H:  {0.6/best['H_ratio']:.2f}× too LOW (need {0.6*TARGETS['H']:.2e} minimum)")
    elif best['H_ratio'] > 1.4:
        print(f"  H:  {best['H_ratio']/1.4:.2f}× too HIGH")
    else:
        print(f"  H:  ✓ WITHIN TARGET RANGE!")

    if best['CH_ratio'] < 0.6:
        print(f"  CH: {0.6/best['CH_ratio']:.2f}× too LOW")
    elif best['CH_ratio'] > 1.4:
        print(f"  CH: {best['CH_ratio']/1.4:.2f}× too HIGH (need {1.4*TARGETS['CH']:.2e} maximum)")
    else:
        print(f"  CH: ✓ WITHIN TARGET RANGE!")

    if best['C2_ratio'] < 0.6:
        print(f"  C2: {0.6/best['C2_ratio']:.2f}× too LOW (need {0.6*TARGETS['C2']:.2e} minimum)")
    elif best['C2_ratio'] > 1.4:
        print(f"  C2: {best['C2_ratio']/1.4:.2f}× too HIGH")
    else:
        print(f"  C2: ✓ WITHIN TARGET RANGE!")

    print()
    print("="*80)

    # Check if we're in target range for all
    all_in_range = (0.6 <= best['H_ratio'] <= 1.4 and
                    0.6 <= best['CH_ratio'] <= 1.4 and
                    0.6 <= best['C2_ratio'] <= 1.4)

    if all_in_range:
        print("🎉 SUCCESS! All species within target range!")
    else:
        print("Still need improvement - not all species in target range")
