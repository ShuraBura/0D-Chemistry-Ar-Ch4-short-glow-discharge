#!/usr/bin/env python3
"""
Sweep to boost C2H2 production by reducing C2H2 losses
Strategy: Reduce stick_C2H2 and loss_C2H2 to increase C2H2 density
         More C2H2 → More C2 production (92% of C2 comes from C2H2 + H)
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

# Get literature ranges
db = get_complete_rate_database()
stick_c2h2_baseline = params_base['rate_values']['stick_C2H2_9_11']  # 944
loss_c2h2_baseline = params_base['rate_values']['loss_C2H2_11_19']  # (need to check)

# Test these reduction factors for C2H2 losses
# Lower factor = less C2H2 loss = more C2H2 available = more C2 production
c2h2_loss_factors = [
    0.2,   # 80% reduction
    0.3,   # 70% reduction
    0.4,   # 60% reduction
    0.5,   # 50% reduction
    0.6,   # 40% reduction
    0.7,   # 30% reduction
    0.8,   # 20% reduction
    1.0,   # Baseline (no reduction)
]

print("="*80)
print("C2H2 BOOST SWEEP - Reduce C2H2 losses to increase C2 production")
print("="*80)
print(f"\nStrategy: 92% of C2 production comes from C2H2 + H → C2 + H2")
print(f"By reducing C2H2 losses, we increase C2H2 density → more C2 production")
print()
print(f"Best known result (f(x)={best['objective']:.1f}):")
print(f"  H:    {best['densities']['H']:.2e} ({best['densities']['H']/TARGETS['H']:.2f}× target)")
print(f"  C2:   {best['densities']['C2']:.2e} ({best['densities']['C2']/TARGETS['C2']:.2f}× target)")
print(f"  CH:   {best['densities']['CH']:.2e} ({best['densities']['CH']/TARGETS['CH']:.2f}× target)")
print(f"  C2H2: {best['densities']['C2H2']:.2e}")
print()
print(f"Baseline C2H2 losses:")
print(f"  stick_C2H2: {stick_c2h2_baseline:.2e}")
print()
print(f"Testing reduction factors: {c2h2_loss_factors}")
print(f"  (multiply both stick_C2H2 and loss_C2H2 by factor)")
print("="*80)
print()

results = []

for factor in c2h2_loss_factors:
    print(f"Testing C2H2 loss factor = {factor:.2f}...", end='', flush=True)

    # Create params with reduced C2H2 losses
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Reduce both C2H2 loss mechanisms
    params['rate_values']['stick_C2H2_9_11'] *= factor
    if 'loss_C2H2_11_19' in params['rate_values']:
        params['rate_values']['loss_C2H2_11_19'] *= factor

    try:
        # Get rate dictionary with modified C2H2 losses
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
            'factor': factor,
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
print(f"{'Factor':>8} | {'f(x)':>10} | {'H ratio':>8} | {'C2 ratio':>8} | {'C2H2':>12} | {'CH ratio':>8}")
print("-"*80)

for r in results:
    marker = " ← baseline" if abs(r['factor'] - 1.0) < 1e-6 else ""
    print(f"{r['factor']:>8.2f} | {r['objective']:>10.1f} | "
          f"{r['H']/TARGETS['H']:>8.2f} | {r['C2']/TARGETS['C2']:>8.2f} | "
          f"{r['C2H2']:>12.2e} | {r['CH']/TARGETS['CH']:>8.2f}{marker}")

# Find best
if results:
    best_result = min(results, key=lambda x: x['objective'])
    print()
    print(f"BEST: Factor = {best_result['factor']:.2f}, f(x) = {best_result['objective']:.1f}")
    print(f"  H:    {best_result['H']:.2e} ({best_result['H']/TARGETS['H']:.2f}× target)")
    print(f"  C2:   {best_result['C2']:.2e} ({best_result['C2']/TARGETS['C2']:.2f}× target)")
    print(f"  C2H2: {best_result['C2H2']:.2e}")
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
        print(f"Improvement vs baseline:")
        print(f"  Objective: {improvement:.2f}× better")
        print(f"  C2 density: {c2_improvement:.2f}× higher")
        print(f"  C2H2 density: {c2h2_improvement:.2f}× higher")
