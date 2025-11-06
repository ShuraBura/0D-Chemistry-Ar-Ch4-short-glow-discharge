#!/usr/bin/env python3
"""
Manual sweep of loss_C2 values to see impact on C2 density
Starting from best known parameters (checkpoint_f3407.json)
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from odefun import PlasmaODE
from define_rates_tunable import define_rates_tunable
from build_reactions import build_reactions

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

# Test these loss_C2 values
loss_c2_values = [100, 200, 500, 800, 1000, 1200, 1500, 2000]

print("="*80)
print("MANUAL SWEEP: loss_C2 parameter")
print("="*80)
print(f"\nBest known result (f(x)={best['objective']:.1f}):")
print(f"  loss_C2 = {params_base['rate_values']['loss_C2_11_3']:.1f}")
print(f"  H:  {best['densities']['H']:.2e} ({best['densities']['H']/TARGETS['H']:.2f}× target)")
print(f"  C2: {best['densities']['C2']:.2e} ({best['densities']['C2']/TARGETS['C2']:.2f}× target)")
print(f"  CH: {best['densities']['CH']:.2e} ({best['densities']['CH']/TARGETS['CH']:.2f}× target)")
print()
print("Testing loss_C2 values:", loss_c2_values)
print("="*80)
print()

results = []

for loss_c2 in loss_c2_values:
    print(f"Testing loss_C2 = {loss_c2:.0f}...", end='', flush=True)

    # Create params with modified loss_C2
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()
    params['rate_values']['loss_C2_11_3'] = loss_c2

    try:
        # Get rate dictionary with modified loss_C2
        k_dict = define_rates_tunable(params)
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
            'loss_c2': loss_c2,
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
print(f"{'loss_C2':>10} | {'f(x)':>10} | {'H ratio':>8} | {'C2 ratio':>8} | {'CH ratio':>8}")
print("-"*80)

for r in results:
    print(f"{r['loss_c2']:>10.0f} | {r['objective']:>10.1f} | "
          f"{r['H']/TARGETS['H']:>8.2f} | {r['C2']/TARGETS['C2']:>8.2f} | {r['CH']/TARGETS['CH']:>8.2f}")

# Find best
if results:
    best_result = min(results, key=lambda x: x['objective'])
    print()
    print(f"BEST: loss_C2 = {best_result['loss_c2']:.0f}, f(x) = {best_result['objective']:.1f}")
    print(f"  H:  {best_result['H']:.2e} ({best_result['H']/TARGETS['H']:.2f}× target)")
    print(f"  C2: {best_result['C2']:.2e} ({best_result['C2']/TARGETS['C2']:.2f}× target)")
    print(f"  CH: {best_result['CH']:.2e} ({best_result['CH']/TARGETS['CH']:.2f}× target)")
