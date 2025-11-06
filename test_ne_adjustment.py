#!/usr/bin/env python3
"""
Test adjusting Ne to ~3e9 to achieve proper sheath-edge charge balance
Expected: Ni ~ 2-6× Ne (for sheath edge region)
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

db = get_complete_rate_database()

# Test these Ne values around 3e9
ne_values = [
    8.44e8,   # Current checkpoint value
    1.5e9,
    2.0e9,
    2.5e9,
    3.0e9,    # Target value (user's educated guess)
    3.5e9,
    4.0e9,
    5.0e9,
]

print("="*80)
print("Ne ADJUSTMENT TEST - Achieve proper sheath-edge charge balance")
print("="*80)
print(f"\nGoal: Ni ~ 2-6× Ne (expected for sheath edge)")
print(f"Current: Ne = 8.44e8, Ni = 7.8e9 → Ni/Ne = 9.2× (too high!)")
print(f"Target:  Ne ~ 3e9, Ni ~ 7.8e9 → Ni/Ne ~ 2.6× (perfect!)")
print()
print(f"Testing Ne values: {[f'{v:.2e}' for v in ne_values]}")
print("="*80)
print()

results = []

for ne in ne_values:
    print(f"Testing Ne = {ne:.2e}...", end='', flush=True)

    # Start from BEST parameters
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply BEST parameters from all previous tests
    params['rate_values']['stick_C2_9_9'] = 1.25e3      # lit min
    params['rate_values']['loss_C2_11_3'] = 100         # near min
    params['rate_values']['C2_H_CH_C_cm3_7_6'] = 8.0e-11  # lit min
    params['rate_values']['stick_C2H2_9_11'] = 472      # 50% reduction

    # Adjust Ne
    params['ne'] = ne

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
        y0[ode.species.index('e')] = ne  # Use adjusted Ne
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
        ne_final = y_final[ode.species.index('e')]

        # Calculate total ion density
        ni_total = 0
        for i, species_name in enumerate(ode.species):
            if species_name.endswith('Plus'):
                ni_total += y_final[i]
            elif species_name.endswith('Minus'):
                ni_total -= y_final[i]

        # Charge balance
        ni_ne_ratio = ni_total / ne_final

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
            'ne': ne,
            'ne_final': ne_final,
            'ni': ni_total,
            'ni_ne_ratio': ni_ne_ratio,
            'H': H,
            'CH': CH,
            'C2': C2,
            'C2H2': C2H2,
            'objective': objective
        })

        # Check if in expected range
        in_range = "✓" if 2.0 <= ni_ne_ratio <= 6.0 else "✗"

        print(f" {in_range} Ni/Ne={ni_ne_ratio:.2f} | f(x)={objective:.1f} | H={H:.2e} ({H/TARGETS['H']:.2f}×) | C2={C2:.2e} ({C2/TARGETS['C2']:.2f}×) | CH={CH:.2e} ({CH/TARGETS['CH']:.2f}×)")

    except Exception as e:
        import traceback
        print(f" ERROR: {e}")
        traceback.print_exc()
        continue

print()
print("="*80)
print("SUMMARY")
print("="*80)
print(f"{'Ne':>12} | {'Ni':>12} | {'Ni/Ne':>8} | {'Range':>6} | {'f(x)':>10} | {'H ratio':>8} | {'C2 ratio':>8} | {'CH ratio':>8}")
print("-"*80)

for r in results:
    in_range = "✓" if 2.0 <= r['ni_ne_ratio'] <= 6.0 else "✗"
    marker = ""
    if abs(r['ne'] - 8.44e8) < 1e7:
        marker = " ← checkpoint"
    elif abs(r['ne'] - 3.0e9) < 1e8:
        marker = " ← target"

    print(f"{r['ne']:>12.2e} | {r['ni']:>12.2e} | {r['ni_ne_ratio']:>8.2f} | {in_range:>6} | {r['objective']:>10.1f} | "
          f"{r['H']/TARGETS['H']:>8.2f} | {r['C2']/TARGETS['C2']:>8.2f} | {r['CH']/TARGETS['CH']:>8.2f}{marker}")

print()
print("Target Ni/Ne range: [2.0, 6.0] for sheath edge")
print()

# Find best in target range
in_range_results = [r for r in results if 2.0 <= r['ni_ne_ratio'] <= 6.0]
if in_range_results:
    best = min(in_range_results, key=lambda x: x['objective'])
    print(f"BEST with proper charge balance (Ni/Ne in [2, 6]):")
    print(f"  Ne:     {best['ne']:.2e} cm⁻³")
    print(f"  Ni:     {best['ni']:.2e} cm⁻³")
    print(f"  Ni/Ne:  {best['ni_ne_ratio']:.2f} ✓")
    print(f"  f(x):   {best['objective']:.1f}")
    print(f"  H:      {best['H']:.2e} ({best['H']/TARGETS['H']:.2f}× target)")
    print(f"  C2:     {best['C2']:.2e} ({best['C2']/TARGETS['C2']:.2f}× target)")
    print(f"  CH:     {best['CH']:.2e} ({best['CH']/TARGETS['CH']:.2f}× target)")
    print(f"  C2H2:   {best['C2H2']:.2e}")
else:
    print("No results found in target Ni/Ne range [2, 6]")
