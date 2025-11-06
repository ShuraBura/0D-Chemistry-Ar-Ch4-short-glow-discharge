#!/usr/bin/env python3
"""
Ne sweep at SWEET SPOT conditions (Te=1.25 eV, E=100 V/cm)
Goal: See how Ne affects densities and Ni/Ne ratio at optimal Te/E
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

# FIXED sweet spot conditions
TE_SWEET = 1.25  # eV
E_SWEET = 100    # V/cm

# Ne values to test - vary around 3e9
ne_values = [
    5e8, 8e8, 1e9, 1.5e9, 2e9, 2.5e9,
    3e9, 3.5e9, 4e9, 5e9, 6e9, 8e9, 1e10
]

print("="*80)
print("Ne SWEEP AT SWEET SPOT CONDITIONS")
print("="*80)
print()
print(f"FIXED: Te = {TE_SWEET} eV, E = {E_SWEET} V/cm")
print(f"VARYING: Ne from {ne_values[0]:.1e} to {ne_values[-1]:.1e} cm⁻³")
print()
print("Goal: Find Ne that optimizes densities while maintaining good Ni/Ne ratio")
print("="*80)
print()

results = []

for i, ne in enumerate(ne_values):
    print(f"[{i+1}/{len(ne_values)}] Testing Ne={ne:.2e} cm⁻³...", end='', flush=True)

    # Start from BEST parameters
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply BEST chemical parameters
    params['rate_values']['stick_C2_9_9'] = 1.25e3
    params['rate_values']['loss_C2_11_3'] = 100
    params['rate_values']['C2_H_CH_C_cm3_7_6'] = 8.0e-11
    params['rate_values']['stick_C2H2_9_11'] = 472

    # Apply SWEET SPOT conditions
    params['Te'] = TE_SWEET
    params['E_field'] = E_SWEET
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
        ne_final = y_final[ode.species.index('e')]

        # Calculate total ion density
        ni_total = 0
        for i_sp, species_name in enumerate(ode.species):
            if species_name.endswith('Plus'):
                ni_total += y_final[i_sp]
            elif species_name.endswith('Minus'):
                ni_total -= y_final[i_sp]

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
            'objective': objective,
            'H_ratio': H / TARGETS['H'],
            'CH_ratio': CH / TARGETS['CH'],
            'C2_ratio': C2 / TARGETS['C2'],
        })

        # Check if in target ranges
        ni_ne_ok = "✓" if 2.0 <= ni_ne_ratio <= 6.0 else " "

        print(f" Ni/Ne={ni_ne_ratio:5.2f}{ni_ne_ok} | f(x)={objective:6.0f} | "
              f"H={H/TARGETS['H']:.2f}× | CH={CH/TARGETS['CH']:.1f}× | C2={C2/TARGETS['C2']:.2f}×")

    except Exception as e:
        print(f" ERROR: {e}")
        continue

print()
print("="*80)
print("ANALYSIS OF Ne IMPACT")
print("="*80)
print()

if results:
    # Best overall objective
    best_obj = min(results, key=lambda x: x['objective'])
    print("BEST OBJECTIVE FUNCTION:")
    print(f"  Ne={best_obj['ne']:.2e} cm⁻³")
    print(f"  f(x)={best_obj['objective']:.1f}")
    print(f"  Ni/Ne={best_obj['ni_ne_ratio']:.2f} | H={best_obj['H_ratio']:.2f}× | "
          f"CH={best_obj['CH_ratio']:.1f}× | C2={best_obj['C2_ratio']:.2f}×")
    print()

    # Best Ni/Ne ratio (closest to target 2-6×)
    def ni_ne_score(r):
        ratio = r['ni_ne_ratio']
        if 2.0 <= ratio <= 6.0:
            return abs(ratio - 3.5)  # Prefer middle of range
        elif ratio < 2.0:
            return 2.0 - ratio + 10
        else:
            return ratio - 6.0 + 10

    best_ni_ne = min(results, key=ni_ne_score)
    print("BEST Ni/Ne RATIO (closest to target 2-6×):")
    print(f"  Ne={best_ni_ne['ne']:.2e} cm⁻³")
    print(f"  Ni/Ne={best_ni_ne['ni_ne_ratio']:.2f}")
    print(f"  f(x)={best_ni_ne['objective']:.1f} | H={best_ni_ne['H_ratio']:.2f}× | "
          f"CH={best_ni_ne['CH_ratio']:.1f}× | C2={best_ni_ne['C2_ratio']:.2f}×")
    print()

    # Best C2 (closest to target)
    best_c2 = min(results, key=lambda x: abs(x['C2_ratio'] - 1.0))
    print("BEST C2 (closest to target):")
    print(f"  Ne={best_c2['ne']:.2e} cm⁻³")
    print(f"  C2={best_c2['C2_ratio']:.2f}×")
    print(f"  Ni/Ne={best_c2['ni_ne_ratio']:.2f} | f(x)={best_c2['objective']:.1f} | "
          f"H={best_c2['H_ratio']:.2f}× | CH={best_c2['CH_ratio']:.1f}×")
    print()

    # Best CH (closest to target)
    best_ch = min(results, key=lambda x: abs(x['CH_ratio'] - 1.0))
    print("BEST CH (closest to target):")
    print(f"  Ne={best_ch['ne']:.2e} cm⁻³")
    print(f"  CH={best_ch['CH_ratio']:.2f}×")
    print(f"  Ni/Ne={best_ch['ni_ne_ratio']:.2f} | f(x)={best_ch['objective']:.1f} | "
          f"H={best_ch['H_ratio']:.2f}× | C2={best_ch['C2_ratio']:.2f}×")
    print()

print()
print("="*80)
print("FULL RESULTS TABLE")
print("="*80)
print(f"{'Ne':>10} | {'Ni/Ne':>7} | {'f(x)':>8} | {'H':>6} | {'CH':>7} | {'C2':>6}")
print("-"*80)

for r in results:
    ni_ne_ok = "✓" if 2.0 <= r['ni_ne_ratio'] <= 6.0 else " "
    print(f"{r['ne']:>10.2e} | {r['ni_ne_ratio']:>6.2f}{ni_ne_ok} | {r['objective']:>8.0f} | "
          f"{r['H_ratio']:>6.2f} | {r['CH_ratio']:>6.1f}× | {r['C2_ratio']:>6.2f}")

print()
print("="*80)

# Save detailed results
with open('ne_sweep_sweet_spot_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Detailed results saved to: ne_sweep_sweet_spot_results.json")
