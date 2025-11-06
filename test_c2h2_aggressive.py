#!/usr/bin/env python3
"""
Aggressive C2H2 optimization at NEW best conditions (with stick_H reduced)
Test all available C2H2 parameters:
- stick_C2H2 (already tested, but re-verify)
- loss_C2H2 (NOT tested yet!)
- loss_C2H2Star (NOT tested yet!)
- C2H2 production reactions (boosted before, but at different conditions)
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

# BEST conditions (with stick_H reduced)
TE_BEST = 1.25
E_BEST = 100
NE_BEST = 5.0e8

print("="*80)
print("AGGRESSIVE C2H2 OPTIMIZATION AT NEW BEST CONDITIONS")
print("="*80)
print()
print(f"Conditions: Te={TE_BEST} eV, E={E_BEST} V/cm, Ne={NE_BEST:.1e} cm⁻³")
print(f"            stick_H = 3.89e2 (reduced!)")
print()
print("C2H2 Parameter Ranges:")
print(f"  stick_C2H2:      {db['stick_C2H2_9_11'].min:.0f} to {db['stick_C2H2_9_11'].max:.0f}")
print(f"  loss_C2H2:       {db['loss_C2H2_11_19'].min:.0f} to {db['loss_C2H2_11_19'].max:.0f}")
print(f"  loss_C2H2Star:   {db['loss_C2H2Star_11_25'].min:.0f} to {db['loss_C2H2Star_11_25'].max:.0f}")
print("="*80)
print()

def run_simulation(rate_values_override=None):
    """Run simulation and return results."""
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply BEST chemical parameters
    params['rate_values']['stick_C2_9_9'] = 1.25e3
    params['rate_values']['loss_C2_11_3'] = 100
    params['rate_values']['C2_H_CH_C_cm3_7_6'] = 8.0e-11
    params['rate_values']['stick_H_9_1'] = 3.89e2  # NEW: reduced stick_H

    # Apply any overrides
    if rate_values_override:
        params['rate_values'].update(rate_values_override)

    # Apply BEST conditions
    params['Te'] = TE_BEST
    params['E_field'] = E_BEST
    params['ne'] = NE_BEST

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
        return None

    # Extract final densities
    y_final = sol.y[:, -1]
    H_idx = ode.species.index('H')
    CH_idx = ode.species.index('CH')
    C2_idx = ode.species.index('C2')
    C2H2_idx = ode.species.index('C2H2')

    H = y_final[H_idx]
    CH = y_final[CH_idx]
    C2 = y_final[C2_idx]
    C2H2 = y_final[C2H2_idx]

    # Calculate total ion density
    ni_total = 0
    for i_sp, species_name in enumerate(ode.species):
        if species_name.endswith('Plus'):
            ni_total += y_final[i_sp]
        elif species_name.endswith('Minus'):
            ni_total -= y_final[i_sp]

    ni_ne_ratio = ni_total / y_final[ode.species.index('e')]

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

    return {
        'H': H,
        'CH': CH,
        'C2': C2,
        'C2H2': C2H2,
        'ni_ne_ratio': ni_ne_ratio,
        'objective': objective,
    }

# Test configurations
configs = [
    {
        'name': 'Baseline (stick_H reduced)',
        'rate_values': {},
        'description': 'stick_H=3.89e2, default C2H2 parameters'
    },
    {
        'name': 'Minimize stick_C2H2',
        'rate_values': {
            'stick_C2H2_9_11': 500,  # Literature minimum
        },
        'description': 'stick_C2H2 = 500 (minimum)'
    },
    {
        'name': 'Minimize loss_C2H2',
        'rate_values': {
            'loss_C2H2_11_19': 1000,  # Already at minimum
        },
        'description': 'loss_C2H2 = 1000 (already at minimum)'
    },
    {
        'name': 'Minimize loss_C2H2Star',
        'rate_values': {
            'loss_C2H2Star_11_25': 100,  # Minimum
        },
        'description': 'loss_C2H2Star = 100 (minimum)'
    },
    {
        'name': 'Boost C2H2 production',
        'rate_values': {
            'CH_CH2_C2H2_H_cm3_7_7': 1.20e-10,      # max
            'CH2_CH2_C2H2_H2_cm3_7_15': 1.20e-10,   # max
            'CH3_CH_C2H2_H2_cm3_7_16': 1.20e-10,    # max
            'CH2_C_C2H2_cm3_7_17': 1.20e-10,        # max
        },
        'description': 'All C2H2 production reactions to max'
    },
    {
        'name': 'Reduce C2H2 → C3 consumption',
        'rate_values': {
            'CH_C2H2_C3H2_H_cm3_7_22': 8.0e-11,     # min
            'CH_C2H2_C3H_H2_cm3_7_27': 8.0e-11,     # min
            'CH_C2H2_C2H_CH2_cm3_7_29': 8.0e-11,    # min
        },
        'description': 'C2H2 → C3 consumption to minimum'
    },
    {
        'name': 'ALL C2H2 optimizations COMBINED',
        'rate_values': {
            # Wall losses to minimum
            'stick_C2H2_9_11': 500,
            'loss_C2H2_11_19': 1000,
            'loss_C2H2Star_11_25': 100,
            # Production to maximum
            'CH_CH2_C2H2_H_cm3_7_7': 1.20e-10,
            'CH2_CH2_C2H2_H2_cm3_7_15': 1.20e-10,
            'CH3_CH_C2H2_H2_cm3_7_16': 1.20e-10,
            'CH2_C_C2H2_cm3_7_17': 1.20e-10,
            # C3 consumption to minimum
            'CH_C2H2_C3H2_H_cm3_7_22': 8.0e-11,
            'CH_C2H2_C3H_H2_cm3_7_27': 8.0e-11,
            'CH_C2H2_C2H_CH2_cm3_7_29': 8.0e-11,
        },
        'description': 'ALL C2H2 parameters optimized simultaneously'
    },
]

results = []
baseline_result = None

for i, config in enumerate(configs):
    print(f"[{i+1}/{len(configs)}] Testing: {config['name']}")
    print(f"  {config['description']}")
    print("  ", end='', flush=True)

    result = run_simulation(config['rate_values'])

    if result is None:
        print("FAILED")
        continue

    if i == 0:
        baseline_result = result

    results.append({
        'name': config['name'],
        **result
    })

    ni_ne_ok = "✓" if 2.0 <= result['ni_ne_ratio'] <= 6.0 else " "

    # Calculate changes vs baseline
    if baseline_result:
        c2_change = (result['C2'] / baseline_result['C2'] - 1) * 100
        c2h2_change = (result['C2H2'] / baseline_result['C2H2'] - 1) * 100
        change_str = f" | Δ C2{c2_change:+.1f}% C2H2{c2h2_change:+.1f}%"
    else:
        change_str = ""

    print(f"Ni/Ne={result['ni_ne_ratio']:.2f}{ni_ne_ok} | f(x)={result['objective']:.1f} | "
          f"H={result['H']/TARGETS['H']:.2f}× | CH={result['CH']/TARGETS['CH']:.1f}× | "
          f"C2={result['C2']/TARGETS['C2']:.2f}× | C2H2={result['C2H2']:.2e}{change_str}")
    print()

print()
print("="*80)
print("COMPARISON OF RESULTS")
print("="*80)
print()

if results:
    print(f"{'Config':<40} | {'Ni/Ne':>7} | {'f(x)':>7} | {'C2':>7} | {'C2H2':>11} | {'ΔC2':>7}")
    print("-"*105)

    for r in results:
        ni_ne_ok = "✓" if 2.0 <= r['ni_ne_ratio'] <= 6.0 else " "
        c2_ratio = r['C2'] / TARGETS['C2']

        if baseline_result:
            c2_change = (r['C2'] / baseline_result['C2'] - 1) * 100
            change_str = f"{c2_change:+6.1f}%"
        else:
            change_str = "    -  "

        config_name = r['name'][:38]
        print(f"{config_name:<40} | {r['ni_ne_ratio']:>6.2f}{ni_ne_ok} | {r['objective']:>7.1f} | "
              f"{c2_ratio:>6.2f}× | {r['C2H2']:>11.2e} | {change_str}")

    print()
    print("="*80)
    print("BEST RESULTS:")
    print("="*80)

    # Best C2
    best_c2 = max(results, key=lambda x: x['C2'])
    print()
    print("BEST C2:")
    print(f"  {best_c2['name']}")
    print(f"  C2 = {best_c2['C2']:.2e} cm⁻³ ({best_c2['C2']/TARGETS['C2']:.2f}× target)")
    print(f"  C2H2 = {best_c2['C2H2']:.2e} cm⁻³")
    print(f"  H = {best_c2['H']:.2e} cm⁻³ ({best_c2['H']/TARGETS['H']:.2f}×)")
    print(f"  CH = {best_c2['CH']:.2e} cm⁻³ ({best_c2['CH']/TARGETS['CH']:.1f}×)")
    print(f"  Ni/Ne = {best_c2['ni_ne_ratio']:.2f}")
    print(f"  f(x) = {best_c2['objective']:.1f}")

    if baseline_result:
        print(f"  C2 improvement: {(best_c2['C2']/baseline_result['C2'] - 1)*100:+.1f}%")
        print(f"  C2H2 improvement: {(best_c2['C2H2']/baseline_result['C2H2'] - 1)*100:+.1f}%")

    # Best objective
    best_obj = min(results, key=lambda x: x['objective'])
    print()
    print("BEST OBJECTIVE:")
    print(f"  {best_obj['name']}")
    print(f"  f(x) = {best_obj['objective']:.1f}")
    print(f"  H={best_obj['H']/TARGETS['H']:.2f}× | CH={best_obj['CH']/TARGETS['CH']:.1f}× | "
          f"C2={best_obj['C2']/TARGETS['C2']:.2f}×")

print()
print("="*80)

# Save results
with open('c2h2_aggressive_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("Results saved to: c2h2_aggressive_results.json")
