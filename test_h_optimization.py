#!/usr/bin/env python3
"""
Test reducing H consumption pathways to boost H density
1. Reduce H + CH4 → CH3 + H2 reaction to minimum
2. Check H wall sticking (stick_H) - may already be at minimum
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

# BEST conditions
TE_BEST = 1.25
E_BEST = 100
NE_BEST = 5.0e8

print("="*80)
print("H OPTIMIZATION - REDUCE H CONSUMPTION")
print("="*80)
print()
print(f"Conditions: Te={TE_BEST} eV, E={E_BEST} V/cm, Ne={NE_BEST:.1e} cm⁻³")
print()
print("Literature ranges:")
print(f"  H_CH4_CH3_H2: {db['H_CH4_CH3_H2_cm3_7_25'].min:.2e} to {db['H_CH4_CH3_H2_cm3_7_25'].max:.2e}")
print(f"  stick_H:      {db['stick_H_9_1'].min:.2e} to {db['stick_H_9_1'].max:.2e}")
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
    params['rate_values']['stick_C2H2_9_11'] = 472

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

    H = y_final[H_idx]
    CH = y_final[CH_idx]
    C2 = y_final[C2_idx]

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
        'ni_ne_ratio': ni_ne_ratio,
        'objective': objective,
        'k_dict': k_dict,
    }

# Test configurations
configs = [
    {
        'name': 'Baseline (current best)',
        'rate_values': {},
        'description': 'Current best parameters'
    },
    {
        'name': 'Reduce H + CH4 to minimum',
        'rate_values': {
            'H_CH4_CH3_H2_cm3_7_25': 4.0e-12,  # Min (from 6.0e-12)
        },
        'description': 'H_CH4_CH3_H2 = 4.0e-12 (literature minimum)'
    },
    {
        'name': 'Reduce stick_H to minimum',
        'rate_values': {
            'stick_H_9_1': 3.89e2,  # Already at minimum
        },
        'description': 'stick_H = 3.89e2 (already at minimum)'
    },
    {
        'name': 'COMBINED: Both at minimum',
        'rate_values': {
            'H_CH4_CH3_H2_cm3_7_25': 4.0e-12,
            'stick_H_9_1': 3.89e2,
        },
        'description': 'Both H + CH4 and stick_H at minimum'
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
        h_change = (result['H'] / baseline_result['H'] - 1) * 100
        ch_change = (result['CH'] / baseline_result['CH'] - 1) * 100
        c2_change = (result['C2'] / baseline_result['C2'] - 1) * 100
        change_str = f" | Δ: H{h_change:+.1f}% CH{ch_change:+.1f}% C2{c2_change:+.1f}%"
    else:
        change_str = ""

    print(f"Ni/Ne={result['ni_ne_ratio']:.2f}{ni_ne_ok} | f(x)={result['objective']:.1f} | "
          f"H={result['H']/TARGETS['H']:.2f}× | CH={result['CH']/TARGETS['CH']:.1f}× | "
          f"C2={result['C2']/TARGETS['C2']:.2f}×{change_str}")
    print()

print()
print("="*80)
print("COMPARISON OF RESULTS")
print("="*80)
print()

if results and baseline_result:
    print(f"{'Config':<35} | {'Ni/Ne':>7} | {'f(x)':>7} | {'H':>7} | {'CH':>7} | {'C2':>7}")
    print("-"*95)

    for r in results:
        ni_ne_ok = "✓" if 2.0 <= r['ni_ne_ratio'] <= 6.0 else " "
        h_ratio = r['H'] / TARGETS['H']
        ch_ratio = r['CH'] / TARGETS['CH']
        c2_ratio = r['C2'] / TARGETS['C2']

        config_name = r['name'][:33]
        print(f"{config_name:<35} | {r['ni_ne_ratio']:>6.2f}{ni_ne_ok} | {r['objective']:>7.1f} | "
              f"{h_ratio:>6.2f}× | {ch_ratio:>6.1f}× | {c2_ratio:>6.2f}×")

    print()
    print("="*80)
    print("BEST RESULT:")
    print("="*80)

    best = min(results, key=lambda x: x['objective'])
    print(f"  {best['name']}")
    print(f"  f(x) = {best['objective']:.1f}")
    print(f"  Ni/Ne = {best['ni_ne_ratio']:.2f}")
    print(f"  H  = {best['H']:.2e} cm⁻³ ({best['H']/TARGETS['H']:.2f}× target)")
    print(f"  CH = {best['CH']:.2e} cm⁻³ ({best['CH']/TARGETS['CH']:.1f}× target)")
    print(f"  C2 = {best['C2']:.2e} cm⁻³ ({best['C2']/TARGETS['C2']:.2f}× target)")
    print()

    if best != results[0]:
        print(f"Improvement vs baseline:")
        print(f"  H:  {(best['H']/baseline_result['H'] - 1)*100:+.1f}%")
        print(f"  CH: {(best['CH']/baseline_result['CH'] - 1)*100:+.1f}%")
        print(f"  C2: {(best['C2']/baseline_result['C2'] - 1)*100:+.1f}%")
        print(f"  f(x): {(best['objective']/baseline_result['objective'] - 1)*100:+.1f}%")

print()
print("="*80)

# Save results
with open('h_optimization_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("Results saved to: h_optimization_results.json")
