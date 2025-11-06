#!/usr/bin/env python3
"""
Test maximizing CH consumption reactions
Goal: Reduce CH density by increasing consumption rates
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
TE_BEST = 1.20
E_BEST = 75
NE_BEST = 1.0e8

print("="*80)
print("MAXIMIZE CH CONSUMPTION TEST")
print("="*80)
print()
print(f"Conditions: Te={TE_BEST} eV, E={E_BEST} V/cm, Ne={NE_BEST:.1e} cm⁻³")
print()
print("Testing strategies to increase CH consumption:")
print("  1. Maximize loss_CH (volumetric loss) - +229% potential")
print("  2. Maximize stick_CH (wall sticking) - +10% potential")
print("  3. Maximize CH reaction rates - various increases")
print("  4. Combined ALL consumption maximizations")
print("="*80)
print()

def run_simulation(rate_values_override=None):
    """Run simulation and return results."""
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply ALL CURRENT OPTIMIZATIONS (baseline)
    params['rate_values']['stick_H_9_1'] = 3.89e2
    params['rate_values']['stick_C2_9_9'] = 1.25e3
    params['rate_values']['loss_C2_11_3'] = 100
    params['rate_values']['C2_H_CH_C_cm3_7_6'] = 8.0e-11
    params['rate_values']['stick_C2H2_9_11'] = 500
    params['rate_values']['loss_C2H2_11_19'] = 1000
    params['rate_values']['loss_C2H2Star_11_25'] = 100
    params['rate_values']['CH_CH2_C2H2_H_cm3_7_7'] = 1.20e-10
    params['rate_values']['CH2_CH2_C2H2_H2_cm3_7_15'] = 1.20e-10
    params['rate_values']['CH3_CH_C2H2_H2_cm3_7_16'] = 1.20e-10
    params['rate_values']['CH2_C_C2H2_cm3_7_17'] = 1.20e-10
    params['rate_values']['CH_C2H2_C3H2_H_cm3_7_22'] = 8.0e-11
    params['rate_values']['CH_C2H2_C3H_H2_cm3_7_27'] = 8.0e-11
    params['rate_values']['CH_C2H2_C2H_CH2_cm3_7_29'] = 8.0e-11
    params['rate_values']['e_CH4_CH_H_H2_cm3_1_11'] = 2.0e-11
    params['rate_values']['e_CH4_CH_H2_H_vib_cm3_1_3'] = 2.0e-11

    # Apply any CH consumption maximizations
    if rate_values_override:
        params['rate_values'].update(rate_values_override)

    params['Te'] = TE_BEST
    params['E_field'] = E_BEST
    params['ne'] = NE_BEST

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
            'H_ratio': H / TARGETS['H'],
            'CH_ratio': CH / TARGETS['CH'],
            'C2_ratio': C2 / TARGETS['C2'],
        }

    except Exception as e:
        return None

# Test configurations
configs = [
    {
        'name': 'Baseline (current best)',
        'rate_values': {},
        'description': 'Current best with CH production minimized'
    },
    {
        'name': 'Maximize loss_CH',
        'rate_values': {
            'loss_CH_11_9': 1.00e4,  # Maximum (from 3.04e3, +229%!)
        },
        'description': 'loss_CH to maximum (1.00e4)'
    },
    {
        'name': 'Maximize stick_CH',
        'rate_values': {
            'stick_CH_9_3': 6.25e3,  # Maximum (from 5.69e3, +10%)
        },
        'description': 'stick_CH to maximum (6.25e3)'
    },
    {
        'name': 'Maximize CH reaction rates',
        'rate_values': {
            'CH_CH4_CH2_CH3_cm3_7_39': 1.20e-11,  # Max (+20%)
            'CH_H_CH2_cm3_7_21': 1.20e-10,        # Max (+20%)
            'CH_CH3_C2H3_H_cm3_7_10': 1.20e-10,   # Max (+50%)
            'CH_CH3_C2H2_H2_cm3_7_23': 1.20e-10,  # Max (+20%)
            'CH_H2_CH2_H_cm3_7_30': 1.20e-11,     # Max (+20%)
            'CH_H2_CH2_H_cm3_7_51': 1.20e-11,     # Max (+20%)
            'C2H2_CH_C3_H2_H_cm3_7_56': 1.20e-10, # Max (+20%)
            'CH_C2H2_C3H2_H_cm3_7_22': 1.20e-10,  # Max (+50%, was minimized!)
        },
        'description': 'All tunable CH reaction rates to maximum'
    },
    {
        'name': 'COMBINED: All CH consumption maximized',
        'rate_values': {
            # Loss terms
            'loss_CH_11_9': 1.00e4,
            'stick_CH_9_3': 6.25e3,
            # Reaction rates
            'CH_CH4_CH2_CH3_cm3_7_39': 1.20e-11,
            'CH_H_CH2_cm3_7_21': 1.20e-10,
            'CH_CH3_C2H3_H_cm3_7_10': 1.20e-10,
            'CH_CH3_C2H2_H2_cm3_7_23': 1.20e-10,
            'CH_H2_CH2_H_cm3_7_30': 1.20e-11,
            'CH_H2_CH2_H_cm3_7_51': 1.20e-11,
            'C2H2_CH_C3_H2_H_cm3_7_56': 1.20e-10,
            'CH_C2H2_C3H2_H_cm3_7_22': 1.20e-10,
        },
        'description': 'ALL CH consumption parameters maximized'
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
        ch_change = (result['CH'] / baseline_result['CH'] - 1) * 100
        h_change = (result['H'] / baseline_result['H'] - 1) * 100
        c2_change = (result['C2'] / baseline_result['C2'] - 1) * 100
        change_str = f" | Δ: CH{ch_change:+.1f}% H{h_change:+.1f}% C2{c2_change:+.1f}%"
    else:
        change_str = ""

    print(f"Ni/Ne={result['ni_ne_ratio']:.2f}{ni_ne_ok} | f(x)={result['objective']:.1f} | "
          f"H={result['H_ratio']:.2f}× | CH={result['CH_ratio']:.2f}× | "
          f"C2={result['C2_ratio']:.2f}×{change_str}")
    print()

print()
print("="*80)
print("COMPARISON OF RESULTS")
print("="*80)
print()

if results and baseline_result:
    print(f"{'Config':<40} | {'Ni/Ne':>7} | {'f(x)':>7} | {'H':>6} | {'CH':>6} | {'C2':>6}")
    print("-"*95)

    for r in results:
        ni_ne_ok = "✓" if 2.0 <= r['ni_ne_ratio'] <= 6.0 else " "
        h_ratio = r['H'] / TARGETS['H']
        ch_ratio = r['CH'] / TARGETS['CH']
        c2_ratio = r['C2'] / TARGETS['C2']

        config_name = r['name'][:38]
        print(f"{config_name:<40} | {r['ni_ne_ratio']:>6.2f}{ni_ne_ok} | {r['objective']:>7.1f} | "
              f"{h_ratio:>6.2f} | {ch_ratio:>6.2f} | {c2_ratio:>6.2f}")

    print()
    print("="*80)
    print("BEST RESULTS:")
    print("="*80)

    # Best CH (lowest)
    best_ch = min(results, key=lambda x: x['CH'])
    print()
    print("BEST (LOWEST) CH:")
    print(f"  {best_ch['name']}")
    print(f"  CH = {best_ch['CH']:.2e} cm⁻³ ({best_ch['CH_ratio']:.2f}×)")
    print(f"  H = {best_ch['H']:.2e} cm⁻³ ({best_ch['H_ratio']:.2f}×)")
    print(f"  C2 = {best_ch['C2']:.2e} cm⁻³ ({best_ch['C2_ratio']:.2f}×)")
    print(f"  Ni/Ne = {best_ch['ni_ne_ratio']:.2f}")
    print(f"  f(x) = {best_ch['objective']:.1f}")

    if baseline_result:
        ch_reduction = (1 - best_ch['CH'] / baseline_result['CH']) * 100
        print(f"  CH reduction: {ch_reduction:.1f}%")

    print()

    # Best objective
    best_obj = min(results, key=lambda x: x['objective'])
    print("BEST OBJECTIVE:")
    print(f"  {best_obj['name']}")
    print(f"  f(x) = {best_obj['objective']:.1f}")
    print(f"  H={best_obj['H_ratio']:.2f}× | CH={best_obj['CH_ratio']:.2f}× | C2={best_obj['C2_ratio']:.2f}×")
    print(f"  Ni/Ne = {best_obj['ni_ne_ratio']:.2f}")

    # Check if any results achieve ALL targets
    print()
    perfect_results = [r for r in results if
                      0.6 <= r['H_ratio'] <= 1.4 and
                      0.6 <= r['CH_ratio'] <= 1.4 and
                      0.6 <= r['C2_ratio'] <= 1.4 and
                      2.0 <= r['ni_ne_ratio'] <= 6.0]

    if perfect_results:
        print("🎯🎯🎯 PERFECT RESULTS (ALL 4 TARGETS IN RANGE!):")
        print("-"*80)
        for r in sorted(perfect_results, key=lambda x: x['objective']):
            print(f"  {r['name']}")
            print(f"    H={r['H_ratio']:.2f}× | CH={r['CH_ratio']:.2f}× | C2={r['C2_ratio']:.2f}× | Ni/Ne={r['ni_ne_ratio']:.2f}")
            print(f"    f(x)={r['objective']:.1f}")
            print()
    else:
        print("No results with ALL 4 targets in range yet.")

print()
print("="*80)

# Save results
with open('ch_consumption_maximized_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("Results saved to: ch_consumption_maximized_results.json")
