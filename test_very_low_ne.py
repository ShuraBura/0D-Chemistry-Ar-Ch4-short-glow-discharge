#!/usr/bin/env python3
"""
Very low Ne sweep - Test below 3e8 to find minimum CH
Goal: Push CH as close to target as possible
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

# BEST Te/E from previous sweeps
TE_BEST = 1.20
E_BEST = 75

# Very low Ne values to test
ne_values = [
    1e8, 1.5e8, 2e8, 2.5e8, 3e8, 3.5e8, 4e8, 5e8
]

print("="*80)
print("VERY LOW Ne SWEEP - Pushing CH to minimum")
print("="*80)
print()
print(f"Fixed: Te={TE_BEST} eV, E={E_BEST} V/cm")
print(f"Testing VERY LOW Ne: {ne_values[0]:.1e} to {ne_values[-1]:.1e} cm⁻³")
print()
print("Goal: Find lowest achievable CH while maintaining H, C2, and Ni/Ne")
print("="*80)
print()

def run_simulation(ne):
    """Run simulation and return results."""
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply ALL OPTIMIZATIONS
    params['rate_values']['stick_H_9_1'] = 3.89e2
    params['rate_values']['stick_C2_9_9'] = 1.25e3
    params['rate_values']['loss_C2_11_3'] = 100
    params['rate_values']['C2_H_CH_C_cm3_7_6'] = 8.0e-11
    # C2H2 optimizations
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
    # CH production minimized
    params['rate_values']['e_CH4_CH_H_H2_cm3_1_11'] = 2.0e-11
    params['rate_values']['e_CH4_CH_H2_H_vib_cm3_1_3'] = 2.0e-11

    # Apply Te, E, Ne
    params['Te'] = TE_BEST
    params['E_field'] = E_BEST
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
            'ne': ne,
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

results = []

for i, ne in enumerate(ne_values):
    print(f"[{i+1}/{len(ne_values)}] Testing Ne={ne:.2e} cm⁻³...", end='', flush=True)

    result = run_simulation(ne)

    if result is None:
        print(" FAILED")
        continue

    results.append(result)

    ni_ne_ok = "✓" if 2.0 <= result['ni_ne_ratio'] <= 6.0 else " "

    print(f" Ni/Ne={result['ni_ne_ratio']:5.2f}{ni_ne_ok} | f(x)={result['objective']:6.1f} | "
          f"H={result['H_ratio']:.2f}× | CH={result['CH_ratio']:.2f}× | C2={result['C2_ratio']:.2f}×")

print()
print("="*80)
print("ANALYSIS OF RESULTS")
print("="*80)
print()

if results:
    # Best objective
    best_obj = min(results, key=lambda x: x['objective'])
    print("BEST OBJECTIVE:")
    print(f"  Ne={best_obj['ne']:.2e} cm⁻³")
    print(f"  f(x)={best_obj['objective']:.1f} | Ni/Ne={best_obj['ni_ne_ratio']:.2f}")
    print(f"  H={best_obj['H']:.2e} cm⁻³ ({best_obj['H_ratio']:.2f}×)")
    print(f"  CH={best_obj['CH']:.2e} cm⁻³ ({best_obj['CH_ratio']:.2f}×)")
    print(f"  C2={best_obj['C2']:.2e} cm⁻³ ({best_obj['C2_ratio']:.2f}×)")
    print()

    # Best CH (closest to 1.0)
    best_ch = min(results, key=lambda x: abs(x['CH_ratio'] - 1.0))
    print("BEST CH (closest to target 1.0×):")
    print(f"  Ne={best_ch['ne']:.2e} cm⁻³")
    print(f"  CH={best_ch['CH']:.2e} cm⁻³ ({best_ch['CH_ratio']:.2f}×)")
    print(f"  H={best_ch['H_ratio']:.2f}× | C2={best_ch['C2_ratio']:.2f}× | Ni/Ne={best_ch['ni_ne_ratio']:.2f}")
    print(f"  f(x)={best_ch['objective']:.1f}")
    print()

    # Lowest CH
    lowest_ch = min(results, key=lambda x: x['CH'])
    print("LOWEST CH ACHIEVED:")
    print(f"  Ne={lowest_ch['ne']:.2e} cm⁻³")
    print(f"  CH={lowest_ch['CH']:.2e} cm⁻³ ({lowest_ch['CH_ratio']:.2f}×)")
    print(f"  H={lowest_ch['H_ratio']:.2f}× | C2={lowest_ch['C2_ratio']:.2f}× | Ni/Ne={lowest_ch['ni_ne_ratio']:.2f}")
    print(f"  f(x)={lowest_ch['objective']:.1f}")
    print()

    # Check if any results are in ALL target ranges
    perfect_results = [r for r in results if
                      0.6 <= r['H_ratio'] <= 1.4 and
                      0.6 <= r['CH_ratio'] <= 1.4 and
                      0.6 <= r['C2_ratio'] <= 1.4 and
                      2.0 <= r['ni_ne_ratio'] <= 6.0]

    if perfect_results:
        print("🎯 PERFECT RESULTS (ALL targets in range!):")
        print("-"*80)
        for r in sorted(perfect_results, key=lambda x: x['objective']):
            print(f"  Ne={r['ne']:.2e} | H={r['H_ratio']:.2f}× | CH={r['CH_ratio']:.2f}× | "
                  f"C2={r['C2_ratio']:.2f}× | Ni/Ne={r['ni_ne_ratio']:.2f} | f(x)={r['objective']:.1f}")
        print()
    else:
        print("No results with ALL targets in range.")
        print()

    print("="*80)
    print("FULL RESULTS TABLE")
    print("="*80)
    print(f"{'Ne':>10} | {'Ni/Ne':>7} | {'f(x)':>7} | {'H':>6} | {'CH':>6} | {'C2':>6}")
    print("-"*70)
    for r in results:
        ni_ne_ok = "✓" if 2.0 <= r['ni_ne_ratio'] <= 6.0 else " "
        print(f"{r['ne']:>10.2e} | {r['ni_ne_ratio']:>6.2f}{ni_ne_ok} | {r['objective']:>7.1f} | "
              f"{r['H_ratio']:>6.2f} | {r['CH_ratio']:>6.2f} | {r['C2_ratio']:>6.2f}")

print()
print("="*80)

# Save detailed results
with open('very_low_ne_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Detailed results saved to: very_low_ne_results.json")
