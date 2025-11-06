#!/usr/bin/env python3
"""
Te & E_field sweep with ALL current optimizations
Goal: Find conditions that reduce CH while maintaining H and C2
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

# Fixed Ne (found to be optimal)
NE_BEST = 5.0e8

# Test ranges
te_values = [1.0, 1.1, 1.2, 1.25, 1.3, 1.35, 1.4, 1.5]  # eV
e_field_values = [50, 75, 100, 125, 150, 175, 200, 250]  # V/cm

print("="*80)
print("Te & E_field SWEEP WITH OPTIMIZED CHEMISTRY")
print("="*80)
print()
print(f"Fixed Ne: {NE_BEST:.1e} cm⁻³")
print(f"Testing Te: {te_values} eV")
print(f"Testing E:  {e_field_values} V/cm")
print(f"Total combinations: {len(te_values) * len(e_field_values)}")
print()
print("Applied optimizations:")
print("  - stick_H = 3.89e2 (reduced)")
print("  - stick_C2 = 1.25e3 (reduced)")
print("  - loss_C2 = 100 (reduced)")
print("  - C2H2 losses minimized")
print("  - CH production minimized (e+CH4→CH)")
print("="*80)
print()

def run_simulation(te, e_field):
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
    params['Te'] = te
    params['E_field'] = e_field
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
            'te': te,
            'e_field': e_field,
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
test_count = 0
total_tests = len(te_values) * len(e_field_values)

for te in te_values:
    for e_field in e_field_values:
        test_count += 1
        print(f"[{test_count}/{total_tests}] Te={te:.2f} eV, E={e_field:3d} V/cm...", end='', flush=True)

        result = run_simulation(te, e_field)

        if result is None:
            print(" FAILED")
            continue

        results.append(result)

        ni_ne_ok = "✓" if 2.0 <= result['ni_ne_ratio'] <= 6.0 else " "

        print(f" Ni/Ne={result['ni_ne_ratio']:5.2f}{ni_ne_ok} | f(x)={result['objective']:6.1f} | "
              f"H={result['H_ratio']:.2f}× | CH={result['CH_ratio']:.1f}× | C2={result['C2_ratio']:.2f}×")

print()
print("="*80)
print("ANALYSIS OF RESULTS")
print("="*80)
print()

if results:
    # Best objective
    best_obj = min(results, key=lambda x: x['objective'])
    print("BEST OBJECTIVE:")
    print(f"  Te={best_obj['te']:.2f} eV, E={best_obj['e_field']} V/cm")
    print(f"  f(x)={best_obj['objective']:.1f} | Ni/Ne={best_obj['ni_ne_ratio']:.2f}")
    print(f"  H={best_obj['H_ratio']:.2f}× | CH={best_obj['CH_ratio']:.1f}× | C2={best_obj['C2_ratio']:.2f}×")
    print()

    # Best CH (closest to 1.0)
    best_ch = min(results, key=lambda x: abs(x['CH_ratio'] - 1.0))
    print("BEST CH (closest to target):")
    print(f"  Te={best_ch['te']:.2f} eV, E={best_ch['e_field']} V/cm")
    print(f"  CH={best_ch['CH_ratio']:.2f}× | Ni/Ne={best_ch['ni_ne_ratio']:.2f}")
    print(f"  H={best_ch['H_ratio']:.2f}× | C2={best_ch['C2_ratio']:.2f}× | f(x)={best_ch['objective']:.1f}")
    print()

    # Best C2 (highest)
    best_c2 = max(results, key=lambda x: x['C2'])
    print("BEST C2 (highest):")
    print(f"  Te={best_c2['te']:.2f} eV, E={best_c2['e_field']} V/cm")
    print(f"  C2={best_c2['C2_ratio']:.2f}× | Ni/Ne={best_c2['ni_ne_ratio']:.2f}")
    print(f"  H={best_c2['H_ratio']:.2f}× | CH={best_c2['CH_ratio']:.1f}× | f(x)={best_c2['objective']:.1f}")
    print()

    # Best Ni/Ne (closest to 3.5)
    def ni_ne_score(r):
        ratio = r['ni_ne_ratio']
        if 2.0 <= ratio <= 6.0:
            return abs(ratio - 3.5)
        elif ratio < 2.0:
            return 2.0 - ratio + 10
        else:
            return ratio - 6.0 + 10

    best_ni_ne = min(results, key=ni_ne_score)
    print("BEST Ni/Ne (closest to 2-6× range):")
    print(f"  Te={best_ni_ne['te']:.2f} eV, E={best_ni_ne['e_field']} V/cm")
    print(f"  Ni/Ne={best_ni_ne['ni_ne_ratio']:.2f} | f(x)={best_ni_ne['objective']:.1f}")
    print(f"  H={best_ni_ne['H_ratio']:.2f}× | CH={best_ni_ne['CH_ratio']:.1f}× | C2={best_ni_ne['C2_ratio']:.2f}×")
    print()

    # Results with Ni/Ne in range
    in_range = [r for r in results if 2.0 <= r['ni_ne_ratio'] <= 6.0]
    if in_range:
        print(f"Results with Ni/Ne in target range (2-6×): {len(in_range)}/{len(results)}")
        print()
        print("TOP 10 by objective (with good Ni/Ne):")
        print(f"{'Te':>6} | {'E':>5} | {'Ni/Ne':>7} | {'f(x)':>7} | {'H':>6} | {'CH':>7} | {'C2':>6}")
        print("-"*70)
        for r in sorted(in_range, key=lambda x: x['objective'])[:10]:
            print(f"{r['te']:>6.2f} | {r['e_field']:>5} | {r['ni_ne_ratio']:>7.2f} | "
                  f"{r['objective']:>7.1f} | {r['H_ratio']:>6.2f} | {r['CH_ratio']:>6.1f}× | {r['C2_ratio']:>6.2f}")

print()
print("="*80)

# Save detailed results
with open('te_e_sweep_optimized_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Detailed results saved to: te_e_sweep_optimized_results.json")
