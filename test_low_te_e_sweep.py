#!/usr/bin/env python3
"""
Focused sweep at LOW Te (0.8-1.0 eV) with varying E_field
Strategy: At low Te, ionization is weak. Need to find E_field that:
  1. Gives positive Ni/Ne ratio (not negative!)
  2. Keeps Ni/Ne in reasonable range (ideally 2-6×)
  3. Minimizes CH production (low Te helps this)
  4. See if we can hit target densities
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

# Create test matrix: LOW Te and wide E_field range
te_values = [0.8, 0.9, 1.0]  # eV - LOW Te region
e_field_values = [50, 100, 150, 200, 250, 300, 400, 500]  # V/cm - wide range

print("="*80)
print("LOW Te & E_field SWEEP - Exploring low electron temperature regime")
print("="*80)
print()
print("Strategy:")
print("  1. Low Te → minimal e + CH4 → CH production (good for CH target!)")
print("  2. Low Te → weak ionization → need low E to avoid negative Ni/Ne")
print("  3. Find E_field that gives positive, reasonable Ni/Ne ratio")
print()
print(f"Testing Te values: {te_values} eV")
print(f"Testing E_field values: {e_field_values} V/cm")
print(f"Fixed Ne: 3.0e9 cm⁻³")
print(f"Target Ni/Ne: 2-6×")
print("="*80)
print()

results = []
test_count = 0
total_tests = len(te_values) * len(e_field_values)

for te in te_values:
    for e_field in e_field_values:
        test_count += 1
        print(f"[{test_count}/{total_tests}] Testing Te={te:.2f} eV, E={e_field} V/cm...", end='', flush=True)

        # Start from BEST parameters
        params = params_base.copy()
        params['rate_values'] = params_base['rate_values'].copy()

        # Apply BEST chemical parameters
        params['rate_values']['stick_C2_9_9'] = 1.25e3
        params['rate_values']['loss_C2_11_3'] = 100
        params['rate_values']['C2_H_CH_C_cm3_7_6'] = 8.0e-11
        params['rate_values']['stick_C2H2_9_11'] = 472

        # Apply Te, E_field, Ne
        params['Te'] = te
        params['E_field'] = e_field
        params['ne'] = 3.0e9  # User's educated guess

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
                'te': te,
                'e_field': e_field,
                'ne': ne_final,
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

            # Check if in target range
            ni_ne_ok = "✓" if 0 < ni_ne_ratio <= 10.0 else "✗"
            in_range = "✓" if 2.0 <= ni_ne_ratio <= 6.0 else " "

            print(f" {ni_ne_ok} Ni/Ne={ni_ne_ratio:+.2f}{in_range} | f(x)={objective:.0f} | H={H/TARGETS['H']:.2f}× | C2={C2/TARGETS['C2']:.2f}× | CH={CH/TARGETS['CH']:.1f}×")

        except Exception as e:
            print(f" ERROR: {e}")
            continue

print()
print("="*80)
print("BEST RESULTS BY CRITERIA")
print("="*80)
print()

if results:
    # Filter for positive Ni/Ne only
    positive_ni_results = [r for r in results if r['ni_ne_ratio'] > 0]

    if positive_ni_results:
        # Best overall objective with positive Ni/Ne
        best_obj = min(positive_ni_results, key=lambda x: x['objective'])
        print("BEST OBJECTIVE FUNCTION (positive Ni/Ne only):")
        print(f"  Te={best_obj['te']:.2f} eV, E={best_obj['e_field']} V/cm")
        print(f"  f(x)={best_obj['objective']:.1f} | Ni/Ne={best_obj['ni_ne_ratio']:.2f}")
        print(f"  H={best_obj['H_ratio']:.2f}×, C2={best_obj['C2_ratio']:.2f}×, CH={best_obj['CH_ratio']:.1f}×")
        print()

        # Best Ni/Ne ratio (closest to target 2-6, positive only)
        def ni_ne_score(r):
            ratio = r['ni_ne_ratio']
            if 2.0 <= ratio <= 6.0:
                return abs(ratio - 3.5)  # Prefer middle of range
            elif ratio < 2.0:
                return 2.0 - ratio + 10  # Penalize below range
            else:
                return ratio - 6.0 + 10  # Penalize above range

        best_ni_ne = min(positive_ni_results, key=ni_ne_score)
        print("BEST Ni/Ne RATIO (closest to target 2-6×):")
        print(f"  Te={best_ni_ne['te']:.2f} eV, E={best_ni_ne['e_field']} V/cm")
        print(f"  Ni/Ne={best_ni_ne['ni_ne_ratio']:.2f} | f(x)={best_ni_ne['objective']:.1f}")
        print(f"  H={best_ni_ne['H_ratio']:.2f}×, C2={best_ni_ne['C2_ratio']:.2f}×, CH={best_ni_ne['CH_ratio']:.1f}×")
        print()

        # Best CH (closest to 1.0×)
        best_ch = min(positive_ni_results, key=lambda x: abs(x['CH_ratio'] - 1.0))
        print("BEST CH (closest to target):")
        print(f"  Te={best_ch['te']:.2f} eV, E={best_ch['e_field']} V/cm")
        print(f"  CH={best_ch['CH_ratio']:.2f}× | Ni/Ne={best_ch['ni_ne_ratio']:.2f}")
        print(f"  H={best_ch['H_ratio']:.2f}×, C2={best_ch['C2_ratio']:.2f}×, f(x)={best_ch['objective']:.1f}")
        print()
    else:
        print("WARNING: No results with positive Ni/Ne ratio found!")
        print()

    print("="*80)
    print("SUMMARY TABLE (all results)")
    print("="*80)
    print(f"{'Te':>6} | {'E':>5} | {'Ni/Ne':>8} | {'f(x)':>8} | {'H':>6} | {'C2':>6} | {'CH':>7}")
    print("-"*80)

    for r in sorted(results, key=lambda x: (x['te'], x['e_field'])):
        ni_ok = "✓" if r['ni_ne_ratio'] > 0 else "✗"
        in_range = "✓" if 2.0 <= r['ni_ne_ratio'] <= 6.0 else " "
        print(f"{r['te']:>6.2f} | {r['e_field']:>5} | {r['ni_ne_ratio']:>7.2f}{ni_ok}{in_range} | {r['objective']:>8.0f} | "
              f"{r['H_ratio']:>6.2f} | {r['C2_ratio']:>6.2f} | {r['CH_ratio']:>7.1f}")

print()
print("="*80)

# Save detailed results
with open('low_te_e_sweep_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Detailed results saved to: low_te_e_sweep_results.json")
