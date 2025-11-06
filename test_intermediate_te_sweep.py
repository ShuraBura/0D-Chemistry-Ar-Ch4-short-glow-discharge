#!/usr/bin/env python3
"""
INTERMEDIATE Te sweep (1.1-1.4 eV) - Looking for the sweet spot!
Strategy: Find Te where:
  1. Ionization strong enough for Ni/Ne = 2-6×
  2. But CH production still controlled (not exploding like at 1.8+ eV)
  3. Best chance to balance all requirements
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

# Create test matrix: INTERMEDIATE Te and E_field range
te_values = [1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4]  # eV - intermediate range
e_field_values = [100, 150, 200, 250, 300, 350, 400, 450, 500]  # V/cm

print("="*80)
print("INTERMEDIATE Te & E_field SWEEP - Searching for the sweet spot!")
print("="*80)
print()
print("Goal: Find conditions where:")
print("  1. Ni/Ne = 2-6× (proper sheath-edge physics)")
print("  2. CH stays controlled (< 10× target)")
print("  3. H and C2 as close to target as possible")
print()
print("Known regimes:")
print("  Te < 1.0 eV: Ni/Ne too LOW (0.2-0.4)")
print("  Te = 1.47 eV: Ni/Ne good (4-6×) but CH=16×")
print("  Te > 1.8 eV: CH explodes (64-143×), Ni/Ne too HIGH")
print()
print(f"Testing Te values: {[f'{t:.2f}' for t in te_values]} eV")
print(f"Testing E_field values: {e_field_values} V/cm")
print(f"Fixed Ne: 3.0e9 cm⁻³")
print("="*80)
print()

results = []
test_count = 0
total_tests = len(te_values) * len(e_field_values)

for te in te_values:
    for e_field in e_field_values:
        test_count += 1
        print(f"[{test_count}/{total_tests}] Te={te:.2f} eV, E={e_field:3d} V/cm...", end='', flush=True)

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
        params['ne'] = 3.0e9

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

            # Check if in target ranges
            ni_ne_ok = "✓" if 2.0 <= ni_ne_ratio <= 6.0 else " "
            ch_ok = "✓" if ni_ne_ratio <= 10.0 else " "

            print(f" Ni/Ne={ni_ne_ratio:5.2f}{ni_ne_ok} | CH={CH/TARGETS['CH']:5.1f}×{ch_ok} | "
                  f"f(x)={objective:5.0f} | H={H/TARGETS['H']:.2f}× | C2={C2/TARGETS['C2']:.2f}×")

        except Exception as e:
            print(f" ERROR: {e}")
            continue

print()
print("="*80)
print("BEST RESULTS BY CRITERIA")
print("="*80)
print()

if results:
    # Filter for results with BOTH good Ni/Ne AND reasonable CH
    sweet_spot_results = [r for r in results
                          if 2.0 <= r['ni_ne_ratio'] <= 6.0 and r['CH_ratio'] < 10.0]

    if sweet_spot_results:
        print("🎯 SWEET SPOT FOUND! (Ni/Ne in 2-6× AND CH < 10×):")
        print("-" * 80)
        best_sweet = min(sweet_spot_results, key=lambda x: x['objective'])
        print(f"  Te={best_sweet['te']:.2f} eV, E={best_sweet['e_field']} V/cm")
        print(f"  Ni/Ne={best_sweet['ni_ne_ratio']:.2f} ✓ | CH={best_sweet['CH_ratio']:.1f}× ✓")
        print(f"  f(x)={best_sweet['objective']:.1f}")
        print(f"  H={best_sweet['H_ratio']:.2f}×, C2={best_sweet['C2_ratio']:.2f}×")
        print()
    else:
        print("No sweet spot found (Ni/Ne=2-6× AND CH<10× simultaneously)")
        print()

    # Best overall objective
    best_obj = min(results, key=lambda x: x['objective'])
    print("BEST OBJECTIVE FUNCTION:")
    print(f"  Te={best_obj['te']:.2f} eV, E={best_obj['e_field']} V/cm")
    print(f"  f(x)={best_obj['objective']:.1f}")
    print(f"  Ni/Ne={best_obj['ni_ne_ratio']:.2f} | CH={best_obj['CH_ratio']:.1f}×")
    print(f"  H={best_obj['H_ratio']:.2f}×, C2={best_obj['C2_ratio']:.2f}×")
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
    print(f"  Te={best_ni_ne['te']:.2f} eV, E={best_ni_ne['e_field']} V/cm")
    print(f"  Ni/Ne={best_ni_ne['ni_ne_ratio']:.2f} | f(x)={best_ni_ne['objective']:.1f}")
    print(f"  H={best_ni_ne['H_ratio']:.2f}×, C2={best_ni_ne['C2_ratio']:.2f}×, CH={best_ni_ne['CH_ratio']:.1f}×")
    print()

    # Best CH control
    best_ch = min(results, key=lambda x: abs(x['CH_ratio'] - 1.0))
    print("BEST CH CONTROL (closest to 1.0×):")
    print(f"  Te={best_ch['te']:.2f} eV, E={best_ch['e_field']} V/cm")
    print(f"  CH={best_ch['CH_ratio']:.2f}× | Ni/Ne={best_ch['ni_ne_ratio']:.2f}")
    print(f"  H={best_ch['H_ratio']:.2f}×, C2={best_ch['C2_ratio']:.2f}×, f(x)={best_ch['objective']:.1f}")
    print()

print()
print("="*80)
print("SUMMARY: Conditions with Ni/Ne in target range (2-6×)")
print("="*80)

target_ni_results = [r for r in results if 2.0 <= r['ni_ne_ratio'] <= 6.0]
if target_ni_results:
    print(f"{'Te':>6} | {'E':>5} | {'Ni/Ne':>7} | {'CH':>7} | {'f(x)':>8} | {'H':>6} | {'C2':>6}")
    print("-"*80)
    for r in sorted(target_ni_results, key=lambda x: (x['te'], x['e_field'])):
        ch_mark = "✓" if r['CH_ratio'] < 10.0 else "✗"
        print(f"{r['te']:>6.2f} | {r['e_field']:>5} | {r['ni_ne_ratio']:>7.2f} | "
              f"{r['CH_ratio']:>6.1f}×{ch_mark} | {r['objective']:>8.0f} | "
              f"{r['H_ratio']:>6.2f} | {r['C2_ratio']:>6.2f}")
else:
    print("No results with Ni/Ne in target range (2-6×)")

print()
print("="*80)

# Save detailed results
with open('intermediate_te_sweep_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Detailed results saved to: intermediate_te_sweep_results.json")
