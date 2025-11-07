#!/usr/bin/env python3
"""
Refine Te and E sweep from current best sweet spot
Current best: Te=1.20 eV, E=75 V/cm, Ne=1.0e8
Goal: Find conditions that further reduce CH while maintaining other targets
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

# FIXED condition
NE_BEST = 1.0e8

# Parameter sweeps - refined around sweet spot
te_values = [
    1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6
]

e_field_values = [
    50, 60, 70, 75, 80, 90, 100, 125, 150, 175, 200
]

print("="*80)
print("TE AND E REFINEMENT SWEEP")
print("="*80)
print()
print(f"Fixed: Ne={NE_BEST:.1e} cm⁻³")
print(f"Te sweep: {len(te_values)} values from {min(te_values)} to {max(te_values)} eV")
print(f"E sweep: {len(e_field_values)} values from {min(e_field_values)} to {max(e_field_values)} V/cm")
print(f"Total combinations: {len(te_values) * len(e_field_values)}")
print()
print("All optimizations applied (including CH consumption maximized)")
print("="*80)
print()

def run_simulation(te, e_field):
    """Run simulation and return results."""
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply ALL OPTIMIZATIONS (including CH consumption maximized)
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

    # CH consumption maximized
    params['rate_values']['loss_CH_11_9'] = 1.00e4
    params['rate_values']['stick_CH_9_3'] = 6.25e3
    params['rate_values']['CH_CH4_CH2_CH3_cm3_7_39'] = 1.20e-11
    params['rate_values']['CH_H_CH2_cm3_7_21'] = 1.20e-10
    params['rate_values']['CH_CH3_C2H3_H_cm3_7_10'] = 1.20e-10
    params['rate_values']['CH_CH3_C2H2_H2_cm3_7_23'] = 1.20e-10
    params['rate_values']['CH_H2_CH2_H_cm3_7_30'] = 1.20e-11
    params['rate_values']['CH_H2_CH2_H_cm3_7_51'] = 1.20e-11
    params['rate_values']['C2H2_CH_C3_H2_H_cm3_7_56'] = 1.20e-10
    params['rate_values']['CH_C2H2_C3H2_H_cm3_7_22'] = 1.20e-10

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
            'Te': te,
            'E_field': e_field,
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

# Run sweep
results = []
total = len(te_values) * len(e_field_values)
count = 0

for te in te_values:
    for e_field in e_field_values:
        count += 1
        print(f"[{count}/{total}] Te={te:.2f} eV, E={e_field} V/cm ... ", end='', flush=True)

        result = run_simulation(te, e_field)

        if result is None:
            print("FAILED")
            continue

        results.append(result)

        ni_ne_ok = "✓" if 2.0 <= result['ni_ne_ratio'] <= 6.0 else " "
        print(f"Ni/Ne={result['ni_ne_ratio']:.2f}{ni_ne_ok} | f(x)={result['objective']:.1f} | "
              f"H={result['H_ratio']:.2f}× | CH={result['CH_ratio']:.2f}× | C2={result['C2_ratio']:.2f}×")

print()
print("="*80)
print("ANALYSIS OF RESULTS")
print("="*80)
print()

if results:
    # Best CH (lowest)
    best_ch = min(results, key=lambda x: x['CH'])
    print("BEST (LOWEST) CH:")
    print(f"  Te={best_ch['Te']:.2f} eV, E={best_ch['E_field']} V/cm")
    print(f"  CH = {best_ch['CH']:.2e} cm⁻³ ({best_ch['CH_ratio']:.2f}×)")
    print(f"  H = {best_ch['H']:.2e} cm⁻³ ({best_ch['H_ratio']:.2f}×)")
    print(f"  C2 = {best_ch['C2']:.2e} cm⁻³ ({best_ch['C2_ratio']:.2f}×)")
    print(f"  Ni/Ne = {best_ch['ni_ne_ratio']:.2f}")
    print(f"  f(x) = {best_ch['objective']:.1f}")
    print()

    # Best objective
    best_obj = min(results, key=lambda x: x['objective'])
    print("BEST OBJECTIVE:")
    print(f"  Te={best_obj['Te']:.2f} eV, E={best_obj['E_field']} V/cm")
    print(f"  f(x) = {best_obj['objective']:.1f}")
    print(f"  H={best_obj['H_ratio']:.2f}× | CH={best_obj['CH_ratio']:.2f}× | C2={best_obj['C2_ratio']:.2f}×")
    print(f"  Ni/Ne = {best_obj['ni_ne_ratio']:.2f}")
    print()

    # Check if any achieve ALL targets
    perfect_results = [r for r in results if
                      0.6 <= r['H_ratio'] <= 1.4 and
                      0.6 <= r['CH_ratio'] <= 1.4 and
                      0.6 <= r['C2_ratio'] <= 1.4 and
                      2.0 <= r['ni_ne_ratio'] <= 6.0]

    if perfect_results:
        print("🎯🎯🎯 PERFECT RESULTS (ALL 4 TARGETS IN RANGE!):")
        print("-"*80)
        for r in sorted(perfect_results, key=lambda x: x['objective']):
            print(f"  Te={r['Te']:.2f} eV, E={r['E_field']} V/cm")
            print(f"    H={r['H_ratio']:.2f}× | CH={r['CH_ratio']:.2f}× | C2={r['C2_ratio']:.2f}× | Ni/Ne={r['ni_ne_ratio']:.2f}")
            print(f"    f(x)={r['objective']:.1f}")
            print()
    else:
        print("No results with ALL 4 targets in range.")
        print()

        # Show results with 3 targets in range
        three_targets = [r for r in results if
                        sum([0.6 <= r['H_ratio'] <= 1.4,
                             0.6 <= r['CH_ratio'] <= 1.4,
                             0.6 <= r['C2_ratio'] <= 1.4]) >= 2 and
                        2.0 <= r['ni_ne_ratio'] <= 6.0]

        if three_targets:
            print(f"Results with ≥2 species targets + Ni/Ne in range: {len(three_targets)}")
            best_3 = sorted(three_targets, key=lambda x: x['objective'])[:5]
            print("Top 5:")
            for r in best_3:
                h_ok = "✓" if 0.6 <= r['H_ratio'] <= 1.4 else " "
                ch_ok = "✓" if 0.6 <= r['CH_ratio'] <= 1.4 else " "
                c2_ok = "✓" if 0.6 <= r['C2_ratio'] <= 1.4 else " "
                print(f"  Te={r['Te']:.2f} eV, E={r['E_field']:3.0f} V/cm | "
                      f"f(x)={r['objective']:.1f} | "
                      f"H={r['H_ratio']:.2f}×{h_ok} CH={r['CH_ratio']:.2f}×{ch_ok} C2={r['C2_ratio']:.2f}×{c2_ok}")

    # CH statistics
    ch_ratios = [r['CH_ratio'] for r in results]
    print()
    print(f"CH STATISTICS:")
    print(f"  Min:  {min(ch_ratios):.2f}× (Te={best_ch['Te']:.2f} eV, E={best_ch['E_field']} V/cm)")
    print(f"  Max:  {max(ch_ratios):.2f}×")
    print(f"  Mean: {np.mean(ch_ratios):.2f}×")
    print(f"  Std:  {np.std(ch_ratios):.2f}×")
    print(f"  Range: {max(ch_ratios) - min(ch_ratios):.2f}×")
    print()

print("="*80)

# Save results
with open('te_e_refinement_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to: te_e_refinement_results.json")
