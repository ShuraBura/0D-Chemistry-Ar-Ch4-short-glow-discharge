#!/usr/bin/env python3
"""
Test VERY LOW Ne values (below 1e8)
Since lower Ne gives better CH control, explore if trend continues
Goal: Find if we can reduce CH below 4.18× by going even lower
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

# FIXED conditions
TE_BEST = 1.20
E_BEST = 75

# Very low Ne sweep - going below 1e8
ne_values = [
    # Very low range
    5e6, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7, 7e7, 8e7, 9e7,
    # Baseline for comparison
    1.0e8, 1.5e8, 2.0e8, 2.5e8, 3.0e8
]

print("="*80)
print("VERY LOW NE SWEEP (BELOW 1e8)")
print("="*80)
print()
print(f"Fixed: Te={TE_BEST} eV, E={E_BEST} V/cm")
print(f"Ne sweep: {len(ne_values)} values from {min(ne_values):.1e} to {max(ne_values):.1e} cm⁻³")
print()
print("Hypothesis: Lower Ne → Less e+CH4→CH → Lower CH density")
print("Testing if CH can be reduced below 4.18× by going to very low Ne")
print("All optimizations applied (including CH consumption maximized)")
print("="*80)
print()

def run_simulation(ne):
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

# Run sweep
results = []
total = len(ne_values)

for i, ne in enumerate(ne_values):
    print(f"[{i+1}/{total}] Ne={ne:.1e} cm⁻³ ... ", end='', flush=True)

    result = run_simulation(ne)

    if result is None:
        print("FAILED")
        continue

    results.append(result)

    ni_ne_ok = "✓" if 2.0 <= result['ni_ne_ratio'] <= 6.0 else " "

    # Highlight if CH is lower than 4.18×
    ch_marker = " 🎯" if result['CH_ratio'] < 4.18 else ""

    print(f"Ni/Ne={result['ni_ne_ratio']:.2f}{ni_ne_ok} | f(x)={result['objective']:.1f} | "
          f"H={result['H_ratio']:.2f}× | CH={result['CH_ratio']:.2f}×{ch_marker} | C2={result['C2_ratio']:.2f}×")

print()
print("="*80)
print("ANALYSIS OF RESULTS")
print("="*80)
print()

if results:
    # Best CH (lowest)
    best_ch = min(results, key=lambda x: x['CH'])
    print("BEST (LOWEST) CH:")
    print(f"  Ne={best_ch['ne']:.1e} cm⁻³")
    print(f"  CH = {best_ch['CH']:.2e} cm⁻³ ({best_ch['CH_ratio']:.2f}×)")
    print(f"  H = {best_ch['H']:.2e} cm⁻³ ({best_ch['H_ratio']:.2f}×)")
    print(f"  C2 = {best_ch['C2']:.2e} cm⁻³ ({best_ch['C2_ratio']:.2f}×)")
    print(f"  Ni/Ne = {best_ch['ni_ne_ratio']:.2f}")
    print(f"  f(x) = {best_ch['objective']:.1f}")

    # Compare to baseline 1e8
    baseline = [r for r in results if r['ne'] == 1.0e8][0] if any(r['ne'] == 1.0e8 for r in results) else None
    if baseline and best_ch['ne'] != 1.0e8:
        ch_improvement = (1 - best_ch['CH'] / baseline['CH']) * 100
        print(f"  Improvement vs Ne=1e8: {ch_improvement:+.1f}% CH reduction")
    print()

    # Best objective
    best_obj = min(results, key=lambda x: x['objective'])
    print("BEST OBJECTIVE:")
    print(f"  Ne={best_obj['ne']:.1e} cm⁻³")
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
            print(f"  Ne={r['ne']:.1e} cm⁻³")
            print(f"    H={r['H_ratio']:.2f}× | CH={r['CH_ratio']:.2f}× | C2={r['C2_ratio']:.2f}× | Ni/Ne={r['ni_ne_ratio']:.2f}")
            print(f"    f(x)={r['objective']:.1f}")
            print()
    else:
        print("No results with ALL 4 targets in range.")
        print()

    # CH vs Ne trend
    print("DETAILED CH TREND WITH NE:")
    print(f"{'Ne (cm⁻³)':>12} | {'CH ratio':>9} | {'Ni/Ne':>6} | {'f(x)':>7} | Notes")
    print("-"*70)
    for r in results:
        ni_ne_ok = "✓" if 2.0 <= r['ni_ne_ratio'] <= 6.0 else "✗"
        ch_best = "← BEST CH" if r['CH'] == best_ch['CH'] else ""
        print(f"{r['ne']:12.1e} | {r['CH_ratio']:9.2f} | {r['ni_ne_ratio']:5.2f}{ni_ne_ok} | {r['objective']:7.1f} | {ch_best}")

    # Statistics
    ch_ratios = [r['CH_ratio'] for r in results]
    print()
    print("CH STATISTICS:")
    print(f"  Min:  {min(ch_ratios):.2f}× at Ne={best_ch['ne']:.1e} cm⁻³")
    print(f"  Max:  {max(ch_ratios):.2f}×")
    print(f"  Range: {max(ch_ratios) - min(ch_ratios):.2f}×")
    print()

    # Check if trend continues
    print("TREND ANALYSIS:")
    if best_ch['ne'] < 1.0e8:
        print(f"  ✓ Lower Ne continues to reduce CH!")
        print(f"  ✓ New best at Ne={best_ch['ne']:.1e}: CH={best_ch['CH_ratio']:.2f}×")
    else:
        print(f"  ✗ Trend stops at Ne=1e8")
        print(f"  Going lower does not improve CH")
    print()

    # Ni/Ne check
    ni_ne_ratios = [r['ni_ne_ratio'] for r in results]
    print("NI/NE STATISTICS:")
    print(f"  Min:  {min(ni_ne_ratios):.2f}")
    print(f"  Max:  {max(ni_ne_ratios):.2f}")
    print(f"  Mean: {np.mean(ni_ne_ratios):.2f}")
    out_of_range = sum(1 for r in ni_ne_ratios if r < 2.0 or r > 6.0)
    if out_of_range > 0:
        print(f"  WARNING: {out_of_range}/{len(results)} results have Ni/Ne out of range!")
    print()

print("="*80)

# Save results
with open('very_low_ne_extended_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to: very_low_ne_extended_results.json")
