#!/usr/bin/env python3
"""
Test with corrected CH + CH4 → C2H4 + H rate
Current: 1.50e-10 (12.5× above lit max!)
Should be: ≤ 1.20e-11 (literature maximum)

This is the dominant CH consumption pathway (75.7%!)
Correcting it will likely INCREASE CH significantly.
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
db = get_complete_rate_database()

TARGETS = {'H': 5.18e13, 'CH': 1.0e9, 'C2': 1.3e11}

print("="*80)
print("TEST: CORRECTED CH + CH4 RATE")
print("="*80)
print()
print("Issue: CH + CH4 → C2H4 + H is at 1.50e-10 cm³/s")
print("       Literature maximum: 1.20e-11 cm³/s (12.5× lower!)")
print("       This is 75.7% of ALL CH consumption!")
print()
print("Baulch et al. (2005) combustion chemistry database")
print(f"Gas temperature: {checkpoint['params']['Tgas']} K")
print()
print("Testing scenarios:")
print("  1. Current (UNCONSTRAINED): 1.50e-10")
print("  2. Literature MAXIMUM: 1.20e-11")
print("  3. Literature MINIMUM: 8.00e-12")
print("  4. Check if Tg-dependent formula exists")
print("="*80)
print()

def run_simulation(ch_ch4_rate, label):
    """Run simulation with specified CH+CH4 rate."""
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply ALL optimizations
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

    # KEY: Set the CH + CH4 rate
    params['rate_values']['CH_CH4_C2H4_H_cm3_7_20'] = ch_ch4_rate

    params['Te'] = 1.20
    params['E_field'] = 75
    params['ne'] = 1.0e8

    try:
        k_dict = define_rates_tunable(params)

        # Apply rate_values overrides
        for name, val in params.get('rate_values', {}).items():
            if name in k_dict:
                if name in db:
                    k_dict[name] = np.clip(val, db[name].min, db[name].max)
                else:
                    k_dict[name] = val

        # FORCE override the CH+CH4 rate (it's hard-coded in define_rates_tunable!)
        k_dict['CH_CH4_C2H4_H_cm3_7_20'] = ch_ch4_rate

        params['k'] = k_dict
        reactions, tags = build_reactions(params)
        params['R'] = reactions
        params['tags'] = tags

        ode = PlasmaODE(params)

        y0 = np.zeros(ode.ns)
        y0[ode.species.index('e')] = params['ne']
        y0[ode.species.index('Ar')] = 1.29e16
        y0[ode.species.index('CH4')] = 1.29e15

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

        y_final = sol.y[:, -1]
        H_idx = ode.species.index('H')
        CH_idx = ode.species.index('CH')
        C2_idx = ode.species.index('C2')

        H = y_final[H_idx]
        CH = y_final[CH_idx]
        C2 = y_final[C2_idx]

        ni_total = 0
        for i_sp, species_name in enumerate(ode.species):
            if species_name.endswith('Plus'):
                ni_total += y_final[i_sp]
            elif species_name.endswith('Minus'):
                ni_total -= y_final[i_sp]

        ni_ne_ratio = ni_total / y_final[ode.species.index('e')]

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
            'label': label,
            'rate': ch_ch4_rate,
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
        print(f"Error: {e}")
        return None

# Test scenarios
scenarios = [
    (1.50e-10, "Current (UNCONSTRAINED)"),
    (1.20e-11, "Literature MAXIMUM"),
    (1.00e-11, "Mid-range"),
    (8.00e-12, "Literature MINIMUM"),
]

results = []

for rate, label in scenarios:
    print(f"Testing: {label} (k = {rate:.2e} cm³/s)")
    result = run_simulation(rate, label)

    if result:
        results.append(result)
        ni_ne_ok = "✓" if 2.0 <= result['ni_ne_ratio'] <= 6.0 else " "
        print(f"  CH = {result['CH']:.2e} cm⁻³ ({result['CH_ratio']:.2f}×)")
        print(f"  H = {result['H']:.2e} cm⁻³ ({result['H_ratio']:.2f}×)")
        print(f"  C2 = {result['C2']:.2e} cm⁻³ ({result['C2_ratio']:.2f}×)")
        print(f"  Ni/Ne = {result['ni_ne_ratio']:.2f}{ni_ne_ok}")
        print(f"  f(x) = {result['objective']:.1f}")
    else:
        print(f"  FAILED!")
    print()

print("="*80)
print("COMPARISON")
print("="*80)
print()

if len(results) >= 2:
    baseline = results[0]  # Current unconstrained
    corrected = results[1]  # Literature max

    print(f"{'Scenario':<30} | {'CH ratio':>9} | {'f(x)':>7} | {'Status'}")
    print("-"*80)
    for r in results:
        ch_ok = "✓" if 0.6 <= r['CH_ratio'] <= 1.4 else "✗"
        print(f"{r['label']:<30} | {r['CH_ratio']:9.2f} | {r['objective']:7.1f} | {ch_ok}")

    print()
    print("IMPACT OF CORRECTION:")
    print("-"*80)
    ch_change = (corrected['CH'] / baseline['CH'] - 1) * 100
    print(f"CH increase: {ch_change:+.1f}%")
    print(f"  From: {baseline['CH']:.2e} ({baseline['CH_ratio']:.2f}×)")
    print(f"  To:   {corrected['CH']:.2e} ({corrected['CH_ratio']:.2f}×)")
    print()

    if ch_change > 100:
        print("⚠️  WARNING: Enforcing literature constraint causes CH to MORE THAN DOUBLE!")
        print("    This means the unconstrained rate was artificially suppressing CH.")
        print("    The 'best' result of CH=4.18× is NOT physically realistic!")
    elif ch_change > 50:
        print("⚠️  CAUTION: CH increases significantly when corrected.")
        print("    The unconstrained rate was masking true CH density.")
    elif ch_change > 20:
        print("⚠️  NOTICE: Moderate CH increase when enforcing lit constraint.")
    else:
        print("✓  Minor impact from correction.")

    print()
    print("CONCLUSION:")
    print("-"*80)
    if corrected['CH_ratio'] > 1.4:
        print(f"With properly constrained CH+CH4 rate:")
        print(f"  CH = {corrected['CH_ratio']:.2f}× (ABOVE 1.4× target)")
        print(f"  The true CH density is HIGHER than reported!")
        print()
        print(f"The previous 'best' result (CH=4.18×) was achieved by using")
        print(f"an unphysically high CH consumption rate (12.5× lit max).")
        print()
        print(f"TRUE best result with literature-constrained rates:")
        print(f"  CH = {corrected['CH']:.2e} cm⁻³ ({corrected['CH_ratio']:.2f}× target)")

print()
print("="*80)

# Save results
with open('corrected_ch_ch4_rate_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Results saved to: corrected_ch_ch4_rate_results.json")
