#!/usr/bin/env python3
"""
Comprehensive C2H2 optimization at BEST conditions (Te=1.25, E=100, Ne=5e8)
Strategy: Boost C2H2 → Boost C2 (since 92% of C2 comes from C2H2 + H → C2 + H2)

Approach:
1. Reduce C2H2 wall losses (stick_C2H2 to literature minimum)
2. Boost C2H2 PRODUCTION reactions (CH + CH2 → C2H2, etc.)
3. Reduce C2H2 consumption that doesn't produce C2 (CH + C2H2 → C3 species)
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

# BEST conditions found so far
TE_BEST = 1.25    # eV
E_BEST = 100      # V/cm
NE_BEST = 5.0e8   # cm⁻³

print("="*80)
print("COMPREHENSIVE C2H2 OPTIMIZATION")
print("="*80)
print()
print("STRATEGY: Boost C2H2 → Boost C2 (C2H2 + H → C2 + H2 is 92% of C2 production)")
print()
print(f"FIXED: Te={TE_BEST} eV, E={E_BEST} V/cm, Ne={NE_BEST:.1e} cm⁻³")
print()
print("Testing combinations:")
print("  1. Baseline (current best)")
print("  2. Further reduce stick_C2H2 to literature minimum (500)")
print("  3. Boost C2H2 production reactions (CH+CH2→C2H2, CH2+CH2→C2H2, etc.)")
print("  4. Reduce C2H2 consumption to C3 species (CH+C2H2→C3H2, C3H, etc.)")
print("  5. Combine all optimizations")
print("="*80)
print()

# Test configurations
configs = [
    {
        'name': 'Baseline (current best)',
        'description': 'stick_C2=1.25e3, loss_C2=100, C2_H_CH_C=8e-11, stick_C2H2=472',
        'rate_values': {
            'stick_C2_9_9': 1.25e3,
            'loss_C2_11_3': 100,
            'C2_H_CH_C_cm3_7_6': 8.0e-11,
            'stick_C2H2_9_11': 472,
        }
    },
    {
        'name': 'Reduce stick_C2H2 to minimum',
        'description': 'stick_C2H2 = 500 (literature minimum)',
        'rate_values': {
            'stick_C2_9_9': 1.25e3,
            'loss_C2_11_3': 100,
            'C2_H_CH_C_cm3_7_6': 8.0e-11,
            'stick_C2H2_9_11': 500,  # Literature minimum
        }
    },
    {
        'name': 'Boost C2H2 production',
        'description': 'CH+CH2→C2H2, CH2+CH2→C2H2, CH3+CH→C2H2 to max',
        'rate_values': {
            'stick_C2_9_9': 1.25e3,
            'loss_C2_11_3': 100,
            'C2_H_CH_C_cm3_7_6': 8.0e-11,
            'stick_C2H2_9_11': 472,
            # Boost C2H2 production to max
            'CH_CH2_C2H2_H_cm3_7_7': 1.20e-10,      # max (baseline)
            'CH2_CH2_C2H2_H2_cm3_7_15': 1.20e-10,   # max (boost from 1e-10)
            'CH3_CH_C2H2_H2_cm3_7_16': 1.20e-10,    # max (baseline)
            'CH2_C_C2H2_cm3_7_17': 1.20e-10,        # max (boost from 1e-10)
        }
    },
    {
        'name': 'Reduce C2H2 → C3 consumption',
        'description': 'CH+C2H2→C3 species to minimum',
        'rate_values': {
            'stick_C2_9_9': 1.25e3,
            'loss_C2_11_3': 100,
            'C2_H_CH_C_cm3_7_6': 8.0e-11,
            'stick_C2H2_9_11': 472,
            # Reduce C2H2 consumption that doesn't produce C2
            'CH_C2H2_C3H2_H_cm3_7_22': 8.0e-11,     # min (from 1e-10)
            'CH_C2H2_C3H_H2_cm3_7_27': 8.0e-11,     # min (from 1e-10)
            'CH_C2H2_C2H_CH2_cm3_7_29': 8.0e-11,    # min (from 1e-10)
        }
    },
    {
        'name': 'COMBINED: All C2H2 optimizations',
        'description': 'stick_C2H2 min + production max + consumption min',
        'rate_values': {
            'stick_C2_9_9': 1.25e3,
            'loss_C2_11_3': 100,
            'C2_H_CH_C_cm3_7_6': 8.0e-11,
            # C2H2 wall loss
            'stick_C2H2_9_11': 500,  # Literature minimum
            # Boost C2H2 production
            'CH_CH2_C2H2_H_cm3_7_7': 1.20e-10,
            'CH2_CH2_C2H2_H2_cm3_7_15': 1.20e-10,
            'CH3_CH_C2H2_H2_cm3_7_16': 1.20e-10,
            'CH2_C_C2H2_cm3_7_17': 1.20e-10,
            # Reduce C2H2 → C3 consumption
            'CH_C2H2_C3H2_H_cm3_7_22': 8.0e-11,
            'CH_C2H2_C3H_H2_cm3_7_27': 8.0e-11,
            'CH_C2H2_C2H_CH2_cm3_7_29': 8.0e-11,
        }
    },
]

results = []

for i, config in enumerate(configs):
    print(f"[{i+1}/{len(configs)}] Testing: {config['name']}")
    print(f"  {config['description']}")
    print("  ", end='', flush=True)

    # Start from base parameters
    params = params_base.copy()
    params['rate_values'] = config['rate_values'].copy()

    # Apply BEST conditions
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
            print(f"FAILED: {sol.message}")
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
        for i_sp, species_name in enumerate(ode.species):
            if species_name.endswith('Plus'):
                ni_total += y_final[i_sp]
            elif species_name.endswith('Minus'):
                ni_total -= y_final[i_sp]

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
            'name': config['name'],
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

        ni_ne_ok = "✓" if 2.0 <= ni_ne_ratio <= 6.0 else " "

        print(f"Ni/Ne={ni_ne_ratio:.2f}{ni_ne_ok} | f(x)={objective:.1f} | "
              f"H={H/TARGETS['H']:.2f}× | CH={CH/TARGETS['CH']:.1f}× | "
              f"C2={C2/TARGETS['C2']:.2f}× | C2H2={C2H2:.2e}")

    except Exception as e:
        print(f"ERROR: {e}")
        continue

    print()

print("="*80)
print("COMPARISON OF RESULTS")
print("="*80)
print()

if results:
    # Best objective
    best_obj = min(results, key=lambda x: x['objective'])
    print("BEST OBJECTIVE:")
    print(f"  {best_obj['name']}")
    print(f"  f(x)={best_obj['objective']:.1f} | Ni/Ne={best_obj['ni_ne_ratio']:.2f}")
    print(f"  H={best_obj['H_ratio']:.2f}× | CH={best_obj['CH_ratio']:.1f}× | "
          f"C2={best_obj['C2_ratio']:.2f}×")
    print(f"  C2H2={best_obj['C2H2']:.2e} cm⁻³")
    print()

    # Best C2
    best_c2 = max(results, key=lambda x: x['C2'])
    print("BEST C2:")
    print(f"  {best_c2['name']}")
    print(f"  C2={best_c2['C2_ratio']:.2f}× | C2H2={best_c2['C2H2']:.2e}")
    print(f"  f(x)={best_c2['objective']:.1f} | Ni/Ne={best_c2['ni_ne_ratio']:.2f}")
    print(f"  H={best_c2['H_ratio']:.2f}× | CH={best_c2['CH_ratio']:.1f}×")
    print()

    # Best C2H2
    best_c2h2 = max(results, key=lambda x: x['C2H2'])
    print("BEST C2H2:")
    print(f"  {best_c2h2['name']}")
    print(f"  C2H2={best_c2h2['C2H2']:.2e} cm⁻³")
    print(f"  C2={best_c2h2['C2_ratio']:.2f}× | f(x)={best_c2h2['objective']:.1f}")
    print(f"  Ni/Ne={best_c2h2['ni_ne_ratio']:.2f} | H={best_c2h2['H_ratio']:.2f}× | "
          f"CH={best_c2h2['CH_ratio']:.1f}×")
    print()

print("="*80)
print("FULL RESULTS TABLE")
print("="*80)
print()
print(f"{'Config':<40} | {'Ni/Ne':>7} | {'f(x)':>7} | {'H':>6} | {'CH':>6} | {'C2':>6} | {'C2H2':>10}")
print("-"*120)

for r in results:
    ni_ne_ok = "✓" if 2.0 <= r['ni_ne_ratio'] <= 6.0 else " "
    config_name = r['name'][:38]  # Truncate long names
    print(f"{config_name:<40} | {r['ni_ne_ratio']:>6.2f}{ni_ne_ok} | {r['objective']:>7.1f} | "
          f"{r['H_ratio']:>6.2f} | {r['CH_ratio']:>5.1f}× | {r['C2_ratio']:>6.2f} | {r['C2H2']:>10.2e}")

print()
print("="*80)

# Save detailed results
with open('c2h2_boost_comprehensive_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Detailed results saved to: c2h2_boost_comprehensive_results.json")
