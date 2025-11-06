#!/usr/bin/env python3
"""
Investigate why CH increased from 2.0× to 5.2× when we optimized H and C2
Compare CH pathways before and after stick_H reduction + C2H2 optimization
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
print("CH INCREASE INVESTIGATION")
print("="*80)
print()
print(f"Conditions: Te={TE_BEST} eV, E={E_BEST} V/cm, Ne={NE_BEST:.1e} cm⁻³")
print()

def run_simulation_with_analysis(rate_values, description):
    """Run simulation and return detailed CH pathway analysis."""
    params = params_base.copy()
    params['rate_values'] = rate_values.copy()

    # Apply conditions
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
    CH_idx = ode.species.index('CH')
    H_idx = ode.species.index('H')
    C2_idx = ode.species.index('C2')

    CH_density = y_final[CH_idx]
    H_density = y_final[H_idx]
    C2_density = y_final[C2_idx]

    # Analyze CH production/consumption
    CH_production = {}
    CH_consumption = {}

    for i, reaction in enumerate(ode.R):
        # Calculate reaction rate
        rate = k_dict[tags[i]]
        for idx, coeff in zip(ode.reactant_indices[i], ode.reactant_coeffs[i]):
            rate *= y_final[idx]**coeff

        # Net CH change in this reaction
        CH_net = reaction.products[CH_idx] - reaction.reactants[CH_idx]

        if CH_net > 0:
            # CH production
            CH_production[tags[i]] = rate * CH_net
        elif CH_net < 0:
            # CH consumption
            CH_consumption[tags[i]] = -rate * CH_net  # Store as positive

    return {
        'description': description,
        'CH_density': CH_density,
        'H_density': H_density,
        'C2_density': C2_density,
        'CH_production': CH_production,
        'CH_consumption': CH_consumption,
        'y_final': y_final,
        'species': ode.species,
    }

# Test 1: LOW CH configuration (before stick_H reduction)
print("[1/2] Running: LOW CH configuration (before optimization)")
low_ch_params = {
    'stick_C2_9_9': 1.25e3,
    'loss_C2_11_3': 100,
    'C2_H_CH_C_cm3_7_6': 8.0e-11,
    'stick_C2H2_9_11': 472,
    # stick_H at DEFAULT (not reduced)
}

result_low_ch = run_simulation_with_analysis(low_ch_params, "LOW CH (before stick_H reduction)")

if result_low_ch:
    print(f"  CH = {result_low_ch['CH_density']:.2e} cm⁻³ ({result_low_ch['CH_density']/TARGETS['CH']:.1f}×)")
    print(f"  H  = {result_low_ch['H_density']:.2e} cm⁻³ ({result_low_ch['H_density']/TARGETS['H']:.2f}×)")
    print(f"  C2 = {result_low_ch['C2_density']:.2e} cm⁻³ ({result_low_ch['C2_density']/TARGETS['C2']:.2f}×)")
else:
    print("  FAILED")
print()

# Test 2: HIGH CH configuration (after full optimization)
print("[2/2] Running: HIGH CH configuration (after optimization)")
high_ch_params = {
    'stick_C2_9_9': 1.25e3,
    'loss_C2_11_3': 100,
    'C2_H_CH_C_cm3_7_6': 8.0e-11,
    'stick_H_9_1': 3.89e2,  # REDUCED
    # C2H2 optimizations
    'stick_C2H2_9_11': 500,
    'loss_C2H2_11_19': 1000,
    'loss_C2H2Star_11_25': 100,
    'CH_CH2_C2H2_H_cm3_7_7': 1.20e-10,
    'CH2_CH2_C2H2_H2_cm3_7_15': 1.20e-10,
    'CH3_CH_C2H2_H2_cm3_7_16': 1.20e-10,
    'CH2_C_C2H2_cm3_7_17': 1.20e-10,
    'CH_C2H2_C3H2_H_cm3_7_22': 8.0e-11,
    'CH_C2H2_C3H_H2_cm3_7_27': 8.0e-11,
    'CH_C2H2_C2H_CH2_cm3_7_29': 8.0e-11,
}

result_high_ch = run_simulation_with_analysis(high_ch_params, "HIGH CH (after optimization)")

if result_high_ch:
    print(f"  CH = {result_high_ch['CH_density']:.2e} cm⁻³ ({result_high_ch['CH_density']/TARGETS['CH']:.1f}×)")
    print(f"  H  = {result_high_ch['H_density']:.2e} cm⁻³ ({result_high_ch['H_density']/TARGETS['H']:.2f}×)")
    print(f"  C2 = {result_high_ch['C2_density']:.2e} cm⁻³ ({result_high_ch['C2_density']/TARGETS['C2']:.2f}×)")
else:
    print("  FAILED")
print()

if not (result_low_ch and result_high_ch):
    print("ERROR: Simulations failed!")
    exit(1)

# Compare CH pathways
print("="*80)
print("CH PRODUCTION COMPARISON")
print("="*80)
print()

total_prod_low = sum(result_low_ch['CH_production'].values())
total_prod_high = sum(result_high_ch['CH_production'].values())

print(f"{'Reaction':<50} | {'Low CH':<15} | {'High CH':<15} | {'Change':<10}")
print("-"*95)

# Get all unique reactions
all_reactions = set(result_low_ch['CH_production'].keys()) | set(result_high_ch['CH_production'].keys())
sorted_reactions = sorted(all_reactions,
                         key=lambda x: result_high_ch['CH_production'].get(x, 0),
                         reverse=True)

for tag in sorted_reactions[:15]:
    rate_low = result_low_ch['CH_production'].get(tag, 0)
    rate_high = result_high_ch['CH_production'].get(tag, 0)

    if rate_low > 0:
        change = (rate_high / rate_low - 1) * 100
        change_str = f"{change:+.0f}%"
    else:
        change_str = "NEW"

    pct_low = 100.0 * rate_low / total_prod_low if total_prod_low > 0 else 0
    pct_high = 100.0 * rate_high / total_prod_high if total_prod_high > 0 else 0

    print(f"{tag[:48]:<50} | {rate_low:<15.2e} | {rate_high:<15.2e} | {change_str:<10}")

print()
print(f"{'TOTAL':<50} | {total_prod_low:<15.2e} | {total_prod_high:<15.2e} | "
      f"{(total_prod_high/total_prod_low - 1)*100:+.0f}%")
print()

print("="*80)
print("CH CONSUMPTION COMPARISON")
print("="*80)
print()

total_cons_low = sum(result_low_ch['CH_consumption'].values())
total_cons_high = sum(result_high_ch['CH_consumption'].values())

print(f"{'Reaction':<50} | {'Low CH':<15} | {'High CH':<15} | {'Change':<10}")
print("-"*95)

all_cons_reactions = set(result_low_ch['CH_consumption'].keys()) | set(result_high_ch['CH_consumption'].keys())
sorted_cons_reactions = sorted(all_cons_reactions,
                              key=lambda x: result_high_ch['CH_consumption'].get(x, 0),
                              reverse=True)

for tag in sorted_cons_reactions[:15]:
    rate_low = result_low_ch['CH_consumption'].get(tag, 0)
    rate_high = result_high_ch['CH_consumption'].get(tag, 0)

    if rate_low > 0:
        change = (rate_high / rate_low - 1) * 100
        change_str = f"{change:+.0f}%"
    else:
        change_str = "NEW"

    print(f"{tag[:48]:<50} | {rate_low:<15.2e} | {rate_high:<15.2e} | {change_str:<10}")

print()
print(f"{'TOTAL':<50} | {total_cons_low:<15.2e} | {total_cons_high:<15.2e} | "
      f"{(total_cons_high/total_cons_low - 1)*100:+.0f}%")
print()

print("="*80)
print("KEY FINDINGS")
print("="*80)
print()

ch_increase = (result_high_ch['CH_density'] / result_low_ch['CH_density'] - 1) * 100
prod_increase = (total_prod_high / total_prod_low - 1) * 100
cons_increase = (total_cons_high / total_cons_low - 1) * 100

print(f"CH density increased by: {ch_increase:+.0f}% ({result_low_ch['CH_density']/TARGETS['CH']:.1f}× → "
      f"{result_high_ch['CH_density']/TARGETS['CH']:.1f}×)")
print(f"CH production increased by: {prod_increase:+.0f}%")
print(f"CH consumption increased by: {cons_increase:+.0f}%")
print()

print("Net effect: Production increased MORE than consumption → higher steady-state CH")
print()

print("="*80)
print("PROPOSED STRATEGIES TO REDUCE CH")
print("="*80)
print()

# Identify top CH production pathways
top_prod = sorted(result_high_ch['CH_production'].items(), key=lambda x: x[1], reverse=True)[:5]

print("Target the TOP CH PRODUCTION pathways:")
for i, (tag, rate) in enumerate(top_prod, 1):
    pct = 100.0 * rate / total_prod_high
    print(f"{i}. {tag}")
    print(f"   Rate: {rate:.2e} cm⁻³/s ({pct:.1f}% of CH production)")

    # Check if tunable
    if tag in db:
        print(f"   Range: {db[tag].min:.2e} to {db[tag].max:.2e}")
        print(f"   → Can reduce to minimum to suppress CH production")
    else:
        print(f"   → Not directly tunable (complex reaction or fixed)")
    print()

print("="*80)

# Save detailed results
results_summary = {
    'low_ch': {
        'CH': result_low_ch['CH_density'],
        'H': result_low_ch['H_density'],
        'C2': result_low_ch['C2_density'],
        'CH_production_total': total_prod_low,
        'CH_consumption_total': total_cons_low,
    },
    'high_ch': {
        'CH': result_high_ch['CH_density'],
        'H': result_high_ch['H_density'],
        'C2': result_high_ch['C2_density'],
        'CH_production_total': total_prod_high,
        'CH_consumption_total': total_cons_high,
    },
    'changes': {
        'CH_increase_pct': ch_increase,
        'production_increase_pct': prod_increase,
        'consumption_increase_pct': cons_increase,
    }
}

with open('ch_investigation_results.json', 'w') as f:
    json.dump(results_summary, f, indent=2)

print("Results saved to: ch_investigation_results.json")
