#!/usr/bin/env python3
"""
Analyze CH consumption pathways at current best conditions
Check if consumption reactions are maximized
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
print("CH CONSUMPTION ANALYSIS")
print("="*80)
print()
print(f"Conditions: Te={TE_BEST} eV, E={E_BEST} V/cm, Ne={NE_BEST:.1e} cm⁻³")
print()

def run_simulation_with_analysis(description):
    """Run simulation and return CH pathway analysis."""
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply ALL OPTIMIZATIONS
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
    CH_density = y_final[CH_idx]

    # Analyze CH consumption
    CH_consumption = {}

    for i, reaction in enumerate(ode.R):
        # Calculate reaction rate
        rate = k_dict[tags[i]]
        for idx, coeff in zip(ode.reactant_indices[i], ode.reactant_coeffs[i]):
            rate *= y_final[idx]**coeff

        # Net CH change in this reaction
        CH_net = reaction.products[CH_idx] - reaction.reactants[CH_idx]

        if CH_net < 0:
            # CH consumption
            CH_consumption[tags[i]] = -rate * CH_net  # Store as positive

    return {
        'description': description,
        'CH_density': CH_density,
        'CH_consumption': CH_consumption,
        'k_dict': k_dict,
    }

# Run current best
result = run_simulation_with_analysis("Current best")

if result is None:
    print("Simulation failed!")
    exit(1)

CH_density = result['CH_density']
CH_consumption = result['CH_consumption']
k_dict = result['k_dict']

print(f"CH density: {CH_density:.2e} cm⁻³ ({CH_density/TARGETS['CH']:.2f}× target)")
print()

# Analyze top CH consumption pathways
print("="*80)
print("TOP 15 CH CONSUMPTION PATHWAYS")
print("="*80)
print()

total_consumption = sum(CH_consumption.values())
sorted_cons = sorted(CH_consumption.items(), key=lambda x: x[1], reverse=True)

print(f"{'Reaction':<50} | {'Rate (cm⁻³/s)':<15} | {'%':<8} | {'Status':<20}")
print("-"*110)

tunable_consumption = []

for tag, rate in sorted_cons[:15]:
    pct = 100.0 * rate / total_consumption

    # Check if this is a tunable rate
    if tag in db:
        current_val = k_dict.get(tag, None)
        if current_val is not None:
            lit_min = db[tag].min
            lit_max = db[tag].max

            # Check if at maximum
            if abs(current_val - lit_max) / lit_max < 0.01:
                status = "AT MAXIMUM ✓"
            elif abs(current_val - lit_min) / lit_min < 0.01:
                status = "AT MINIMUM ✗"
                tunable_consumption.append({
                    'tag': tag,
                    'current': current_val,
                    'min': lit_min,
                    'max': lit_max,
                    'rate': rate,
                    'pct': pct
                })
            else:
                status = f"Mid ({current_val:.2e})"
                tunable_consumption.append({
                    'tag': tag,
                    'current': current_val,
                    'min': lit_min,
                    'max': lit_max,
                    'rate': rate,
                    'pct': pct
                })
        else:
            status = "Not in k_dict"
    else:
        status = "Not tunable"

    print(f"{tag[:48]:<50} | {rate:<15.2e} | {pct:>6.1f}% | {status:<20}")

print()
print(f"Total CH consumption: {total_consumption:.2e} cm⁻³/s")
print()

print("="*80)
print("TUNABLE CH CONSUMPTION REACTIONS (NOT AT MAXIMUM)")
print("="*80)
print()

if tunable_consumption:
    print(f"Found {len(tunable_consumption)} reactions that can be increased!")
    print()
    print(f"{'Reaction':<50} | {'Current':<12} | {'Max':<12} | {'Increase':<10} | {'% of CH cons':<12}")
    print("-"*110)

    for item in tunable_consumption:
        increase = (item['max'] / item['current'] - 1) * 100
        print(f"{item['tag'][:48]:<50} | {item['current']:<12.2e} | {item['max']:<12.2e} | "
              f"{increase:>8.1f}% | {item['pct']:>10.1f}%")

    print()
    print("RECOMMENDATION: Maximize these reactions to increase CH consumption!")
    print()

    # Save for testing
    with open('ch_consumption_to_maximize.json', 'w') as f:
        json.dump(tunable_consumption, f, indent=2)
    print("Tunable reactions saved to: ch_consumption_to_maximize.json")

else:
    print("All major CH consumption reactions are already at maximum!")
    print("No further optimization possible from consumption side.")

print()
print("="*80)
