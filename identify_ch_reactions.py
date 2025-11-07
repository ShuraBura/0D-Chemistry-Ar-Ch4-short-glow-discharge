#!/usr/bin/env python3
"""
Identify the specific CH reactions from the analysis
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

# Best conditions with all optimizations
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

params['Te'] = 1.20
params['E_field'] = 75
params['ne'] = 1.0e8

k_dict = define_rates_tunable(params)
for name, val in params.get('rate_values', {}).items():
    if name in k_dict and name in db:
        k_dict[name] = np.clip(val, db[name].min, db[name].max)

params['k'] = k_dict
reactions, tags = build_reactions(params)
params['R'] = reactions
params['tags'] = tags

ode = PlasmaODE(params)

y0 = np.zeros(ode.ns)
y0[ode.species.index('e')] = params['ne']
y0[ode.species.index('Ar')] = 1.29e16
y0[ode.species.index('CH4')] = 1.29e15

sol = solve_ivp(ode, t_span=[0, 1e-3], y0=y0, method='BDF', rtol=1e-6, atol=1e-10, max_step=1e-5)
y_final = sol.y[:, -1]
CH_idx = ode.species.index('CH')

print("="*80)
print("TOP CH REACTIONS WITH DETAILS")
print("="*80)
print()

# CH Production
print("TOP 5 CH PRODUCTION:")
print("-"*80)
CH_production = {}
for i, reaction in enumerate(ode.R):
    CH_net = reaction.products[CH_idx] - reaction.reactants[CH_idx]
    if CH_net > 0:
        rate = k_dict[tags[i]]
        for idx, coeff in zip(ode.reactant_indices[i], ode.reactant_coeffs[i]):
            rate *= y_final[idx]**coeff
        CH_production[tags[i]] = (rate * CH_net, reaction)

sorted_prod = sorted(CH_production.items(), key=lambda x: x[1][0], reverse=True)
total_prod = sum(v[0] for v in CH_production.values())

for rank, (tag, (rate, rxn)) in enumerate(sorted_prod[:5], 1):
    percent = rate / total_prod * 100

    # Get reaction string - just use str(rxn)
    print(f"{rank}. {str(rxn)}")
    print(f"   Tag: {tag}")
    print(f"   Rate: {rate:.2e} cm⁻³/s ({percent:.1f}% of total)")

    if tag in db:
        current = k_dict.get(tag, None)
        if current:
            print(f"   Current: {current:.2e}, Min: {db[tag].min:.2e}, Max: {db[tag].max:.2e}")
            if abs(current - db[tag].min) / db[tag].min < 0.01:
                print(f"   Status: AT MINIMUM ✓")
            elif abs(current - db[tag].max) / db[tag].max < 0.01:
                print(f"   Status: AT MAXIMUM")
    print()

# CH Consumption
print("="*80)
print("TOP 10 CH CONSUMPTION:")
print("-"*80)
CH_consumption = {}
for i, reaction in enumerate(ode.R):
    CH_net = reaction.products[CH_idx] - reaction.reactants[CH_idx]
    if CH_net < 0:
        rate = k_dict[tags[i]]
        for idx, coeff in zip(ode.reactant_indices[i], ode.reactant_coeffs[i]):
            rate *= y_final[idx]**coeff
        CH_consumption[tags[i]] = (-rate * CH_net, reaction)

sorted_cons = sorted(CH_consumption.items(), key=lambda x: x[1][0], reverse=True)
total_cons = sum(v[0] for v in CH_consumption.values())

for rank, (tag, (rate, rxn)) in enumerate(sorted_cons[:10], 1):
    percent = rate / total_cons * 100

    # Get reaction string - just use str(rxn)
    print(f"{rank}. {str(rxn)}")
    print(f"   Tag: {tag}")
    print(f"   Rate: {rate:.2e} cm⁻³/s ({percent:.1f}% of total)")

    if tag in db:
        current = k_dict.get(tag, None)
        if current:
            print(f"   Current: {current:.2e}, Min: {db[tag].min:.2e}, Max: {db[tag].max:.2e}")
            if abs(current - db[tag].min) / db[tag].min < 0.01:
                print(f"   Status: AT MINIMUM ✗ (should maximize!)")
            elif abs(current - db[tag].max) / db[tag].max < 0.01:
                print(f"   Status: AT MAXIMUM ✓")
            else:
                # Check if we can increase
                if db[tag].max > current:
                    potential = (db[tag].max / current - 1) * 100
                    print(f"   Status: Can increase by {potential:.0f}%! ← OPPORTUNITY")
                else:
                    print(f"   Status: Already above lit max (constrained)")
    print()

print("="*80)
