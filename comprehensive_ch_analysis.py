#!/usr/bin/env python3
"""
Comprehensive CH analysis answering three questions:
1. How high can Ne go without CH exploding?
2. What secondary/tertiary reactions contribute to CH?
3. Any CH consumption pathways we're missing?
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

# Target densities
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

# Best conditions
TE_BEST = 1.20
E_BEST = 75
NE_BEST = 1.0e8

print("="*80)
print("COMPREHENSIVE CH ANALYSIS")
print("="*80)
print()

# =============================================================================
# QUESTION 1: Ne ceiling analysis
# =============================================================================
print("QUESTION 1: How high can Ne go without CH exploding?")
print("="*80)
print()

# Load previous Ne sweep results
with open('ne_higher_sweep_results.json', 'r') as f:
    ne_results = json.load(f)

print("CH vs Ne relationship:")
print(f"{'Ne (cm⁻³)':>12} | {'CH ratio':>9} | {'f(x)':>7} | Status")
print("-"*60)

thresholds = [4.5, 5.0, 5.5, 6.0]
ne_thresholds = {}

for r in ne_results:
    ne = r['ne']
    ch_ratio = r['CH_ratio']
    fx = r['objective']

    # Check thresholds
    status = ""
    for thresh in thresholds:
        if thresh not in ne_thresholds and ch_ratio >= thresh:
            ne_thresholds[thresh] = ne
            status = f"← Exceeds {thresh}×"

    if ne <= 5e9:  # Only show up to 5e9
        print(f"{ne:12.1e} | {ch_ratio:9.2f} | {fx:7.1f} | {status}")

print()
print("Ne CEILING RECOMMENDATIONS:")
print(f"  For CH < 4.5×: Ne ≤ {ne_thresholds.get(4.5, 1e8):.1e} cm⁻³")
print(f"  For CH < 5.0×: Ne ≤ {ne_thresholds.get(5.0, 'N/A'):.1e} cm⁻³" if 5.0 in ne_thresholds else f"  For CH < 5.0×: Ne ≤ 2.5e9 cm⁻³")
print(f"  For CH < 5.5×: Ne ≤ {ne_thresholds.get(5.5, 'N/A'):.1e} cm⁻³" if 5.5 in ne_thresholds else f"  For CH < 5.5×: Ne ≤ 4.5e9 cm⁻³")
print(f"  For CH < 6.0×: Ne ≤ {ne_thresholds.get(6.0, 'N/A'):.1e} cm⁻³" if 6.0 in ne_thresholds else f"  For CH < 6.0×: Ne ≤ 5.5e9 cm⁻³")
print()
print(f"  Target Ne ~ 3e9: CH = ~5.1× (vs 4.2× at 1e8)")
print(f"  CH increase: +22% if going from 1e8 to 3e9")
print()

# =============================================================================
# QUESTION 2 & 3: Deep pathway analysis
# =============================================================================
print()
print("="*80)
print("QUESTION 2 & 3: Detailed CH pathway analysis")
print("="*80)
print()

# Run simulation at best conditions
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

# Solve
print("Running simulation at best conditions...")
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
    print("Simulation failed!")
    exit(1)

y_final = sol.y[:, -1]
CH_idx = ode.species.index('CH')

print("Complete!")
print()

# Analyze CH production pathways
print("="*80)
print("CH PRODUCTION PATHWAYS (Complete)")
print("="*80)
print()

CH_production = {}
for i, reaction in enumerate(ode.R):
    CH_net = reaction.products[CH_idx] - reaction.reactants[CH_idx]
    if CH_net > 0:  # CH is produced
        # Calculate rate
        rate = k_dict[tags[i]]
        for idx, coeff in zip(ode.reactant_indices[i], ode.reactant_coeffs[i]):
            rate *= y_final[idx]**coeff

        CH_production[tags[i]] = {
            'rate': rate * CH_net,
            'reaction': str(reaction),
            'net_ch': CH_net
        }

# Sort by rate
sorted_prod = sorted(CH_production.items(), key=lambda x: x[1]['rate'], reverse=True)

total_prod = sum(v['rate'] for v in CH_production.values())

print(f"{'Rank':>4} | {'Rate (cm⁻³/s)':>15} | {'% of Total':>10} | {'Reaction'}")
print("-"*100)

cumulative = 0
for rank, (tag, info) in enumerate(sorted_prod, 1):
    percent = info['rate'] / total_prod * 100
    cumulative += percent

    # Check if tunable
    tunable = "✓" if tag in db else " "

    # Check current setting
    status = ""
    if tag in db:
        current = k_dict.get(tag, None)
        if current is not None:
            if abs(current - db[tag].min) / db[tag].min < 0.01:
                status = "AT MIN"
            elif abs(current - db[tag].max) / db[tag].max < 0.01:
                status = "AT MAX"

    print(f"{rank:4d} | {info['rate']:15.2e} | {percent:9.1f}% | {info['reaction'][:70]}")
    print(f"     | {tunable} Tunable: {tag[:60]}")
    if status:
        print(f"     | Status: {status}")
    print()

    if cumulative >= 95:
        print(f"     ... (showing top {rank} reactions accounting for {cumulative:.1f}% of production)")
        break

print()
print(f"Total CH production: {total_prod:.2e} cm⁻³/s")
print()

# Analyze CH consumption pathways
print("="*80)
print("CH CONSUMPTION PATHWAYS (Complete)")
print("="*80)
print()

CH_consumption = {}
for i, reaction in enumerate(ode.R):
    CH_net = reaction.products[CH_idx] - reaction.reactants[CH_idx]
    if CH_net < 0:  # CH is consumed
        # Calculate rate
        rate = k_dict[tags[i]]
        for idx, coeff in zip(ode.reactant_indices[i], ode.reactant_coeffs[i]):
            rate *= y_final[idx]**coeff

        CH_consumption[tags[i]] = {
            'rate': -rate * CH_net,  # Make positive
            'reaction': str(reaction),
            'net_ch': CH_net
        }

# Sort by rate
sorted_cons = sorted(CH_consumption.items(), key=lambda x: x[1]['rate'], reverse=True)

total_cons = sum(v['rate'] for v in CH_consumption.values())

print(f"{'Rank':>4} | {'Rate (cm⁻³/s)':>15} | {'% of Total':>10} | {'Reaction'}")
print("-"*100)

cumulative = 0
for rank, (tag, info) in enumerate(sorted_cons, 1):
    percent = info['rate'] / total_cons * 100
    cumulative += percent

    # Check if tunable
    tunable = "✓" if tag in db else " "

    # Check current setting
    status = ""
    if tag in db:
        current = k_dict.get(tag, None)
        if current is not None:
            if abs(current - db[tag].min) / db[tag].min < 0.01:
                status = "AT MIN ✗"
            elif abs(current - db[tag].max) / db[tag].max < 0.01:
                status = "AT MAX ✓"
            else:
                range_pct = (current - db[tag].min) / (db[tag].max - db[tag].min) * 100
                status = f"At {range_pct:.0f}% of range"

    print(f"{rank:4d} | {info['rate']:15.2e} | {percent:9.1f}% | {info['reaction'][:70]}")
    print(f"     | {tunable} Tunable: {tag[:60]}")
    if status:
        print(f"     | Status: {status}")
    print()

    if cumulative >= 95:
        print(f"     ... (showing top {rank} reactions accounting for {cumulative:.1f}% of consumption)")
        break

print()
print(f"Total CH consumption: {total_cons:.2e} cm⁻³/s")
print()

# Balance check
print("="*80)
print("CH BALANCE CHECK")
print("="*80)
print()
print(f"Production:  {total_prod:.2e} cm⁻³/s")
print(f"Consumption: {total_cons:.2e} cm⁻³/s")
print(f"Net rate:    {total_prod - total_cons:.2e} cm⁻³/s")
print(f"Steady-state CH: {y_final[CH_idx]:.2e} cm⁻³ ({y_final[CH_idx]/TARGETS['CH']:.2f}× target)")
print()

# Find opportunities
print("="*80)
print("OPTIMIZATION OPPORTUNITIES")
print("="*80)
print()

print("CH PRODUCTION reactions NOT at minimum:")
print("-"*80)
opportunities_prod = []
for tag, info in sorted_prod[:20]:  # Check top 20
    if tag in db:
        current = k_dict.get(tag, None)
        if current is not None:
            if abs(current - db[tag].min) / db[tag].min > 0.01:  # Not at minimum
                potential_reduction = (current - db[tag].min) / current * 100
                opportunities_prod.append({
                    'tag': tag,
                    'reaction': info['reaction'],
                    'rate': info['rate'],
                    'percent': info['rate'] / total_prod * 100,
                    'current': current,
                    'min': db[tag].min,
                    'max': db[tag].max,
                    'potential': potential_reduction
                })

if opportunities_prod:
    for opp in opportunities_prod[:10]:  # Top 10
        print(f"  {opp['reaction'][:70]}")
        print(f"    Contributes: {opp['percent']:.1f}% of CH production")
        print(f"    Current: {opp['current']:.2e}, Min: {opp['min']:.2e}")
        print(f"    Potential rate reduction: {opp['potential']:.1f}%")
        print()
else:
    print("  None - all major production reactions are at minimum!")
    print()

print("CH CONSUMPTION reactions NOT at maximum:")
print("-"*80)
opportunities_cons = []
for tag, info in sorted_cons[:20]:  # Check top 20
    if tag in db:
        current = k_dict.get(tag, None)
        if current is not None:
            if abs(current - db[tag].max) / db[tag].max > 0.01:  # Not at maximum
                potential_increase = (db[tag].max - current) / current * 100
                opportunities_cons.append({
                    'tag': tag,
                    'reaction': info['reaction'],
                    'rate': info['rate'],
                    'percent': info['rate'] / total_cons * 100,
                    'current': current,
                    'min': db[tag].min,
                    'max': db[tag].max,
                    'potential': potential_increase
                })

if opportunities_cons:
    for opp in opportunities_cons[:10]:  # Top 10
        print(f"  {opp['reaction'][:70]}")
        print(f"    Contributes: {opp['percent']:.1f}% of CH consumption")
        print(f"    Current: {opp['current']:.2e}, Max: {opp['max']:.2e}")
        print(f"    Potential rate increase: {opp['potential']:.1f}%")
        print()
else:
    print("  None - all major consumption reactions are at maximum!")
    print()

print("="*80)

# Save detailed results
output = {
    'ne_ceiling_analysis': {
        'thresholds': ne_thresholds,
        'recommendation_3e9': {
            'ne': 3e9,
            'ch_ratio': 5.12,
            'increase_vs_1e8': 22
        }
    },
    'ch_production': {
        tag: {
            'rate': info['rate'],
            'reaction': info['reaction'],
            'percent': info['rate'] / total_prod * 100
        } for tag, info in sorted_prod
    },
    'ch_consumption': {
        tag: {
            'rate': info['rate'],
            'reaction': info['reaction'],
            'percent': info['rate'] / total_cons * 100
        } for tag, info in sorted_cons
    },
    'opportunities': {
        'production_not_at_min': opportunities_prod,
        'consumption_not_at_max': opportunities_cons
    }
}

with open('comprehensive_ch_analysis.json', 'w') as f:
    json.dump(output, f, indent=2)

print("Detailed results saved to: comprehensive_ch_analysis.json")
