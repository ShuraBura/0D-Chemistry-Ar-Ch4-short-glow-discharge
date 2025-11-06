#!/usr/bin/env python3
"""
Analyze H production and consumption pathways at best conditions
Then test minimizing H + C2 → CH + C reaction
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
print("H PRODUCTION/CONSUMPTION ANALYSIS")
print("="*80)
print()
print(f"Conditions: Te={TE_BEST} eV, E={E_BEST} V/cm, Ne={NE_BEST:.1e} cm⁻³")
print()

def run_simulation(rate_values_override=None):
    """Run simulation and return detailed results."""
    params = params_base.copy()
    params['rate_values'] = params_base['rate_values'].copy()

    # Apply BEST chemical parameters
    params['rate_values']['stick_C2_9_9'] = 1.25e3
    params['rate_values']['loss_C2_11_3'] = 100
    params['rate_values']['C2_H_CH_C_cm3_7_6'] = 8.0e-11
    params['rate_values']['stick_C2H2_9_11'] = 472

    # Apply any overrides
    if rate_values_override:
        params['rate_values'].update(rate_values_override)

    # Apply BEST conditions
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

    # Calculate total ion density
    ni_total = 0
    for i_sp, species_name in enumerate(ode.species):
        if species_name.endswith('Plus'):
            ni_total += y_final[i_sp]
        elif species_name.endswith('Minus'):
            ni_total -= y_final[i_sp]

    # Analyze H production/consumption
    H_idx = ode.species.index('H')
    H_density = y_final[H_idx]

    H_production = {}
    H_consumption = {}

    for i, reaction in enumerate(ode.R):
        # Calculate reaction rate
        rate = k_dict[tags[i]]
        for idx, coeff in zip(ode.reactant_indices[i], ode.reactant_coeffs[i]):
            rate *= y_final[idx]**coeff

        # Net H change in this reaction
        H_net = reaction.products[H_idx] - reaction.reactants[H_idx]

        if H_net > 0:
            # H production
            H_production[tags[i]] = rate * H_net
        elif H_net < 0:
            # H consumption
            H_consumption[tags[i]] = -rate * H_net  # Store as positive

    # Add drift term
    H_production['H_drift_gain'] = ode.H_drift_gain

    return {
        'y_final': y_final,
        'species': ode.species,
        'H_density': H_density,
        'H_production': H_production,
        'H_consumption': H_consumption,
        'ni_total': ni_total,
        'ne': y_final[ode.species.index('e')],
        'k_dict': k_dict,
        'tags': tags,
    }

# Run baseline simulation
print("Running baseline simulation...")
result = run_simulation()

if result is None:
    print("Simulation failed!")
    exit(1)

H_idx = result['species'].index('H')
CH_idx = result['species'].index('CH')
C2_idx = result['species'].index('C2')

H_density = result['H_density']
CH_density = result['y_final'][CH_idx]
C2_density = result['y_final'][C2_idx]

print()
print("BASELINE RESULTS:")
print(f"  H  = {H_density:.2e} cm⁻³ ({H_density/TARGETS['H']:.2f}× target)")
print(f"  CH = {CH_density:.2e} cm⁻³ ({CH_density/TARGETS['CH']:.1f}× target)")
print(f"  C2 = {C2_density:.2e} cm⁻³ ({C2_density/TARGETS['C2']:.2f}× target)")
print(f"  Ni/Ne = {result['ni_total']/result['ne']:.2f}")
print()

# Analyze H production
print("="*80)
print("H PRODUCTION PATHWAYS (Top 10)")
print("="*80)

total_production = sum(result['H_production'].values())
sorted_prod = sorted(result['H_production'].items(), key=lambda x: x[1], reverse=True)

print(f"{'Reaction':<50} | {'Rate (cm⁻³/s)':<15} | {'%':<8}")
print("-"*80)
for tag, rate in sorted_prod[:10]:
    pct = 100.0 * rate / total_production
    print(f"{tag:<50} | {rate:<15.2e} | {pct:>6.1f}%")

print(f"\nTotal H production: {total_production:.2e} cm⁻³/s")
print()

# Analyze H consumption
print("="*80)
print("H CONSUMPTION PATHWAYS (Top 10)")
print("="*80)

total_consumption = sum(result['H_consumption'].values())
sorted_cons = sorted(result['H_consumption'].items(), key=lambda x: x[1], reverse=True)

print(f"{'Reaction':<50} | {'Rate (cm⁻³/s)':<15} | {'%':<8}")
print("-"*80)
for tag, rate in sorted_cons[:10]:
    pct = 100.0 * rate / total_consumption
    print(f"{tag:<50} | {rate:<15.2e} | {pct:>6.1f}%")

print(f"\nTotal H consumption: {total_consumption:.2e} cm⁻³/s")
print()

# Check balance
print(f"Net balance: {total_production - total_consumption:.2e} cm⁻³/s")
print(f"  (Should be ~0 at steady state)")
print()

# Find C2_H_CH_C contribution
c2_h_ch_c_rate = result['H_consumption'].get('C2_H_CH_C_cm3_7_6', 0.0)
c2_h_ch_c_pct = 100.0 * c2_h_ch_c_rate / total_consumption if total_consumption > 0 else 0

print("="*80)
print(f"H + C2 → CH + C pathway:")
print(f"  Current rate coefficient: {result['k_dict']['C2_H_CH_C_cm3_7_6']:.2e} cm³/s")
print(f"  Literature range: {db['C2_H_CH_C_cm3_7_6'].min:.2e} to {db['C2_H_CH_C_cm3_7_6'].max:.2e}")
print(f"  Currently at: MINIMUM")
print(f"  Consumption rate: {c2_h_ch_c_rate:.2e} cm⁻³/s ({c2_h_ch_c_pct:.1f}% of total H consumption)")
print("="*80)
print()

# Now test with further reduction (already at minimum, so this won't help)
print("NOTE: C2_H_CH_C is already at literature MINIMUM (8.0e-11 cm³/s)")
print("Cannot reduce further within literature constraints.")
print()

# But let's verify what happens if we hypothetically could reduce it
print("="*80)
print("HYPOTHETICAL TEST: Reduce C2_H_CH_C to 50% of minimum")
print("(This is OUTSIDE literature range, for analysis only)")
print("="*80)
print()

result_reduced = run_simulation({'C2_H_CH_C_cm3_7_6': 4.0e-11})

if result_reduced:
    H_new = result_reduced['H_density']
    CH_new = result_reduced['y_final'][CH_idx]
    C2_new = result_reduced['y_final'][C2_idx]
    ni_ne_new = result_reduced['ni_total'] / result_reduced['ne']

    print(f"Results with C2_H_CH_C = 4.0e-11 (50% of minimum):")
    print(f"  H  = {H_new:.2e} cm⁻³ ({H_new/TARGETS['H']:.2f}×) | Change: {(H_new/H_density - 1)*100:+.1f}%")
    print(f"  CH = {CH_new:.2e} cm⁻³ ({CH_new/TARGETS['CH']:.1f}×) | Change: {(CH_new/CH_density - 1)*100:+.1f}%")
    print(f"  C2 = {C2_new:.2e} cm⁻³ ({C2_new/TARGETS['C2']:.2f}×) | Change: {(C2_new/C2_density - 1)*100:+.1f}%")
    print(f"  Ni/Ne = {ni_ne_new:.2f}")
    print()
    print("Impact:")
    print(f"  - H would increase by {(H_new/H_density - 1)*100:.1f}%")
    print(f"  - C2 would increase by {(C2_new/C2_density - 1)*100:.1f}%")
    print(f"  - CH would decrease by {(CH_new/CH_density - 1)*100:.1f}%")
else:
    print("Hypothetical simulation failed.")

print()
print("="*80)
print("CONCLUSION")
print("="*80)
print()
print("The H + C2 → CH + C reaction is already at the MINIMUM literature value.")
print(f"It accounts for {c2_h_ch_c_pct:.1f}% of H consumption.")
print("No further reduction possible within validated chemistry constraints.")
