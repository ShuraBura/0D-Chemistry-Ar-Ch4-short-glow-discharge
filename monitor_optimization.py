#!/usr/bin/env python3
"""
Monitor optimization progress and analyze C2 production mechanisms
"""

import json
import os
import glob
from pathlib import Path

# Targets
TARGETS = {'H': 2.52e14, 'CH': 1.0e9, 'C2': 5.6e11}

print("=" * 80)
print("OPTIMIZATION PROGRESS MONITOR")
print("=" * 80)
print()

# Find all result files
results_dir = 'optimization_results_charge_balanced_fixed_weights'
if not os.path.exists(results_dir):
    print(f"Results directory not found: {results_dir}")
    exit(0)

files = sorted(glob.glob(f'{results_dir}/best_f*.json'))

if not files:
    print("No results yet...")
    exit(0)

print(f"Found {len(files)} results so far")
print()

# Parse objectives from filenames
results = []
for f in files:
    fname = os.path.basename(f)
    obj_str = fname.replace('best_f', '').replace('.json', '')
    try:
        obj = float(obj_str)
        results.append((obj, f))
    except:
        continue

# Sort by objective (best first)
results.sort(key=lambda x: x[0])

# Show top 5 results
print("TOP 5 RESULTS:")
print("-" * 80)
for i, (obj, fpath) in enumerate(results[:5], 1):
    with open(fpath, 'r') as f:
        data = json.load(f)

    H = data['target_densities']['H']
    CH = data['target_densities']['CH']
    C2 = data['target_densities']['C2']
    Te = data['Te']
    Ne = data['Ne']
    E = data['E_field']
    Ni_Ne = data.get('Ni_over_Ne', 0)

    print(f"\n{i}. f(x) = {obj:.2f}")
    print(f"   Te={Te:.2f} eV, Ne={Ne:.2e}, E={E:.1f} V/cm, Ni/Ne={Ni_Ne:.2f}")
    print(f"   H:  {H:.2e} ({H/TARGETS['H']*100:5.1f}% of target)")
    print(f"   CH: {CH:.2e} ({CH/TARGETS['CH']*100:5.1f}% of target)")
    print(f"   C2: {C2:.2e} ({C2/TARGETS['C2']*100:5.1f}% of target)")

print()
print("=" * 80)
print("BEST RESULT ANALYSIS")
print("=" * 80)
print()

# Analyze best result in detail
best_obj, best_file = results[0]
with open(best_file, 'r') as f:
    best = json.load(f)

print(f"File: {os.path.basename(best_file)}")
print(f"Objective: {best_obj:.2f}")
print()

# Check C2-producing reaction rates
print("C2-PRODUCING REACTION RATES:")
print("-" * 80)

rate_values = best.get('rate_values', {})
c2_reactions = {
    'CH_CH_C2_H2_cm3_5_4': 'CH + CH → C2 + H2',
    'e_C2H2_C2_H2_cm3_1_16': 'e + C2H2 → C2 + H2',
    'C2HPlus_e_C2_H_cm3_6_18': 'C2H+ + e → C2 + H',
    'C_CH_C2_H_cm3_7_4': 'C + CH → C2 + H',
    'CH_C_C2_H_cm3_7_9': 'CH + C → C2 + H',
    'C2H_H_C2_H2_cm3_7_47': 'C2H + H → C2 + H2',
}

c2_rates_found = False
for rate_name, rxn_desc in c2_reactions.items():
    if rate_name in rate_values:
        k = rate_values[rate_name]
        print(f"{rxn_desc}")
        print(f"  k = {k:.2e} cm³/s")
        c2_rates_found = True

if not c2_rates_found:
    print("⚠ WARNING: No C2-producing reactions in rate_values!")
    print("This means they're using default values (not being tuned)")

print()
print("C2 LOSS RATES:")
print("-" * 80)

c2_loss = {
    'loss_C2_11_3': 'C2 volumetric loss',
    'stick_C2_9_9': 'C2 wall sticking',
}

for rate_name, desc in c2_loss.items():
    if rate_name in rate_values:
        k = rate_values[rate_name]
        print(f"{desc}")
        print(f"  k = {k:.2e} s⁻¹")

print()
print("PRECURSOR DENSITIES:")
print("-" * 80)

all_dens = best.get('all_densities', {})
precursors = {
    'CH': 'CH radicals',
    'C': 'C atoms',
    'C2H2': 'C2H2 (acetylene)',
    'C2H': 'C2H radicals',
    'C2HPlus': 'C2H+ ions',
    'e': 'Electrons',
}

for species, desc in precursors.items():
    if species in all_dens:
        n = all_dens[species]
        print(f"{species:8s} = {n:.2e} cm⁻³  ({desc})")

print()
print("=" * 80)
print("CONVERGENCE CHECK")
print("=" * 80)
print()

# Check if optimizer is converging
if len(results) >= 10:
    recent_10 = [obj for obj, _ in results[:10]]
    best_obj = min(recent_10)
    worst_in_top10 = max(recent_10)
    spread = worst_in_top10 - best_obj

    print(f"Best objective:        {best_obj:.2f}")
    print(f"Worst in top 10:       {worst_in_top10:.2f}")
    print(f"Spread in top 10:      {spread:.2f}")

    if spread < 10.0:
        print()
        print("✓ Optimizer appears to be converging (spread < 10)")
    else:
        print()
        print("⚠ Optimizer still exploring (large spread in top results)")
else:
    print(f"Only {len(results)} results so far - need more data for convergence check")

print()
