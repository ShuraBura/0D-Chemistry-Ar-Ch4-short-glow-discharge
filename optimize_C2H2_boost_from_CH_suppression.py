#!/usr/bin/env python3
"""
C2H2-BOOST OPTIMIZATION FROM CH-SUPPRESSION BASELINE

Starting from CH-suppression Eval 1412 (f=63.7):
- CH: 0.99×✓ (perfect!)
- C2: 0.06× ✗ (way too low)
- H: 0.19× ✗ (too low)

KEY INSIGHT: C2H2 is 97% lost to walls, only 2.5% converts to C2!
- C2H2 + H → C2 pathway is throttled by wall losses
- C2H2: 3.40e11 (0.11× of "target")

Strategy - Boost C2H2 to increase C2:
1. START from CH-suppression solution (CH wall sticking works!)
2. REDUCE C2H2 wall losses (reduce stick_C2H2, loss_C2H2)
3. INCREASE C2H2 density → more C2H2+H→C2 conversion
4. BOOST H density (needed for C2H2→C2 pathway)
5. MAINTAIN CH suppression (keep high CH wall sticking)

New Weights:
- H: 40 (increased from 15 - need more H for C2H2→C2)
- CH: 80 (maintain MAX - keep CH suppressed)
- C2: 60 (maintain - primary target)
- C2H2: 20 (increased from 1 - boost C2H2 to feed C2)

Goal: H, CH, C2 all within 0.6-1.4×, charge balance <100%
"""

import numpy as np
from scipy.optimize import differential_evolution
import json
import os
from datetime import datetime

# Real targets: H, CH, C2 must be within 0.6-1.4× range
TARGETS = {
    'H': 5.18e13,     # ← Must be 0.6-1.4× (3.1e13 to 7.3e13)
    'CH': 1.0e9,      # ← Must be 0.6-1.4× (6e8 to 1.4e9)
    'C2': 1.3e11,     # ← Must be 0.6-1.4× (7.8e10 to 1.8e11)
    'C2H2': 3.0e12,   # ← Boost this to feed C2 production
}

# Create results directory
os.makedirs('optimization_results_C2H2_boost', exist_ok=True)

# Load starting solution (CH-suppression Eval 1412)
with open('optimization_results_CH_suppression/best_f63.7_Hdrift1.05e+17.json', 'r') as f:
    baseline = json.load(f)

print("="*80)
print(" C2H2-BOOST OPTIMIZATION FROM CH-SUPPRESSION BASELINE")
print("="*80)
print()
print("Baseline solution (CH-suppression Eval 1412, f=63.7):")
print(f"  H: 0.19× ✗ (need 0.6-1.4×)")
print(f"  CH: 0.99×✓ (PERFECT! Keep this!)")
print(f"  C2: 0.06× ✗ (need 0.6-1.4×)")
print(f"  Charge: 8.62%✓")
print()
print("Problem: C2H2 is 97% lost to walls, only 2.5% converts to C2")
print("  C2H2: 3.40e11 cm⁻³ (0.11×)")
print("  C2H2 + H → C2: 33.8 THz (only 2.5% of C2H2 loss)")
print()
print("Strategy: Reduce C2H2 wall losses → boost C2H2 → more C2H2+H→C2 conversion")
print("          Also boost H (needed for C2H2→C2 pathway)")
print("          Maintain CH suppression (keep high CH wall sticking)")
print()
print("New weights: H:40 (↑), CH:80 (=MAX), C2:60 (=), C2H2:20 (↑↑)")
print("Charge penalty: 20× (doubled)")
print()
print("Goal: H, CH, C2 all within 0.6-1.4× range, charge balance <100%")
print()

# Parse chemistry file
import re

def parse_chemistry_file(filename):
    """Parse chemistry.txt to extract reaction info"""
    reactions = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            reaction_id = parts[0]
            try:
                rate_value = float(parts[1])
                rate_type = parts[2] if len(parts) > 2 else 'constant'
                reactions[reaction_id] = {
                    'value': rate_value,
                    'type': rate_type,
                    'line': line
                }
            except ValueError:
                continue

    return reactions

# Load chemistry
reactions = parse_chemistry_file('chemistry.txt')

print("Selecting tunable rates...")

# Select key rates to tune
# Focus on C2H2 losses, C2 production, CH control
tunable_rates = {
    # C2H2 wall losses - REDUCE these to boost C2H2
    'loss_C2H2_11_19': (100.0, 1000.0),      # Reduce from baseline
    'stick_C2H2_9_11': (100.0, 1000.0),      # Reduce from baseline

    # C2 losses - control C2 balance
    'loss_C2_11_3': (1e-4, 2000.0),
    'stick_C2_9_9': (1000.0, 6000.0),        # Allow some variation

    # CH wall losses - MAINTAIN high to keep CH suppressed
    'stick_CH_9_3': (3000.0, 5000.0),        # Keep high
    'loss_CH_11_9': (3000.0, 6000.0),        # Keep high

    # CH production - keep suppressed
    'e_CH4_CH_H2_H_vib_cm3_1_3': (1e-11, 1e-10),
    'e_CH4_CH_H_H2_cm3_1_11': (1e-11, 1e-10),
    'C2_H_CH_C_cm3_7_6': (1e-12, 1e-11),     # C2+H→CH pathway

    # C2H2 → C2 pathway - maximize this!
    'C2H2_H_C2_H2_H_cm3_7_50': (1e-12, 1e-10),  # Boost this reaction

    # H production/loss
    'stick_H_9_1': (389.0, 3890.0),
    'loss_H_11_7': (1e-4, 1000.0),

    # Electron-impact reactions
    'e_CH4_CH3_HMinus_cm3_8_1': (1e-15, 1e-12),
    'e_C2H2_C2_H2_cm3_1_16': (1e-12, 1e-10),   # C2H2→C2 electron impact

    # Ion-molecule reactions
    'e_CH4Plus_CH2_H2_cm3_6_9': (1e-7, 3e-6),
    'e_CH4Plus_CH_H2_H_cm3_6_11': (1e-7, 1.2e-6),

    # Metastable reactions
    'ArStar_CH4_CH3_H_cm3_3_1': (5e-11, 5e-10),

    # Higher hydrocarbon losses
    'loss_C2H2Star_11_25': (100.0, 1e4),
    'loss_C3H4_11_18': (100.0, 1500.0),
    'loss_C3H_11_10': (100.0, 1000.0),
    'loss_C3_11_14': (10.0, 500.0),
    'loss_C2H_11_12': (100.0, 1000.0),
    'loss_C4H_11_15': (50.0, 500.0),
    'loss_C3H3_11_13': (100.0, 1000.0),
    'loss_C4H2_11_11': (100.0, 1000.0),
    'loss_C3H6_11_16': (100.0, 1000.0),

    # More wall sticking
    'stick_C_9_10': (1000.0, 6000.0),
    'stick_CH3_9_2': (3000.0, 6000.0),
    'stick_C3H_9_19': (500.0, 2000.0),
    'stick_C4H_9_24': (500.0, 2000.0),
    'stick_C2H5_9_17': (500.0, 2000.0),
    'stick_C2H4_9_12': (300.0, 1000.0),
    'stick_C3H4_9_22': (500.0, 1500.0),
    'stick_C3H3_9_21': (500.0, 1500.0),
    'stick_C3_9_23': (500.0, 1500.0),
    'stick_H2_9_16': (300.0, 1000.0),
    'stick_C4H2_9_20': (500.0, 2000.0),
    'stick_C2H6_9_14': (500.0, 2000.0),

    # Electron impact ionization
    'e_Ar_ArPlus_cm3_2_3': (1e-12, 1e-10),
    'e_CH4_CH4Plus_cm3_2_2': (1e-12, 1e-10),
}

print(f"  Selected {len(tunable_rates)} key rates")
print()
print("Selected rates (top 10):")
print("-"*80)
for i, (rate_id, bounds) in enumerate(list(tunable_rates.items())[:10], 1):
    ratio = bounds[1] / bounds[0]
    print(f"{i:2d}. {rate_id:45s} [{bounds[0]:.2e}, {bounds[1]:.2e}] ({ratio:.1f}x)")
print(f"    ... and {len(tunable_rates)-10} more")
print()

# Fixed parameters from baseline
fixed_params = {
    'pressure_torr': 0.01,
    'power_W': 25.0,
    'volume_cm3': 400.0,
    'area_cm2': 400.0,
    'Ar_fraction': 0.85,
    'CH4_fraction': 0.15,
}

# Optimization parameters
param_names = list(tunable_rates.keys())
bounds = []
for rate_id in param_names:
    bounds.append(tunable_rates[rate_id])

# Add global parameters
bounds.append((100, 400))      # E_field (V/cm)
bounds.append((0.7, 8.0))      # Te (eV)
bounds.append((1e8, 5e9))      # Ne (cm^-3)
bounds.append((1e16, 1e18))    # H_drift_gain (cm^-3/s) - boost H source

print(f" Optimization parameters:")
print(f"  Tunable rates: {len(tunable_rates)}")
print(f"  E field: [100, 400] V/cm")
print(f"  Te: [0.7, 8] eV")
print(f"  Ne: [1e8, 5e9] cm⁻³")
print(f"  H_drift_gain: [1e16, 1e18] cm⁻³/s")
print(f"  Total parameters: {len(bounds)}")
print()

print(f" Real Targets (must be 0.6-1.4×):")
print(f"  H:    {TARGETS['H']:.2e} cm^-3  (weight: 40.0) ← BOOSTED (need for C2H2→C2)")
print(f"  CH:   {TARGETS['CH']:.2e} cm^-3  (weight: 80.0) ← MAXIMUM! (maintain suppression)")
print(f"  C2:   {TARGETS['C2']:.2e} cm^-3  (weight: 60.0) → primary target")
print(f"  C2H2: {TARGETS['C2H2']:.2e} cm^-3  (weight: 20.0) ← BOOSTED (feed C2)")
print()
print(f"  + Charge balance (weight: 20.0) ← DOUBLED")
print(f"  + C2+H→CH flux penalty (weight: 2.0) ← Maintain CH suppression")
print()

# Best solution tracking
best_f = float('inf')
best_x = None
eval_count = [0]

def run_plasma_model(params):
    """Run the 0D plasma model with given parameters"""
    # Extract parameters
    rate_values = {}
    for i, rate_id in enumerate(param_names):
        rate_values[rate_id] = params[i]

    E_field = params[-4]
    Te = params[-3]
    Ne = params[-2]
    H_drift_gain = params[-1]

    # Write modified chemistry file
    with open('chemistry.txt', 'r') as f:
        chem_lines = f.readlines()

    with open('chemistry_temp.txt', 'w') as f:
        for line in chem_lines:
            parts = line.split()
            if len(parts) > 0 and parts[0] in rate_values:
                # Replace this rate
                rate_id = parts[0]
                new_value = rate_values[rate_id]
                # Preserve the rest of the line (rate type, etc.)
                f.write(f"{rate_id} {new_value:.6e} " + " ".join(parts[2:]) + "\n")
            else:
                f.write(line)

    # Write input file
    input_data = {
        'chemistry_file': 'chemistry_temp.txt',
        'E_field': E_field,
        'Te': Te,
        'Ne': Ne,
        'H_drift_gain': H_drift_gain,
        **fixed_params
    }

    with open('input_temp.json', 'w') as f:
        json.dump(input_data, f, indent=2)

    # Run model
    import subprocess
    try:
        result = subprocess.run(
            ['python', '0D_plasma_model.py', 'input_temp.json'],
            capture_output=True,
            text=True,
            timeout=10
        )

        if result.returncode != 0:
            return None

        # Parse output
        output = json.loads(result.stdout)
        return output

    except Exception as e:
        return None

def objective(params):
    """Objective function: weighted error from targets"""
    global best_f, best_x, eval_count
    eval_count[0] += 1

    results = run_plasma_model(params)
    if results is None:
        return 1e10

    # NEW FOCUSED WEIGHTS - boost H and C2H2 to feed C2 production
    # Keep CH high to maintain suppression
    weights = {
        'H': 40.0,     # BOOSTED from 15 - need H for C2H2→C2 pathway
        'CH': 80.0,    # MAXIMUM - maintain CH suppression
        'C2': 60.0,    # High - primary target
        'C2H2': 20.0   # BOOSTED from 1 - need C2H2 to feed C2
    }

    species_error = 0.0
    for species in ['H', 'CH', 'C2', 'C2H2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        species_error += weights[species] * rel_error ** 2

    # Charge balance penalty - INCREASED to ensure good balance
    charge_penalty = 20.0 * results['charge_imbalance']**2

    # EXPLICIT pathway penalty for C2+H→CH flux
    # Continue suppressing this CH production source
    flux_penalty = 2.0 * (results['C2_H_CH_flux'] / 1e14) ** 2

    total_error = species_error + charge_penalty + flux_penalty

    # Track best solution
    if total_error < best_f:
        best_f = total_error
        best_x = params.copy()

        print(f"\n  *** NEW BEST: f(x) = {total_error:.2f} at evaluation {eval_count[0]}")
        print(f"      Te: {params[-3]:.2f} eV")
        print(f"      E-field: {params[-4]:.1f} V/cm")
        print(f"      Ne: {params[-2]:.2e} cm⁻³")
        print(f"      H_drift_gain: {params[-1]:.2e} cm⁻³/s")
        print(f"      Charge imbalance: {results['charge_imbalance']:.2f}%")

        for species in ['H', 'CH', 'C2', 'C2H2']:
            ratio = results[species] / TARGETS[species]
            print(f"      {species}: {results[species]:.2e} (target: {TARGETS[species]:.2e}, ratio: {ratio:.2f}x)")

        print(f"      C2+H→CH flux: {results['C2_H_CH_flux']:.2e} cm⁻³/s")

        # Break down error components
        species_err = sum(weights[sp] * ((results[sp] - TARGETS[sp]) / TARGETS[sp])**2
                         for sp in ['H', 'CH', 'C2', 'C2H2'])
        charge_err = 20.0 * results['charge_imbalance']**2
        flux_err = 2.0 * (results['C2_H_CH_flux'] / 1e14) ** 2
        print(f"      Species: {species_err:.2f}, Charge: {charge_err:.2f}, Flux: {flux_err:.2f}")

        # Save solution
        filename = f"optimization_results_C2H2_boost/best_f{total_error:.1f}_Hdrift{params[-1]:.2e}.json"
        save_data = {
            'Te': float(params[-3]),
            'Ne': float(params[-2]),
            'E_field': float(params[-4]),
            'H_drift_gain': float(params[-1]),
            'rate_values': {rate_id: float(params[i]) for i, rate_id in enumerate(param_names)},
            'target_densities': {sp: float(results[sp]) for sp in ['H', 'CH', 'C2']},
            'charge_balance': {
                'n_i_total': float(results.get('n_i_total', 0)),
                'n_e': float(results['Ne']),
                'imbalance_percent': float(results['charge_imbalance'])
            },
            'C2_H_CH_flux': float(results['C2_H_CH_flux']),
            'all_densities': {k: float(v) for k, v in results.items() if isinstance(v, (int, float))},
            'chemistry_analysis': results.get('chemistry_analysis', {})
        }

        with open(filename, 'w') as f:
            json.dump(save_data, f, indent=2)
        print(f"      Saved to: {filename}")

    # Progress indicator
    if eval_count[0] % 10 == 0:
        print(f"  [{eval_count[0]} evaluations, {eval_count[0]*0.007:.1f} min elapsed]")

    return total_error

print("="*80)
print(" Running Optimization (20 iterations, pop=8)")
print("="*80)
print()
print("Target: H, CH, C2 all within 0.6-1.4× range")
print("Strategy: Boost C2H2 to feed C2 production, maintain CH suppression, boost H")
print()

# Run optimization
result = differential_evolution(
    objective,
    bounds,
    maxiter=20,
    popsize=8,
    seed=42,
    disp=True,
    workers=1,
    updating='deferred',
    atol=0.01,
    tol=0.01
)

print()
print("="*80)
print(" OPTIMIZATION COMPLETE")
print("="*80)
print()
print(f"Best objective value: {best_f:.2f}")
print(f"Total evaluations: {eval_count[0]}")
print()

# Run final solution to get full output
final_results = run_plasma_model(best_x)
if final_results:
    print("Final solution:")
    print(f"  Te: {best_x[-3]:.2f} eV")
    print(f"  E-field: {best_x[-4]:.1f} V/cm")
    print(f"  Ne: {best_x[-2]:.2e} cm⁻³")
    print(f"  H_drift_gain: {best_x[-1]:.2e} cm⁻³/s")
    print()
    print("Target densities vs actual:")
    for species in ['H', 'CH', 'C2', 'C2H2']:
        ratio = final_results[species] / TARGETS[species]
        status = '✓' if 0.6 <= ratio <= 1.4 else '✗'
        print(f"  {species:5s}: {final_results[species]:.2e} / {TARGETS[species]:.2e} = {ratio:.2f}× {status}")

    print()
    print(f"Charge imbalance: {final_results['charge_imbalance']:.2f}%")
    print(f"C2+H→CH flux: {final_results['C2_H_CH_flux']:.2e} cm⁻³/s")
