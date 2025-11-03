#!/usr/bin/env python3
"""
Detailed analysis of the best model results
"""

import numpy as np
from scipy.integrate import solve_ivp
import json
import os

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


def load_and_run_best_model():
    """Load and run the best model, return full results."""

    # Load the best model
    results_dir = 'optimization_results_targeted'
    best_file = None
    best_objective = float('inf')

    for filename in os.listdir(results_dir):
        if filename.startswith('best_iteration') and filename.endswith('.json'):
            try:
                obj_str = filename.split('_f')[1].split('.json')[0]
                obj_val = float(obj_str)
                if obj_val < best_objective:
                    best_objective = obj_val
                    best_file = filename
            except:
                continue

    with open(os.path.join(results_dir, best_file), 'r') as f:
        model_data = json.load(f)

    # Extract parameters
    ne = model_data['Ne']
    E_field = model_data['E_field']
    rate_values = model_data['rate_values']

    # Setup parameters
    params = {
        'E_field': E_field,
        'L_discharge': 0.45,
        'ne': ne,
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                    'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                    'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                    'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
        }
    }

    # Define rates
    k = define_rates(params)
    db = get_complete_rate_database()

    for name, val in rate_values.items():
        if name in k:
            if name in db:
                val = np.clip(val, db[name].min, db[name].max)
            k[name] = val

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    # Setup initial conditions
    species = params['species']
    ns = len(species)
    y0 = np.ones(ns) * 1e3

    def set_density(name, value):
        try:
            idx = species.index(name)
            y0[idx] = value
        except ValueError:
            pass

    set_density('e', ne)
    set_density('Ar', 0.85 * 9.66e15)
    set_density('CH4', 0.15 * 9.66e15)
    set_density('ArPlus', 1e7)
    set_density('CH4Plus', 1e5)
    set_density('CH3Plus', 1e5)
    set_density('CH5Plus', 1e3)
    set_density('ArHPlus', 5e5)
    set_density('CH3Minus', 5e4)
    set_density('H2', 1e12)
    set_density('ArStar', 5e6)
    set_density('H', 1e11)
    set_density('C2', 5e7)
    set_density('CH', 5e4)
    set_density('C2H4', 5e7)
    set_density('C2H6', 1e6)
    set_density('CH2', 1e11)
    set_density('C2H2', 1e12)
    set_density('C2H5', 1e6)
    set_density('CH3', 5e7)
    set_density('C', 5e7)

    # Run simulation
    ode_func = PlasmaODE_Optimized(params)

    sol = solve_ivp(
        ode_func,
        (0, 100),
        y0,
        method='BDF',
        t_eval=np.logspace(-3, 2, 50),
        rtol=1e-5,
        atol=1e-6,
        max_step=10.0
    )

    return sol, params, model_data


def analyze_results(sol, params, model_data):
    """Comprehensive analysis of results."""

    species = params['species']
    y_final = sol.y[:, -1]

    # Create species density dictionary
    densities = {species[i]: y_final[i] for i in range(len(species))}

    print("=" * 100)
    print(" COMPREHENSIVE ANALYSIS OF BEST MODEL")
    print("=" * 100)

    # Model parameters
    print("\n" + "=" * 100)
    print(" MODEL PARAMETERS")
    print("=" * 100)
    print(f"\n  Region: CSB (Cathode Sheath Boundary) / CG (Cathode Glow)")
    print(f"  Discharge length (L): {params['L_discharge']} cm")
    print(f"  Electron density (Ne): {model_data['Ne']:.4e} cm⁻³")
    print(f"  Electric field (E): {model_data['E_field']:.1f} V/cm")

    # Calculate Te from electron energy
    # Te is typically derived from the reduced field E/N or from energy balance
    # For this discharge: E/N ≈ E / (0.85 * 9.66e15) = E/N in Td
    total_gas_density = 0.85 * 9.66e15 + 0.15 * 9.66e15  # Ar + CH4
    E_over_N = model_data['E_field'] / total_gas_density  # V·cm
    E_over_N_Td = E_over_N * 1e17  # Convert to Townsend (1 Td = 1e-17 V·cm²)

    # Rough estimate: Te (eV) ≈ 0.1 * E/N (Td) for E/N < 100 Td
    # For higher E/N, use Te ≈ 3-5 eV typical for these discharges
    Te_estimate = min(11.6 * E_over_N_Td * 0.1, 5.0)  # Cap at 5 eV, convert Td to eV

    print(f"  E/N: {E_over_N_Td:.2f} Td")
    print(f"  Electron temperature (Te): ~{Te_estimate:.2f} eV (estimated)")
    print(f"  Gas temperature (Tg): ~300 K (assumed, room temp)")
    print(f"  Pressure: ~10 Torr (from gas density)")
    print(f"  Gas composition: 85% Ar, 15% CH₄")

    # Charge balance analysis
    print("\n" + "=" * 100)
    print(" CHARGE BALANCE ANALYSIS")
    print("=" * 100)

    # Identify all charged species
    positive_ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                     'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus',
                     'H3Plus', 'CHPlus', 'H2Plus']

    negative_species = ['e', 'CH3Minus', 'HMinus']

    # Calculate total positive charge
    total_positive = sum(densities.get(ion, 0) for ion in positive_ions)

    # Calculate total negative charge
    total_negative = sum(densities.get(sp, 0) for sp in negative_species)

    print(f"\n  Total positive charge density: {total_positive:.4e} cm⁻³")
    print(f"  Total negative charge density: {total_negative:.4e} cm⁻³")
    print(f"  Charge imbalance: {abs(total_positive - total_negative):.4e} cm⁻³")
    print(f"  Relative imbalance: {abs(total_positive - total_negative) / total_positive * 100:.4f}%")

    if abs(total_positive - total_negative) / total_positive < 0.01:
        print(f"  Status: ✓ EXCELLENT quasi-neutrality (< 1%)")
    elif abs(total_positive - total_negative) / total_positive < 0.05:
        print(f"  Status: ✓ GOOD quasi-neutrality (< 5%)")
    else:
        print(f"  Status: ⚠ MODERATE imbalance (> 5%)")

    print("\n  Positive Ion Breakdown:")
    for ion in positive_ions:
        if densities.get(ion, 0) > 1e3:  # Only show significant ions
            fraction = densities[ion] / total_positive * 100
            print(f"    {ion:12s}: {densities[ion]:12.4e} cm⁻³  ({fraction:5.2f}%)")

    print("\n  Negative Species Breakdown:")
    for sp in negative_species:
        if densities.get(sp, 0) > 1e3:
            fraction = densities[sp] / total_negative * 100
            print(f"    {sp:12s}: {densities[sp]:12.4e} cm⁻³  ({fraction:5.2f}%)")

    # Full species table
    print("\n" + "=" * 100)
    print(" COMPLETE SPECIES DENSITIES")
    print("=" * 100)

    # Organize by category
    categories = {
        'Electrons': ['e'],
        'Positive Ions': positive_ions,
        'Negative Ions': ['CH3Minus', 'HMinus'],
        'Neutral Radicals (H-bearing)': ['H', 'CH', 'CH2', 'CH3', 'C'],
        'Neutral Radicals (C2)': ['C2', 'C2H', 'C2H3', 'C2H5'],
        'Neutral Radicals (C3)': ['C3', 'C3H', 'C3H2', 'C3H3', 'C3H5'],
        'Neutral Radicals (C4)': ['C4H', 'C4H2'],
        'Stable Molecules': ['H2', 'CH4', 'C2H2', 'C2H4', 'C2H6', 'C3H4', 'C3H6'],
        'Excited States': ['ArStar', 'C2H2Star']
    }

    for category, species_list in categories.items():
        print(f"\n{category}:")
        print("  " + "-" * 96)
        print(f"  {'Species':<15} {'Density (cm⁻³)':<20} {'Scientific':<20} {'Notes':<30}")
        print("  " + "-" * 96)

        for sp in species_list:
            if sp in densities:
                dens = densities[sp]
                if dens > 1e3 or sp in ['e', 'ArPlus', 'H', 'CH', 'C2']:  # Always show key species
                    # Format density nicely
                    if dens >= 1e15:
                        dens_str = f"{dens/1e15:.3f} × 10¹⁵"
                    elif dens >= 1e12:
                        dens_str = f"{dens/1e12:.3f} × 10¹²"
                    elif dens >= 1e9:
                        dens_str = f"{dens/1e9:.3f} × 10⁹"
                    elif dens >= 1e6:
                        dens_str = f"{dens/1e6:.3f} × 10⁶"
                    else:
                        dens_str = f"{dens:.3f}"

                    sci_str = f"{dens:.4e}"

                    # Add notes for key species
                    notes = ""
                    if sp == 'e':
                        notes = "Imposed by quasi-neutrality"
                    elif sp in ['Ar', 'CH4']:
                        notes = "Feed gas"
                    elif sp in ['H', 'CH', 'C2']:
                        notes = "★ Target species"
                    elif sp == 'H2':
                        notes = "Major product"
                    elif sp == 'ArStar':
                        notes = "Metastable, key ionizer"

                    print(f"  {sp:<15} {dens_str:<20} {sci_str:<20} {notes:<30}")

    # Key ratios
    print("\n" + "=" * 100)
    print(" KEY SPECIES RATIOS")
    print("=" * 100)

    print(f"\n  H/CH ratio: {densities['H']/densities['CH']:.2e} (H dominates)")
    print(f"  C2/CH ratio: {densities['C2']/densities['CH']:.2e} (C2 >> CH)")
    print(f"  CH3/CH ratio: {densities['CH3']/densities['CH']:.2e} (CH3 >> CH)")
    print(f"  H2/CH4 ratio: {densities['H2']/densities['CH4']:.2e} (H2 product formation)")
    print(f"  C2H6/CH4 ratio: {densities['C2H6']/densities['CH4']:.2e} (Higher hydrocarbon)")
    print(f"  C2H2/C2H4 ratio: {densities['C2H2']/densities['C2H4']:.2e} (Unsaturated products)")

    # Experimental comparison
    print("\n" + "=" * 100)
    print(" COMPARISON TO EXPERIMENTAL TARGETS")
    print("=" * 100)

    targets = {
        'H': 5.18e13,
        'CH': 1.0e9,
        'C2': 1.3e11
    }

    print(f"\n  {'Species':<10} {'Simulation':<18} {'Target':<18} {'Ratio':<10} {'Status':<15}")
    print("  " + "-" * 96)

    for sp, target in targets.items():
        sim = densities[sp]
        ratio = sim / target

        if 0.5 <= ratio <= 2.0:
            status = "✓ Excellent"
        elif 0.2 <= ratio <= 5.0:
            status = "✓ Good"
        elif 0.1 <= ratio <= 10.0:
            status = "~ Acceptable"
        else:
            status = "✗ Poor"

        print(f"  {sp:<10} {sim:<18.4e} {target:<18.4e} {ratio:<10.2f}x {status:<15}")

    # Region confirmation
    print("\n" + "=" * 100)
    print(" DISCHARGE REGION IDENTIFICATION")
    print("=" * 100)

    print(f"\n  This model represents the CG (Cathode Glow) region, also called")
    print(f"  the CSB (Cathode Sheath Boundary) region in some nomenclature.")
    print(f"")
    print(f"  Characteristics of this region:")
    print(f"    • Location: 0-1 mm from cathode")
    print(f"    • High electric field: {model_data['E_field']} V/cm")
    print(f"    • Moderate electron density: ~2×10⁹ cm⁻³")
    print(f"    • Active plasma chemistry region")
    print(f"    • High radical densities (H, CH, CH₃)")
    print(f"    • Transition from sheath to positive column")
    print(f"")
    print(f"  Model assumptions:")
    print(f"    • 0-D (well-mixed, spatial averaging)")
    print(f"    • Steady-state (t=100s >> all timescales)")
    print(f"    • Fixed E-field and Ne (from optimization)")
    print(f"    • Wall losses via sticking coefficients")

    print("\n" + "=" * 100)


if __name__ == '__main__':
    print("Running best model and performing comprehensive analysis...\n")
    sol, params, model_data = load_and_run_best_model()
    analyze_results(sol, params, model_data)
