#!/usr/bin/env python3
"""
Ar/CH4 Plasma Simulation Main Script
0-D model for short glow discharge
Converted from MATLAB ArCH4_PlasmaSimulationMainScript.m
"""

import numpy as np
from scipy.integrate import solve_ivp
import time
import matplotlib.pyplot as plt

from define_rates import define_rates
from build_reactions import build_reactions
from odefun import PlasmaODE


def calculate_initial_densities(params):
    """Calculate initial species densities."""
    species = params['species']
    ns = len(species)
    y0 = np.ones(ns) * 1e3  # Default minimum density

    # Helper function to set density
    def set_density(name, value):
        try:
            idx = species.index(name)
            y0[idx] = value
        except ValueError:
            pass

    # Set specific initial densities
    set_density('e', params['ne'])
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
    set_density('H3Plus', 1e4)
    set_density('HMinus', 1e4)
    set_density('C2HPlus', 1e5)
    set_density('CHPlus', 1e4)
    set_density('C3H6', 1e3)
    set_density('C2H2Star', 1e6)

    return y0


def main():
    """Main simulation function."""
    print("=" * 60)
    print("Ar/CH4 Plasma Simulation - 0D Chemistry Model")
    print("=" * 60)

    # Simulation parameters
    params = {
        'P': 0.4,  # Pressure in Torr
        'Tg': 400,  # Gas temperature in K
        'n_tot': 9.66e15,  # Total density in cm^-3
        'ne': 1e10,  # Electron density in cm^-3
        'Te': 1,  # Electron temperature in eV
        'E_field': 50,  # V/cm
        'L_discharge': 0.45,  # cm
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                    'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                    'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                    'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'C2H2Star'],
        'ion_species': ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH3Minus', 'H3Plus', 'CHPlus',
                        'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus'],
        'mobilities': {
            'ArPlus': 3057.28,
            'CH4Plus': 6432,
            'CH3Plus': 4949.6,
            'CH5Plus': 4761.6,
            'ArHPlus': 2969.6,
            'CH2Plus': 4949.6,
            'C2H5Plus': 4949.6,
            'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6,
            'C2HPlus': 5000,
            'H3Plus': 5000,
            'CHPlus': 5000,
            'CH3Minus': 3000,
            'HMinus': 3000
        }
    }

    # Build reaction network
    print("\nBuilding reaction network...")
    params['k'] = define_rates(params)
    params['R'], params['tags'] = build_reactions(params)
    print(f"  Number of species: {len(params['species'])}")
    print(f"  Number of reactions: {len(params['R'])}")

    # Initial conditions
    print("\nSetting initial conditions...")
    y0 = calculate_initial_densities(params)

    # Time span
    t_span = (0, 100)  # 0 to 100 seconds
    t_eval = [0, 0.1, 1, 10, 100]  # Evaluation times

    # Create ODE function
    ode_func = PlasmaODE(params)

    # Solver options
    print("\nStarting simulation...")
    print(f"  Time span: {t_span[0]} to {t_span[1]} s")
    print(f"  Method: BDF (stiff solver)")

    start_time = time.time()

    # Solve ODE system
    sol = solve_ivp(
        ode_func,
        t_span,
        y0,
        method='BDF',  # Backward Differentiation Formula for stiff problems
        t_eval=t_eval,
        rtol=1e-7,
        atol=1e-8,
        dense_output=True
    )

    elapsed_time = time.time() - start_time

    print(f"\nSimulation completed in {elapsed_time:.2f} seconds")
    print(f"  Status: {'Success' if sol.success else 'Failed'}")
    print(f"  Message: {sol.message}")
    print(f"  Number of function evaluations: {sol.nfev}")

    # Display results
    print("\n" + "=" * 60)
    print("Steady-State Densities at Final Time (t = 100 s)")
    print("=" * 60)

    # Neutral species
    neutral_species = ['H', 'CH', 'C2', 'H2', 'C2H2', 'CH3', 'CH2', 'Ar', 'CH4',
                       'ArStar', 'C2H4', 'C2H6', 'C2H5', 'C', 'C2H3', 'C3H2',
                       'C3H', 'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6']

    print("\nNeutral Species Densities (cm^-3):")
    for species_name in neutral_species:
        try:
            idx = params['species'].index(species_name)
            print(f"  {species_name:8s}: {sol.y[idx, -1]:.4e}")
        except ValueError:
            pass

    # Ion species
    print("\nIon Species Densities (cm^-3):")
    for species_name in params['ion_species'] + ['e']:
        try:
            idx = params['species'].index(species_name)
            print(f"  {species_name:8s}: {sol.y[idx, -1]:.4e}")
        except ValueError:
            pass

    # Key species
    print("\n" + "=" * 60)
    print("Key Species")
    print("=" * 60)
    key_species = ['H', 'CH', 'C2']
    for species_name in key_species:
        try:
            idx = params['species'].index(species_name)
            print(f"  Final {species_name} density: {sol.y[idx, -1]:.4e} cm^-3")
        except ValueError:
            pass

    # Plot results (simple time evolution of key species)
    print("\nGenerating plots...")
    plot_results(sol, params)

    print("\nSimulation complete!")


def plot_results(sol, params):
    """Plot time evolution of key species."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Key neutral radicals
    key_neutrals = ['H', 'CH', 'C2', 'CH3']
    ax = axes[0, 0]
    for species_name in key_neutrals:
        try:
            idx = params['species'].index(species_name)
            ax.semilogx(sol.t, sol.y[idx, :], 'o-', label=species_name, markersize=4)
        except ValueError:
            pass
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Density (cm$^{-3}$)')
    ax.set_title('Key Neutral Radicals')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Major ions
    key_ions = ['ArPlus', 'CH3Plus', 'CH5Plus']
    ax = axes[0, 1]
    for species_name in key_ions:
        try:
            idx = params['species'].index(species_name)
            ax.semilogx(sol.t, sol.y[idx, :], 'o-', label=species_name, markersize=4)
        except ValueError:
            pass
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Density (cm$^{-3}$)')
    ax.set_title('Major Ions')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Stable molecules
    key_molecules = ['H2', 'CH4', 'C2H2', 'C2H4', 'C2H6']
    ax = axes[1, 0]
    for species_name in key_molecules:
        try:
            idx = params['species'].index(species_name)
            ax.semilogx(sol.t, sol.y[idx, :], 'o-', label=species_name, markersize=4)
        except ValueError:
            pass
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Density (cm$^{-3}$)')
    ax.set_title('Stable Molecules')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Excited/metastable
    key_excited = ['ArStar', 'C2H2Star']
    ax = axes[1, 1]
    for species_name in key_excited:
        try:
            idx = params['species'].index(species_name)
            ax.semilogx(sol.t, sol.y[idx, :], 'o-', label=species_name, markersize=4)
        except ValueError:
            pass
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Density (cm$^{-3}$)')
    ax.set_title('Excited/Metastable Species')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('plasma_simulation_results.png', dpi=150)
    print("  Plot saved to: plasma_simulation_results.png")


if __name__ == '__main__':
    main()
