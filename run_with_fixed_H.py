#!/usr/bin/env python3
"""
Run 0-D chemistry simulation with FIXED H density from experimental profile.

This script demonstrates how to:
1. Load experimental nH(x) profile
2. Calculate spatial average
3. Fix H density (non-evolving) during simulation
4. Compare other species (CH, C2) against targets
"""

import numpy as np
from scipy.integrate import solve_ivp
import time

# Import simulation modules
from species_list import get_species_list
from define_rates_tunable import define_rates_optimized


def load_experimental_H_profile():
    """
    Load your experimental nH(x) profile.

    Replace this with actual data loading from your measurements.

    Returns:
        nH_avg: Spatially-averaged H density (cm⁻³)
    """
    # OPTION A: If you have a data file
    # x, nH = np.loadtxt('your_nH_profile.csv', delimiter=',', unpack=True)
    # nH_avg = np.trapz(nH, x) / (x[-1] - x[0])

    # OPTION B: Use your measured value directly
    # From EXPERIMENTAL_TARGETS.md:
    nH_avg = 8.58e15  # cm⁻³ (CG region, 0-1 mm average)

    return nH_avg


class ODEFunc_FixedH:
    """
    Modified ODE function with FIXED H density.

    H density is held constant at the experimental value.
    """
    def __init__(self, species, k, params):
        self.species = species
        self.k = k
        self.params = params
        self.dydt = np.zeros(len(species))

        # Store indices
        self.e_idx = species.index('e')
        self.Ar_idx = species.index('Ar')
        self.CH4_idx = species.index('CH4')
        self.H_idx = species.index('H')

        # Fixed H density (from experiment)
        self.H_fixed = params.get('H_fixed', 8.58e15)  # cm⁻³

    def __call__(self, t, y):
        """Compute dy/dt with fixed H density."""
        self.dydt.fill(0.0)

        # Enforce fixed H density
        y[self.H_idx] = self.H_fixed

        # Get densities
        species_dict = {sp: y[i] for i, sp in enumerate(self.species)}
        ne = species_dict['e']

        # Calculate reaction rates (using fixed H)
        for reaction_id, rate in self.k.items():
            # Parse reaction: e.g., 'e_CH4_CH3_H_eV_1_1'
            parts = reaction_id.split('_')

            # Simple rate calculation (you'd use your actual rate logic)
            # This is placeholder - use your real reaction processing
            if 'e' in parts[:2]:  # Electron-impact
                R = rate * ne * species_dict.get(parts[1], 0)
            else:  # Neutral-neutral
                R = rate * species_dict.get(parts[0], 0) * species_dict.get(parts[1], 0)

            # Apply to species (simplified - use your actual stoichiometry)
            # This would need your full reaction parsing logic
            pass

        # Fix constant species (including H)
        self.dydt[self.e_idx] = 0
        self.dydt[self.Ar_idx] = 0
        self.dydt[self.CH4_idx] = 0
        self.dydt[self.H_idx] = 0  # KEY: H is fixed!

        return self.dydt


def run_simulation_fixed_H(params, H_fixed):
    """
    Run simulation with fixed H density.

    Args:
        params: Dictionary of plasma parameters
        H_fixed: Fixed H density from experiment (cm⁻³)

    Returns:
        result: Integration result with final densities
    """
    # Add fixed H to params
    params['H_fixed'] = H_fixed

    # Get species list
    species_list = get_species_list()

    # Initialize densities
    y0 = np.zeros(len(species_list))

    def set_density(name, value):
        y0[species_list.index(name)] = value

    # Set initial conditions
    set_density('e', params['ne'])
    set_density('Ar', 0.85 * 9.66e15)
    set_density('CH4', 0.15 * 9.66e15)
    set_density('H', H_fixed)  # Use experimental value

    # Set other species initial guesses
    set_density('ArPlus', params['ne'] * 0.5)
    set_density('CH4Plus', params['ne'] * 0.05)
    set_density('CH3Plus', params['ne'] * 0.05)
    set_density('CH5Plus', params['ne'] * 0.01)
    set_density('ArHPlus', params['ne'] * 0.1)
    set_density('CH3Minus', params['ne'] * 0.01)
    set_density('H2', 5e12)
    set_density('ArStar', params['ne'] * 0.5)
    set_density('C2', 1e11)
    set_density('CH', 1e8)
    set_density('C2H4', 5e7)
    set_density('C2H6', 1e6)
    set_density('CH2', 1e11)
    set_density('C2H2', 1e12)
    set_density('C2H5', 1e6)
    set_density('CH3', 5e7)
    set_density('C', 5e7)

    # Define rate constants
    k = define_rates_optimized(params)

    # Create ODE function (would need full implementation)
    # For now, use your existing odefun_optimized with modification
    from odefun_optimized import ODEFunc
    ode_func = ODEFunc(species_list, k, params)

    # MODIFY: Add fixed H constraint
    original_call = ode_func.__call__

    def fixed_H_wrapper(t, y):
        # Enforce fixed H
        y[ode_func.H_idx] = H_fixed
        dydt = original_call(t, y)
        dydt[ode_func.H_idx] = 0  # Don't let H evolve
        return dydt

    ode_func.__call__ = fixed_H_wrapper

    # Time span
    t_span = (0, 5e-3)  # 5 ms
    t_eval = np.logspace(-7, np.log10(5e-3), 100)

    print("Starting integration with fixed H...")
    t0 = time.time()

    result = solve_ivp(
        ode_func,
        t_span,
        y0,
        method='BDF',
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-10
    )

    t1 = time.time()
    print(f"Integration completed in {t1-t0:.2f} seconds")

    return result, species_list


def analyze_results(result, species_list, H_fixed):
    """
    Analyze results with fixed H.

    Args:
        result: Integration result
        species_list: List of species names
        H_fixed: Fixed H density (cm⁻³)
    """
    # Get final densities
    y_final = result.y[:, -1]

    def get_density(name):
        return y_final[species_list.index(name)]

    # Target densities (CG region)
    target_H = 8.58e15
    target_CH = 4.6e8
    target_C2 = 1.44e11

    # Predicted densities
    pred_H = H_fixed  # Fixed by design
    pred_CH = get_density('CH')
    pred_C2 = get_density('C2')

    # Calculate errors
    error_H = abs(pred_H - target_H) / target_H * 100
    error_CH = abs(pred_CH - target_CH) / target_CH * 100
    error_C2 = abs(pred_C2 - target_C2) / target_C2 * 100

    print("\n" + "="*60)
    print("RESULTS WITH FIXED H DENSITY")
    print("="*60)
    print(f"\n{'Species':<10} {'Target':<12} {'Predicted':<12} {'Error':<10}")
    print("-"*60)
    print(f"{'H':<10} {target_H:<12.2e} {pred_H:<12.2e} {error_H:<10.1f}%")
    print(f"{'CH':<10} {target_CH:<12.2e} {pred_CH:<12.2e} {error_CH:<10.1f}%")
    print(f"{'C2':<10} {target_C2:<12.2e} {pred_C2:<12.2e} {error_C2:<10.1f}%")
    print("-"*60)

    # Ratios
    ratio_H_CH = pred_H / pred_CH
    ratio_C2_CH = pred_C2 / pred_CH
    target_ratio_H_CH = target_H / target_CH
    target_ratio_C2_CH = target_C2 / target_CH

    print(f"\n{'Ratio':<15} {'Target':<12} {'Predicted':<12} {'Match':<10}")
    print("-"*60)
    print(f"{'H/CH':<15} {target_ratio_H_CH:<12.2e} {ratio_H_CH:<12.2e} "
          f"{ratio_H_CH/target_ratio_H_CH:<10.2f}x")
    print(f"{'C2/CH':<15} {target_ratio_C2_CH:<12.2e} {ratio_C2_CH:<12.2e} "
          f"{ratio_C2_CH/target_ratio_C2_CH:<10.2f}x")
    print("="*60)

    return {
        'H': pred_H,
        'CH': pred_CH,
        'C2': pred_C2,
        'error_CH': error_CH,
        'error_C2': error_C2
    }


def main():
    """Main execution."""
    print("="*60)
    print("0-D CHEMISTRY WITH FIXED H DENSITY FROM EXPERIMENT")
    print("="*60)

    # Load experimental H profile
    H_fixed = load_experimental_H_profile()
    print(f"\nUsing fixed H density: {H_fixed:.2e} cm⁻³")
    print("(from spatial average of experimental nH(x) profile)")

    # Define plasma parameters
    params = {
        'P': 0.4,          # Torr
        'Tg': 570,         # K
        'ne': 5e8,         # cm⁻³
        'Te': 3.0,         # eV
        'E_field': 600,    # V/cm
        'L_diff': 0.1,     # cm (measured)
        'L_discharge': 0.45,  # cm
        'gamma_H': 0.01,   # Wall recombination coefficient
    }

    print("\nPlasma parameters:")
    for key, val in params.items():
        print(f"  {key:<15} = {val}")

    # Run simulation
    result, species_list = run_simulation_fixed_H(params, H_fixed)

    # Analyze results
    results_dict = analyze_results(result, species_list, H_fixed)

    # Interpretation
    print("\nINTERPRETATION:")
    print("-" * 60)
    if results_dict['error_CH'] < 50 and results_dict['error_C2'] < 50:
        print("✓ GOOD MATCH for CH and C2 with fixed H!")
        print("  → Chemistry model validates well")
        print("  → H density constraint is reasonable")
    elif results_dict['error_CH'] < 200 and results_dict['error_C2'] < 200:
        print("~ MODERATE MATCH for CH and C2")
        print("  → May need parameter tuning (Te, ne, E_field)")
        print("  → Consider running parameter sweep")
    else:
        print("✗ POOR MATCH for CH and/or C2")
        print("  → Fixed H may be too constraining")
        print("  → Consider Option 1: Let H evolve dynamically")
        print("  → Or adjust H drift/loss mechanisms")

    print("\nNEXT STEPS:")
    print("1. If match is poor, try varying Te, ne, E_field")
    print("2. Run parameter_sweep_cg.py to find optimal parameters")
    print("3. Compare fixed-H vs dynamic-H approaches")
    print("="*60)


if __name__ == '__main__':
    main()
