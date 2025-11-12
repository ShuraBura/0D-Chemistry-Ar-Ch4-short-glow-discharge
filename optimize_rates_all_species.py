"""
Optimize multiple rate constants to match ALL THREE target densities

Use scipy optimization to automatically adjust rate constants within
physically reasonable bounds to minimize error from H, CH, C2 targets

Targets (all from same experiment):
- H:  2.52×10¹⁴ cm⁻³
- CH: 1.00×10⁹ cm⁻³
- C2: 5.60×10¹¹ cm⁻³
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, differential_evolution

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

def pressure_to_density(pressure_mTorr, T_K=300):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline_data = json.load(f)

targets = {'H': 2.52e14, 'CH': 1.0e9, 'C2': 5.6e11}
baseline_Te = baseline_data['Te']
baseline_E_field = baseline_data['E_field']
baseline_rate_values = baseline_data['rate_values']

# Parameters to optimize with bounds [min, max]
# Format: (param_name, baseline_value, min_multiplier, max_multiplier, description)
OPTIMIZATION_PARAMS = [
    ('ne', 2.3e9, 0.01, 10.0, 'Electron density (cm⁻³)'),
    ('C2_H_CH_C_cm3_7_6', 1.0, 0.001, 1.0, 'C2 + H → CH + C rate multiplier'),
    ('dust_density', 0.0, 0.0, 1e9, 'Dust density (cm⁻³)'),
    ('dust_sticking', 0.5, 0.1, 1.0, 'Dust sticking coefficient'),
    ('CH_H_C_H2_cm3_7_3', 1.0, 0.1, 10.0, 'CH + H → C + H2 rate multiplier'),
    ('stick_CH_9_3', 1.0, 0.1, 10.0, 'CH wall sticking multiplier'),
]

def run_simulation_with_params(params_array):
    """Run simulation with given parameter values"""

    # Unpack parameters
    ne_value = params_array[0]
    c2h_mult = params_array[1]
    dust_dens = params_array[2]
    dust_stick = params_array[3]
    ch_h_mult = params_array[4]
    stick_ch_mult = params_array[5]

    P = 500.0
    n_total = pressure_to_density(P)

    params = {
        'P': P,
        'Te': baseline_Te,
        'ne': ne_value,
        'E_field': baseline_E_field,
        'L_discharge': 0.45,
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
        },
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus',
                    'ArHPlus', 'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4',
                    'C2H6', 'CH2', 'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3',
                    'C3H2', 'CHPlus', 'C3H', 'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3',
                    'C3H5', 'C4H', 'C3H6', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                    'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
    }

    # Enable dust if density > 0
    if dust_dens > 1e6:
        params['enable_dust_loss'] = True
        params['dust_density'] = dust_dens
        params['dust_radius'] = 50e-7
        params['dust_sticking'] = dust_stick

    k = define_rates(params)

    # Apply baseline rates
    for rate_name, rate_value in baseline_rate_values.items():
        if rate_name in k:
            k[rate_name] = rate_value

    # Apply standard multipliers
    rate_mults = {
        'e_CH4_CH3_HMinus_cm3_8_1': 10.0,
        'ArStar_CH4_CH3_H_cm3_3_1': 1.4,
        'e_CH4Plus_CH3_H_cm3_6_4': 1.5,
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
    }

    # Apply optimized multipliers
    rate_mults['C2_H_CH_C_cm3_7_6'] = c2h_mult
    rate_mults['CH_H_C_H2_cm3_7_3'] = ch_h_mult
    rate_mults['stick_CH_9_3'] = stick_ch_mult

    for rate_name, mult in rate_mults.items():
        if rate_name in k:
            k[rate_name] *= mult

    params['k'] = k
    params['R'], params['tags'] = build_reactions(params)

    species = params['species']
    y0 = np.ones(len(species)) * 1e3
    y0[species.index('Ar')] = n_total * 0.85
    y0[species.index('CH4')] = n_total * 0.15
    y0[species.index('e')] = ne_value

    try:
        ode_func = PlasmaODE_Optimized(params)
        sol = solve_ivp(ode_func, (0, 500), y0, method='BDF', rtol=1e-7, atol=1e-9, max_step=1.0)

        if not sol.success:
            return None

        y_final = sol.y[:, -1]

        H = y_final[species.index('H')]
        CH = y_final[species.index('CH')]
        C2 = y_final[species.index('C2')]

        # Calculate Ni/Ne ratio
        ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH2Plus',
                'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus', 'H3Plus', 'CHPlus', 'H2Plus']
        Ni_total = sum(y_final[species.index(ion)] for ion in ions if ion in species)
        Ne = y_final[species.index('e')]

        if Ni_total / Ne > 10:
            return None  # Runaway

        return {
            'H': H,
            'CH': CH,
            'C2': C2,
            'Ni_Ne': Ni_total / Ne,
        }

    except Exception:
        return None


def objective_function(params_array, verbose=False):
    """
    Objective function to minimize
    Returns RMS percentage error from all three targets
    """
    result = run_simulation_with_params(params_array)

    if result is None:
        return 1e10  # Penalty for failure

    # Calculate percentage errors
    h_err = abs(result['H'] / targets['H'] - 1.0) * 100
    ch_err = abs(result['CH'] / targets['CH'] - 1.0) * 100
    c2_err = abs(result['C2'] / targets['C2'] - 1.0) * 100

    # RMS error
    rms_err = np.sqrt((h_err**2 + ch_err**2 + c2_err**2) / 3)

    if verbose:
        print(f"  H: {result['H']:.2e} ({result['H']/targets['H']*100:.1f}%), "
              f"CH: {result['CH']:.2e} ({result['CH']/targets['CH']*100:.1f}%), "
              f"C2: {result['C2']:.2e} ({result['C2']/targets['C2']*100:.1f}%), "
              f"RMS: {rms_err:.1f}%")

    return rms_err


if __name__ == '__main__':
    print("="*80)
    print("OPTIMIZE RATE CONSTANTS TO MATCH ALL THREE TARGETS")
    print("="*80)
    print(f"\nTargets:")
    print(f"  H:  {targets['H']:.2e} cm⁻³")
    print(f"  CH: {targets['CH']:.2e} cm⁻³")
    print(f"  C2: {targets['C2']:.2e} cm⁻³")

    print(f"\n{'='*80}")
    print("OPTIMIZATION PARAMETERS:")
    print(f"{'='*80}")
    for i, (name, baseline, min_mult, max_mult, desc) in enumerate(OPTIMIZATION_PARAMS):
        if 'ne' in name or 'dust_density' in name:
            print(f"{i+1}. {desc}")
            print(f"   Range: [{baseline*min_mult:.2e}, {baseline*max_mult:.2e}]")
        else:
            print(f"{i+1}. {desc}")
            print(f"   Multiplier range: [{min_mult}, {max_mult}]")

    # Set up bounds
    bounds = []
    x0 = []
    for name, baseline, min_mult, max_mult, desc in OPTIMIZATION_PARAMS:
        bounds.append((baseline * min_mult, baseline * max_mult))
        x0.append(baseline)  # Start from baseline

    print(f"\n{'='*80}")
    print("BASELINE EVALUATION:")
    print(f"{'='*80}")
    print("Starting from current best parameters...")
    baseline_error = objective_function(np.array(x0), verbose=True)
    print(f"\nBaseline RMS error: {baseline_error:.1f}%")

    print(f"\n{'='*80}")
    print("RUNNING DIFFERENTIAL EVOLUTION OPTIMIZER...")
    print(f"{'='*80}")
    print("This may take several minutes...")
    print("Strategy: Global optimization with population-based search\n")

    # Use differential evolution for global optimization
    result = differential_evolution(
        objective_function,
        bounds,
        strategy='best1bin',
        maxiter=50,
        popsize=15,
        tol=0.01,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=42,
        callback=lambda xk, convergence: print(f"Iteration complete, best score: {convergence:.1f}%"),
        polish=True,
        workers=1,
        updating='deferred',
    )

    print(f"\n{'='*80}")
    print("OPTIMIZATION RESULTS:")
    print(f"{'='*80}")

    if result.success or result.fun < baseline_error:
        print(f"\n✓ Optimization {'converged' if result.success else 'improved'}!")
        print(f"  Initial RMS error: {baseline_error:.1f}%")
        print(f"  Final RMS error:   {result.fun:.1f}%")
        print(f"  Improvement:       {baseline_error - result.fun:.1f}%")

        print(f"\nOptimal parameters:")
        for i, (name, baseline, min_mult, max_mult, desc) in enumerate(OPTIMIZATION_PARAMS):
            if 'ne' in name or 'dust_density' in name:
                print(f"  {desc}: {result.x[i]:.3e} (baseline: {baseline:.3e})")
            else:
                print(f"  {desc}: {result.x[i]:.3f}× (baseline: {baseline:.3f}×)")

        # Evaluate final result
        print(f"\n{'='*80}")
        print("FINAL SIMULATION RESULTS:")
        print(f"{'='*80}")
        final_result = run_simulation_with_params(result.x)

        if final_result:
            print(f"\nDensities:")
            print(f"  H:  {final_result['H']:.3e} cm⁻³ ({final_result['H']/targets['H']*100:6.2f}% of target)")
            print(f"  CH: {final_result['CH']:.3e} cm⁻³ ({final_result['CH']/targets['CH']*100:6.2f}% of target)")
            print(f"  C2: {final_result['C2']:.3e} cm⁻³ ({final_result['C2']/targets['C2']*100:6.2f}% of target)")
            print(f"  Ni/Ne: {final_result['Ni_Ne']:.2f}")

            print(f"\n{'='*80}")
            print("ASSESSMENT:")
            print(f"{'='*80}")

            if result.fun < 20:
                print(f"\n✓ EXCELLENT! All species within {result.fun:.0f}% RMS of targets")
                print("  Model successfully matches experimental data!")
            elif result.fun < 50:
                print(f"\n✓ GOOD! All species within {result.fun:.0f}% RMS of targets")
                print("  Reasonable agreement with experimental data")
            elif result.fun < baseline_error:
                print(f"\n⚠️  IMPROVED! RMS error reduced from {baseline_error:.0f}% to {result.fun:.0f}%")
                print("  But still significant discrepancy remains")
            else:
                print(f"\n✗ No improvement over baseline")

            # Save results
            print(f"\nSaving optimal parameters to optimized_rates.json...")
            opt_params = {
                'ne': result.x[0],
                'rate_multipliers': {
                    'C2_H_CH_C_cm3_7_6': result.x[1],
                    'CH_H_C_H2_cm3_7_3': result.x[4],
                    'stick_CH_9_3': result.x[5],
                },
                'dust_params': {
                    'enable_dust_loss': result.x[2] > 1e6,
                    'dust_density': result.x[2],
                    'dust_radius': 50e-7,
                    'dust_sticking': result.x[3],
                },
                'results': {
                    'H': float(final_result['H']),
                    'CH': float(final_result['CH']),
                    'C2': float(final_result['C2']),
                    'H_pct': float(final_result['H']/targets['H']*100),
                    'CH_pct': float(final_result['CH']/targets['CH']*100),
                    'C2_pct': float(final_result['C2']/targets['C2']*100),
                    'Ni_Ne': float(final_result['Ni_Ne']),
                    'rms_error': float(result.fun),
                }
            }

            with open('optimized_rates.json', 'w') as f:
                json.dump(opt_params, f, indent=2)

            print("  Saved to optimized_rates.json")

    else:
        print(f"\n✗ Optimization failed")
        print(f"  Message: {result.message}")
