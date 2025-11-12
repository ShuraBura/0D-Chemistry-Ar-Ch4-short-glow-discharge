"""
Optimize dust parameters to match ALL experimental targets simultaneously

Goal: Find dust parameters where:
- CH ≈ 100% of target (not 6810%!)
- C2 ≈ 100% of target (maintain our success!)
- H ≈ 100% of target (not too low!)

Strategy: Sweep dust density and radius to find optimal balance
"""

import numpy as np
import json
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

from define_rates import define_rates
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized

def pressure_to_density(pressure_mTorr, T_K=300):
    pressure_Pa = pressure_mTorr * 0.133322
    k_B = 1.380649e-23
    n_m3 = pressure_Pa / (k_B * T_K)
    return n_m3 * 1e-6

with open('optimization_results_comprehensive_1e12/best_f70.3.json') as f:
    baseline = json.load(f)

targets = {'H': 2.52e14, 'CH': 1.0e9, 'C2': 5.6e11}

def run_simulation(dust_params, ne_value=2.3e9, verbose=False):
    """Run simulation with given dust parameters"""

    P = 500.0
    n_total = pressure_to_density(P)

    params = {
        'P': P,
        'Te': baseline['Te'],
        'ne': ne_value,
        'E_field': baseline['E_field'],
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
        'enable_dust_loss': True,
    }

    params.update(dust_params)

    k = define_rates(params)

    # Apply baseline rates
    for rate_name, rate_value in baseline['rate_values'].items():
        if rate_name in k:
            k[rate_name] = rate_value

    # Apply physically realistic multipliers
    rate_mults = {
        'e_CH4_CH3_HMinus_cm3_8_1': 10.0,
        'ArStar_CH4_CH3_H_cm3_3_1': 1.4,
        'e_CH4Plus_CH3_H_cm3_6_4': 1.5,
        'stick_CH3_9_2': 0.01,
        'stick_C2H2_9_11': 0.01,
        'loss_C2H2_11_19': 0.01,
    }

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
            'H_pct': H / targets['H'] * 100,
            'CH_pct': CH / targets['CH'] * 100,
            'C2_pct': C2 / targets['C2'] * 100,
        }

    except Exception as e:
        if verbose:
            print(f"Error: {e}")
        return None


def objective_function(dust_density_log, dust_radius, dust_sticking, weight_CH=1.0, weight_C2=1.0, weight_H=0.5):
    """
    Objective function to minimize
    Goal: All species near 100%
    """
    dust_density = 10**dust_density_log

    dust_params = {
        'dust_density': dust_density,
        'dust_radius': dust_radius * 1e-7,  # Convert nm to cm
        'dust_sticking': dust_sticking,
    }

    result = run_simulation(dust_params)

    if result is None:
        return 1e10  # Penalty for failure

    # Calculate weighted error from 100% target
    error_CH = abs(result['CH_pct'] - 100.0) * weight_CH
    error_C2 = abs(result['C2_pct'] - 100.0) * weight_C2
    error_H = abs(result['H_pct'] - 100.0) * weight_H

    total_error = error_CH + error_C2 + error_H

    return total_error


if __name__ == '__main__':
    print("="*80)
    print("OPTIMIZE DUST PARAMETERS FOR MULTI-SPECIES TARGETS")
    print("="*80)

    print("\n" + "="*80)
    print("PARAMETER SWEEP: Find optimal dust configuration")
    print("="*80 + "\n")

    # Sweep dust density at fixed radius and sticking
    print("Strategy: Sweep dust density (n_dust) with r=50nm, α=0.5")
    print("-"*80)

    # Logarithmic sweep from 1e7 to 1e10 cm⁻³
    dust_densities = np.logspace(7, 10, 20)
    dust_radius = 50e-7  # cm (50 nm)
    dust_sticking = 0.5

    results = []

    print(f"{'n_dust (cm⁻³)':<15} {'CH %':<10} {'C2 %':<10} {'H %':<10} {'Score':<12} {'Status'}")
    print("-"*80)

    best_score = 1e10
    best_params = None
    best_result = None

    for n_dust in dust_densities:
        dust_params = {
            'dust_density': n_dust,
            'dust_radius': dust_radius,
            'dust_sticking': dust_sticking,
        }

        result = run_simulation(dust_params)

        if result is None:
            print(f"{n_dust:<15.2e} {'FAILED':<10} {'FAILED':<10} {'FAILED':<10} {'---':<12} ✗")
            continue

        # Score: sum of absolute deviations from 100%
        score = abs(result['CH_pct'] - 100) + abs(result['C2_pct'] - 100) + 0.5 * abs(result['H_pct'] - 100)

        status = "✓" if result['Ni_Ne'] < 10 else "✗"

        print(f"{n_dust:<15.2e} {result['CH_pct']:<10.1f} {result['C2_pct']:<10.1f} {result['H_pct']:<10.1f} {score:<12.1f} {status}")

        results.append((n_dust, result, score))

        if score < best_score:
            best_score = score
            best_params = dust_params
            best_result = result

    print("\n" + "="*80)
    print("BEST CONFIGURATION FROM SWEEP:")
    print("="*80)

    if best_params:
        print(f"\nDust parameters:")
        print(f"  n_dust = {best_params['dust_density']:.2e} cm⁻³")
        print(f"  r_dust = {best_params['dust_radius']*1e7:.0f} nm")
        print(f"  α_dust = {best_params['dust_sticking']:.2f}")

        print(f"\nResults:")
        print(f"  H:  {best_result['H_pct']:6.1f}% of target")
        print(f"  CH: {best_result['CH_pct']:6.1f}% of target")
        print(f"  C2: {best_result['C2_pct']:6.1f}% of target")
        print(f"  Score: {best_score:.1f}")
        print(f"  Ni/Ne: {best_result['Ni_Ne']:.2f}")

        # Calculate improvement
        baseline_score = abs(6810 - 100) + abs(107 - 100) + 0.5 * abs(104 - 100)
        improvement = (baseline_score - best_score) / baseline_score * 100
        print(f"\nImprovement over baseline: {improvement:.1f}%")

    print("\n" + "="*80)
    print("FINE-TUNE AROUND BEST PARAMETERS:")
    print("="*80 + "\n")

    # Fine-tune around best density
    if best_params:
        best_n = best_params['dust_density']
        print(f"Fine-tuning around n_dust = {best_n:.2e} cm⁻³")
        print("-"*80)

        # Sweep ±30% around best
        fine_densities = np.linspace(best_n * 0.7, best_n * 1.3, 15)

        print(f"{'n_dust (cm⁻³)':<15} {'CH %':<10} {'C2 %':<10} {'H %':<10} {'Score':<12}")
        print("-"*80)

        for n_dust in fine_densities:
            dust_params = {
                'dust_density': n_dust,
                'dust_radius': dust_radius,
                'dust_sticking': dust_sticking,
            }

            result = run_simulation(dust_params)

            if result is None:
                continue

            score = abs(result['CH_pct'] - 100) + abs(result['C2_pct'] - 100) + 0.5 * abs(result['H_pct'] - 100)

            marker = " ← BEST" if score < best_score else ""

            print(f"{n_dust:<15.2e} {result['CH_pct']:<10.1f} {result['C2_pct']:<10.1f} {result['H_pct']:<10.1f} {score:<12.1f}{marker}")

            if score < best_score:
                best_score = score
                best_params = dust_params
                best_result = result

    print("\n" + "="*80)
    print("OPTIMIZE RADIUS AND STICKING COEFFICIENT:")
    print("="*80 + "\n")

    # Try different radii
    radii = [30, 40, 50, 60, 70]  # nm
    stickings = [0.3, 0.5, 0.7]

    print(f"Testing combinations at n_dust = {best_params['dust_density']:.2e} cm⁻³")
    print("-"*80)
    print(f"{'r (nm)':<8} {'α':<6} {'CH %':<10} {'C2 %':<10} {'H %':<10} {'Score':<12}")
    print("-"*80)

    for r in radii:
        for alpha in stickings:
            dust_params = {
                'dust_density': best_params['dust_density'],
                'dust_radius': r * 1e-7,
                'dust_sticking': alpha,
            }

            result = run_simulation(dust_params)

            if result is None:
                continue

            score = abs(result['CH_pct'] - 100) + abs(result['C2_pct'] - 100) + 0.5 * abs(result['H_pct'] - 100)

            marker = " ← BEST" if score < best_score else ""

            print(f"{r:<8} {alpha:<6.1f} {result['CH_pct']:<10.1f} {result['C2_pct']:<10.1f} {result['H_pct']:<10.1f} {score:<12.1f}{marker}")

            if score < best_score:
                best_score = score
                best_params = dust_params
                best_result = result

    print("\n" + "="*80)
    print("FINAL OPTIMIZED CONFIGURATION:")
    print("="*80)

    print(f"\nOptimal dust parameters:")
    print(f"  n_dust = {best_params['dust_density']:.3e} cm⁻³")
    print(f"  r_dust = {best_params['dust_radius']*1e7:.1f} nm")
    print(f"  α_dust = {best_params['dust_sticking']:.2f}")

    print(f"\nFinal results:")
    print(f"  H:  {best_result['H']:.2e} cm⁻³ ({best_result['H_pct']:6.2f}% of target)")
    print(f"  CH: {best_result['CH']:.2e} cm⁻³ ({best_result['CH_pct']:6.2f}% of target)")
    print(f"  C2: {best_result['C2']:.2e} cm⁻³ ({best_result['C2_pct']:6.2f}% of target)")
    print(f"  Ni/Ne: {best_result['Ni_Ne']:.2f}")

    print(f"\nOptimization score: {best_score:.1f}")
    print(f"  (Lower is better; perfect score = 0)")

    # Compare to baseline
    print(f"\nComparison to baseline (no dust):")
    print(f"  Baseline CH: 6810% → Optimized CH: {best_result['CH_pct']:.1f}%")
    print(f"  Baseline C2:  107% → Optimized C2: {best_result['C2_pct']:.1f}%")
    print(f"  Baseline H:   104% → Optimized H:  {best_result['H_pct']:.1f}%")

    print(f"\n{'='*80}")
    print("CONCLUSION:")
    print(f"{'='*80}\n")

    if best_result['CH_pct'] < 200 and best_result['C2_pct'] > 50 and best_result['H_pct'] > 20:
        print("✓ SUCCESS! Found dust parameters that balance all three species!")
        print(f"  All targets within reasonable range of 100%")
    elif best_result['CH_pct'] < 500:
        print("⚠️  PARTIAL SUCCESS: CH significantly reduced but not all targets met")
    else:
        print("✗ Could not find parameters that satisfy all constraints")
        print("   Trade-off between CH reduction and C2/H preservation remains")
