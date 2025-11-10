#!/usr/bin/env python3
"""
COMPREHENSIVE OPTIMIZATION WITH CHARGE BALANCE CONSTRAINT

Same as full optimization but with:
1. Charge balance penalty in objective function (Ni/Ne in [2, 7] range)
2. Ne constrained to 2.3e9 ± 30% = [1.61e9, 2.99e9] cm⁻³

This ensures physically realistic plasma with proper ionization balance.
User clarified: "The balance does not need to be perfect it should be Ni/Ne ~2-7"
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import time
import signal
from contextlib import contextmanager
import os
import json

from define_rates import define_rates
from rate_database_complete import get_complete_rate_database, get_tunable_rates_for_target
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


# Experimental targets (updated 2025-11-10 based on spatial profile)
TARGETS = {
    'H': 2.52e14,     # Radially averaged (r=-6 to 6mm) at y=4mm from experimental profile
    'CH': 1.0e9,      # Unchanged
    'C2': 5.6e11,     # Updated from 1.3e11
}

# High pressure setting
PRESSURE_MTORR = 500.0

# Ne constraint - force higher Ne for H production (full CH4 dissociation)
NE_TARGET = 2.5e9
NE_MIN = 2.0e9  # Force higher Ne to boost H production
NE_MAX = 3.0e9  # Allow up to 3e9

# Create results directory
os.makedirs('optimization_results_charge_balanced', exist_ok=True)

# Track best result
best_result = {'objective': 1e10, 'params': None, 'densities': None}


class TimeoutException(Exception):
    pass


@contextmanager
def time_limit(seconds):
    """Context manager to timeout function calls."""
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


def pressure_to_density(pressure_mTorr, T_K=400):
    """Convert pressure to number density."""
    kB = 1.38064852e-23  # J/K
    Torr_to_Pa = 133.322
    P_Pa = pressure_mTorr * 1e-3 * Torr_to_Pa
    n_m3 = P_Pa / (kB * T_K)
    n_cm3 = n_m3 * 1e-6
    return n_cm3


def select_tunable_rates():
    """Select ~40 most important tunable rates."""
    db = get_complete_rate_database()

    # Critical reactions from chemistry analysis
    critical_reactions = [
        'CH3_CH3_C2H2_H2_H2_cm3_7_49',
        'C_CH3_C2_H2_H_cm3_7_8',
        'CH3_C2H5_C2H2_CH3_H2_cm3_7_61',
        'CH2_CH3_C2H2_H_H2_cm3_7_62',
        'CH2_C2H3_C2H2_CH3_cm3_7_53',
        'CH2_CH2_C2H2_H2_cm3_7_15',
        'C2H2_H_C2_H2_H_cm3_7_50',
        'C2H2_C_C2_CH2_cm3_7_19',
        'C2_H_CH_C_cm3_7_6',
        'C_H_CH_cm3_7_12',
        'CH2_H_CH_H2_cm3_7_1',
        'CH_CH4_C2H4_H_cm3_7_20',
        'CH_CH4_CH2_CH3_cm3_7_39',
        'CH_H_C_H2_cm3_7_3',
        'CH_CH3_C2H4_cm3_7_5',
        'e_CH4_CH3_H_cm3_1_1',
        'CH3_H_CH4_cm3_7_35',
        'e_CH4_CH2_H2_cm3_1_2',
    ]

    h_rates = set(get_tunable_rates_for_target('H').keys())
    ch_rates = set(get_tunable_rates_for_target('CH').keys())
    c2_rates = set(get_tunable_rates_for_target('C2').keys())

    target_rates = h_rates | ch_rates | c2_rates

    large_range_rates = {
        name for name, rate in db.items()
        if (rate.max / rate.min > 3.0) and (rate.min > 0)
    }

    flagged_rates = {name for name, rate in db.items() if rate.flag}

    selected = set(critical_reactions) | target_rates | large_range_rates | flagged_rates
    selected = {name for name in selected if name in db}

    selected_with_range = [
        (name, db[name].max / db[name].min if db[name].min > 0 else 1.0)
        for name in selected
    ]
    selected_with_range.sort(key=lambda x: x[1], reverse=True)

    selected_names = [name for name, _ in selected_with_range[:40]]

    return selected_names


def run_simulation(rate_values, Te, ne, E_field, params_base, log_file=None):
    """Run simulation and calculate charge balance."""

    try:
        n_total = pressure_to_density(PRESSURE_MTORR)

        params = params_base.copy()
        params['E_field'] = E_field
        params['ne'] = ne
        params['Te'] = Te

        k = define_rates(params)
        db = get_complete_rate_database()

        for name, val in rate_values.items():
            if name in k and name in db:
                val = np.clip(val, db[name].min, db[name].max)
                k[name] = val

        for name, rate_db in db.items():
            if name in k:
                if k[name] < rate_db.min:
                    k[name] = rate_db.min
                elif k[name] > rate_db.max:
                    k[name] = rate_db.max

        params['k'] = k
        params['R'], params['tags'] = build_reactions(params)

        species = params['species']
        ns = len(species)
        y0 = np.ones(ns) * 1e3

        def set_density(name, value):
            try:
                idx = species.index(name)
                y0[idx] = value
            except ValueError:
                pass

        n_Ar = 0.85 * n_total
        n_CH4 = 0.15 * n_total

        set_density('e', ne)
        set_density('Ar', n_Ar)
        set_density('CH4', n_CH4)
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

        ode_func = PlasmaODE_Optimized(params)

        try:
            with time_limit(120):  # Increased from 30s to allow full convergence
                sol = solve_ivp(
                    ode_func,
                    (0, 100),
                    y0,
                    method='BDF',
                    rtol=1e-5,
                    atol=1e-6,
                    max_step=10.0
                )
        except TimeoutException:
            return None

        if not sol.success:
            return None

        # DEBUG: Print integration endpoint
        print(f"  [DEBUG] Integration: t_final={sol.t[-1]:.2f}s, nsteps={len(sol.t)}, H_final={sol.y[species.index('H'), -1]:.2e}")

        y_final = sol.y[:, -1]

        def get_density(name):
            try:
                idx = species.index(name)
                return y_final[idx]
            except ValueError:
                return 0.0

        # Calculate total positive ion density
        positive_ions = [
            'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
            'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus',
            'H3Plus', 'CHPlus', 'H2Plus'
        ]

        n_i_total = sum(max(0, get_density(ion)) for ion in positive_ions)

        results = {
            'H': get_density('H'),
            'CH': get_density('CH'),
            'C2': get_density('C2'),
            'C2H2': get_density('C2H2'),
            'n_i_total': n_i_total,
            'Ni_Ne_ratio': n_i_total / ne,
            'charge_imbalance': abs(n_i_total - ne) / ne,
        }

        if log_file:
            all_densities = {species[i]: float(y_final[i]) for i in range(ns)}
            log_data = {
                'pressure_mTorr': PRESSURE_MTORR,
                'n_total': float(n_total),
                'Te': Te,
                'Ne': ne,
                'E_field': E_field,
                'n_i_total': float(n_i_total),
                'Ni_over_Ne': float(n_i_total / ne),
                'charge_imbalance_pct': float(results['charge_imbalance'] * 100),
                'rate_values': {k: float(v) for k, v in rate_values.items()},
                'target_densities': {k: float(results[k]) for k in ['H', 'CH', 'C2']},
                'all_densities': all_densities,
            }

            with open(log_file, 'w') as f:
                json.dump(log_data, f, indent=2)

        return results

    except Exception as e:
        return None


def objective_function(x, param_names, params_base):
    """Objective function WITH CHARGE BALANCE PENALTY."""
    global best_result

    if not hasattr(objective_function, 'counter'):
        objective_function.counter = 0
        objective_function.start_time = time.time()

    objective_function.counter += 1

    if objective_function.counter % 20 == 0:
        elapsed = time.time() - objective_function.start_time
        print(f"  [{objective_function.counter} evaluations, {elapsed/60:.1f} min elapsed]")

    rate_values = {name: val for name, val in zip(param_names[:-3], x[:-3])}
    Te = x[-3]
    ne = x[-2]
    E_field = x[-1]

    results = run_simulation(rate_values, Te, ne, E_field, params_base)

    if results is None:
        return 1e10

    # Species target errors
    weights = {
        'H': 1.0,
        'CH': 20.0,
        'C2': 3.0
    }

    species_error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        species_error += weights[species] * rel_error ** 2

    # CHARGE BALANCE PENALTY - Target Ni/Ne in [2, 7] range
    # Only penalize if outside this range (user clarified: "Ni/Ne ~2-7")
    Ni_Ne_ratio = results['Ni_Ne_ratio']
    if Ni_Ne_ratio < 2.0:
        charge_penalty = 50.0 * (2.0 - Ni_Ne_ratio)**2
    elif Ni_Ne_ratio > 7.0:
        charge_penalty = 50.0 * (Ni_Ne_ratio - 7.0)**2
    else:
        charge_penalty = 0.0  # In acceptable range [2, 7]

    # Total error
    total_error = species_error + charge_penalty

    # Track best result
    if total_error < best_result['objective']:
        best_result['objective'] = total_error
        best_result['params'] = {
            'Te': Te,
            'Ne': ne,
            'E_field': E_field,
            'pressure_mTorr': PRESSURE_MTORR,
            'rates': dict(rate_values)
        }
        best_result['densities'] = results

        log_file = f'optimization_results_charge_balanced/best_f{total_error:.1f}.json'
        run_simulation(rate_values, Te, ne, E_field, params_base, log_file)

        print(f"\n  *** NEW BEST: f(x) = {total_error:.2f} at evaluation {objective_function.counter}")
        print(f"      P: {PRESSURE_MTORR:.0f} mTorr, Te: {Te:.2f} eV, Ne: {ne:.2e}, E: {E_field:.1f} V/cm")
        print(f"      Ni/Ne: {Ni_Ne_ratio:.3f} (target: 2-7)")
        print(f"      H: {results['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {results['H']/TARGETS['H']:.2f}x)")
        print(f"      CH: {results['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {results['CH']/TARGETS['CH']:.2f}x)")
        print(f"      C2: {results['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {results['C2']/TARGETS['C2']:.2f}x)")
        print(f"      Species error: {species_error:.2f}, Charge penalty: {charge_penalty:.2f}")
        print(f"      Saved to: {log_file}\n")

    return total_error


def main():
    print("=" * 80)
    print(" OPTIMIZATION WITH CHARGE BALANCE CONSTRAINT")
    print("=" * 80)
    print(f"\nPressure: {PRESSURE_MTORR} mTorr (fixed)")
    print(f"Ne constraint: {NE_MIN:.2e} to {NE_MAX:.2e} cm⁻³ (2.3e9 ± 30%)")
    print("Charge balance: Ni/Ne in [2, 7] range enforced via penalty term")
    print("\nOptimizing ~43 parameters:")
    print("  - Te: 1.0-1.5 eV (literature range for this discharge type)")
    print(f"  - Ne: {NE_MIN:.2e}-{NE_MAX:.2e} cm⁻³")
    print("  - E-field: 10-500 V/cm")
    print("  - ~40 reaction rates")

    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()
    print(f"  Selected {len(tunable_rates)} key rates")

    db = get_complete_rate_database()

    params_base = {
        'E_field': 50,
        'L_discharge': 0.45,
        'ne': 2.3e9,
        'Te': 1.0,
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

    bounds = []
    param_names = []

    for name in tunable_rates:
        if name in db:
            bounds.append((db[name].min, db[name].max))
            param_names.append(name)

    bounds.append((1.0, 1.5))           # Te (eV) - literature range for this discharge type
    param_names.append('Te')

    bounds.append((NE_MIN, NE_MAX))     # Ne constrained!
    param_names.append('Ne')

    bounds.append((10.0, 500.0))        # E-field (V/cm)
    param_names.append('E_field')

    print(f"\n Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm⁻³")
    print(f"  CH: {TARGETS['CH']:.2e} cm⁻³")
    print(f"  C2: {TARGETS['C2']:.2e} cm⁻³")
    print(f"  + Charge balance (Ni ≈ Ne)")

    print("\n" + "=" * 80)
    print(" Running Optimization (20 iterations, pop=12)")
    print("=" * 80)
    print("\nThis enforces charge balance while optimizing species targets!\n")

    start_time = time.time()

    result = differential_evolution(
        objective_function,
        bounds,
        args=(param_names, params_base),
        strategy='best1bin',
        maxiter=20,
        popsize=12,
        tol=0.01,
        atol=0.0,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=42,
        disp=True,
        workers=1,
        updating='deferred',
        polish=True
    )

    elapsed = time.time() - start_time

    print("\n" + "=" * 80)
    print(f" OPTIMIZATION COMPLETE ({elapsed/60:.1f} minutes)")
    print("=" * 80)

    print(f"\nFinal objective: {result.fun:.6f}")
    print(f"Iterations: {result.nit}")
    print(f"Evaluations: {result.nfev}")

    Te_opt = result.x[-3]
    Ne_opt = result.x[-2]
    E_opt = result.x[-1]

    print(f"\nOptimized parameters:")
    print(f"  Te: {Te_opt:.2f} eV")
    print(f"  Ne: {Ne_opt:.2e} cm⁻³")
    print(f"  E-field: {E_opt:.1f} V/cm")

    final_data = {
        'objective': float(result.fun),
        'parameters': {
            'pressure_mTorr': float(PRESSURE_MTORR),
            'Te': float(Te_opt),
            'Ne': float(Ne_opt),
            'E_field': float(E_opt),
            'rates': {name: float(result.x[i]) for i, name in enumerate(param_names[:-3])}
        },
        'success': result.success,
        'message': result.message,
        'iterations': int(result.nit),
        'evaluations': int(result.nfev)
    }

    with open('optimization_results_charge_balanced/FINAL_RESULT.json', 'w') as f:
        json.dump(final_data, f, indent=2)

    print(f"\n✓ Results saved to optimization_results_charge_balanced/")


if __name__ == '__main__':
    main()
