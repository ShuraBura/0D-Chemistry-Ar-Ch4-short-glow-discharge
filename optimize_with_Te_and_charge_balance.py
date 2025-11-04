#!/usr/bin/env python3
"""
OPTIMIZATION WITH Te-DEPENDENT RATES AND CHARGE BALANCE

Key improvements:
1. Uses define_rates_tunable.py (Te-dependent rates)
2. Te as free optimization parameter (0.7-8 eV)
3. Charge balance penalty in objective function
4. Ne allowed to vary [1e8, 5e9] cm⁻³ (relaxed lower bound)

Expected result:
- Te ~ 5-7 eV (from E/N)
- n_i ~ 2e9 (matches Ne, consistent with sheath width)
- Charge balance < 10%
- Targets matched
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import time
import signal
from contextlib import contextmanager
import os
import json

from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database, get_tunable_rates_for_target
from build_reactions import build_reactions
from odefun_optimized import PlasmaODE_Optimized


# Experimental targets
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

# Create results directory
os.makedirs('optimization_results_Te_low', exist_ok=True)

# Global counter for iterations
iteration_counter = 0
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


def select_tunable_rates():
    """
    Select tunable rates - same as before but for use with define_rates_tunable.
    """
    # Critical neutral reactions
    critical_neutral = [
        'CH2_H_CH_H2_cm3_7_1',          # CH2+H→CH
        'C2_H_CH_C_cm3_7_6',            # C2+H→CH ← KEY!
        'C_H_CH_cm3_7_12',              # C+H→CH
        'C2H2_H_C2_H2_H_cm3_7_50',     # C2H2+H→C2
        'CH_CH4_C2H4_H_cm3_7_20',      # CH+CH4→
        'CH_CH4_CH2_CH3_cm3_7_39',     # CH+CH4→
        'CH_H_C_H2_cm3_7_3',            # CH+H→C+H2
    ]

    db = get_complete_rate_database()

    h_rates = set(get_tunable_rates_for_target('H').keys())
    ch_rates = set(get_tunable_rates_for_target('CH').keys())
    c2_rates = set(get_tunable_rates_for_target('C2').keys())

    target_rates = h_rates | ch_rates | c2_rates

    large_range_rates = {
        name for name, rate in db.items()
        if (rate.max / rate.min > 3.0) and (rate.min > 0)
    }

    flagged_rates = {name for name, rate in db.items() if rate.flag}

    # Combine selections
    selected = set(critical_neutral) | target_rates | large_range_rates | flagged_rates
    selected = {name for name in selected if name in db}

    selected_with_range = [
        (name, db[name].max / db[name].min if db[name].min > 0 else 1.0)
        for name in selected
    ]
    selected_with_range.sort(key=lambda x: x[1], reverse=True)

    # Top 40 rates
    selected_names = [name for name, _ in selected_with_range[:40]]

    return selected_names


def analyze_chemistry(y_final, species, params):
    """Detailed chemistry analysis for important species."""

    # Get reaction rates
    rate_constants = np.array([params['k'][tag] for tag in params['tags']])
    reactions = params['R']

    reaction_rates = rate_constants.copy()
    for rxn_idx, reaction in enumerate(reactions):
        react_species = np.where(reaction.reactants > 0)[0]
        for sp_idx in react_species:
            coeff = reaction.reactants[sp_idx]
            reaction_rates[rxn_idx] *= y_final[sp_idx] ** coeff

    # Analyze important species
    important_species = ['H', 'CH', 'C2', 'C2H2', 'C2H4', 'CH3', 'CH4',
                        'ArStar', 'e', 'ArPlus', 'CH3Plus', 'CH5Plus',
                        'HMinus', 'CH3Minus']

    analysis = {}

    for target_sp in important_species:
        try:
            sp_idx = species.index(target_sp)
        except ValueError:
            continue

        production_rxns = []
        loss_rxns = []

        for rxn_idx, reaction in enumerate(reactions):
            net_change = reaction.products[sp_idx] - reaction.reactants[sp_idx]

            if net_change > 0:
                contribution = net_change * reaction_rates[rxn_idx]
                production_rxns.append((params['tags'][rxn_idx], contribution))
            elif net_change < 0:
                contribution = abs(net_change) * reaction_rates[rxn_idx]
                loss_rxns.append((params['tags'][rxn_idx], contribution))

        production_rxns.sort(key=lambda x: x[1], reverse=True)
        loss_rxns.sort(key=lambda x: x[1], reverse=True)

        total_production = sum(x[1] for x in production_rxns)
        total_loss = sum(x[1] for x in loss_rxns)

        analysis[target_sp] = {
            'density': float(y_final[sp_idx]),
            'production': total_production,
            'loss': total_loss,
            'top_production': [(tag, float(val)) for tag, val in production_rxns[:10]],
            'top_loss': [(tag, float(val)) for tag, val in loss_rxns[:10]]
        }

    return analysis


def run_simulation_with_logging(rate_values, E_field, ne, Te, params_base, log_file=None):
    """Run simulation with Te-dependent rates."""

    try:
        params = params_base.copy()
        params['E_field'] = E_field
        params['ne'] = ne
        params['Te'] = Te  # ← NEW: Pass Te for rate calculation

        k = define_rates_tunable(params)  # ← Use Te-dependent version!
        db = get_complete_rate_database()

        for name, val in rate_values.items():
            if name in k:
                if name in db:
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

        ode_func = PlasmaODE_Optimized(params)

        try:
            with time_limit(30):
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

        y_final = sol.y[:, -1]

        def get_density(name):
            try:
                idx = species.index(name)
                return y_final[idx]
            except ValueError:
                return 0.0

        results = {
            'H': get_density('H'),
            'CH': get_density('CH'),
            'C2': get_density('C2')
        }

        # Calculate ion densities for charge balance
        positive_ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                        'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus',
                        'H3Plus', 'CHPlus', 'H2Plus']

        n_i_total = sum(get_density(ion) for ion in positive_ions)
        results['n_i_total'] = n_i_total
        results['charge_imbalance'] = abs(n_i_total - ne) / ne

        # If logging requested, save detailed analysis
        if log_file:
            all_densities = {species[i]: float(y_final[i]) for i in range(ns)}
            chemistry = analyze_chemistry(y_final, species, params)

            log_data = {
                'Te': Te,
                'Ne': ne,
                'E_field': E_field,
                'rate_values': {k: float(v) for k, v in rate_values.items()},
                'target_densities': {k: results[k] for k in ['H', 'CH', 'C2']},
                'charge_balance': {
                    'n_i_total': float(n_i_total),
                    'n_e': float(ne),
                    'imbalance_percent': float(results['charge_imbalance'] * 100)
                },
                'all_densities': all_densities,
                'chemistry_analysis': chemistry
            }

            with open(log_file, 'w') as f:
                json.dump(log_data, f, indent=2)

        return results

    except Exception as e:
        return None


def objective_function(x, param_names, params_base):
    """Objective function with charge balance penalty."""
    global iteration_counter, best_result

    if not hasattr(objective_function, 'counter'):
        objective_function.counter = 0
        objective_function.start_time = time.time()

    objective_function.counter += 1

    if objective_function.counter % 10 == 0:
        elapsed = time.time() - objective_function.start_time
        print(f"  [{objective_function.counter} evaluations, {elapsed/60:.1f} min elapsed]")

    # Extract parameters
    rate_values = {name: val for name, val in zip(param_names[:-3], x[:-3])}
    E_field = x[-3]
    Te = x[-2]  # ← NEW!
    ne = x[-1]

    params_updated = params_base.copy()
    params_updated['ne'] = ne

    results = run_simulation_with_logging(rate_values, E_field, ne, Te, params_updated)

    if results is None:
        return 1e10

    # Species target errors (same weights as before)
    weights = {
        'H': 1.0,
        'CH': 20.0,
        'C2': 3.0
    }

    species_error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        species_error += weights[species] * rel_error ** 2

    # Charge balance penalty (NEW!)
    charge_penalty = 10.0 * results['charge_imbalance']**2

    # Total error
    total_error = species_error + charge_penalty

    # Save best result with detailed logging
    if total_error < best_result['objective']:
        best_result['objective'] = total_error
        best_result['params'] = {
            'Te': Te,
            'Ne': ne,
            'E_field': E_field,
            'rates': dict(rate_values)
        }
        best_result['densities'] = results

        # Save detailed log for best result
        log_file = f'optimization_results_Te_low/best_f{total_error:.1f}_Te{Te:.2f}.json'
        run_simulation_with_logging(rate_values, E_field, ne, Te, params_updated, log_file)

        print(f"\n  *** NEW BEST: f(x) = {total_error:.2f} at evaluation {objective_function.counter}")
        print(f"      Te: {Te:.2f} eV")
        print(f"      E-field: {E_field:.1f} V/cm")
        print(f"      Ne: {ne:.2e} cm⁻³")
        print(f"      Charge imbalance: {results['charge_imbalance']*100:.2f}%")
        print(f"      H: {results['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {results['H']/TARGETS['H']:.2f}x)")
        print(f"      CH: {results['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {results['CH']/TARGETS['CH']:.2f}x)")
        print(f"      C2: {results['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {results['C2']/TARGETS['C2']:.2f}x)")
        print(f"      Species error: {species_error:.2f}, Charge penalty: {charge_penalty:.2f}")
        print(f"      Saved to: {log_file}\n")

    return total_error


def main():
    global iteration_counter

    print("=" * 80)
    print(" OPTIMIZATION WITH Te-DEPENDENT RATES AND CHARGE BALANCE")
    print("=" * 80)
    print("\nKey improvements:")
    print("  1. Te-dependent rate coefficients")
    print("  2. Te as free parameter (0.7-8 eV)")
    print("  3. Charge balance penalty in objective")
    print("  4. Ne allowed to vary [1e8, 5e9] cm⁻³ (relaxed lower bound)")

    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()
    print(f"  Selected {len(tunable_rates)} key rates")

    db = get_complete_rate_database()

    print("\nSelected rates (top 10):")
    print("-" * 80)
    for i, name in enumerate(tunable_rates[:10], 1):
        if name in db:
            rate = db[name]
            range_factor = rate.max / rate.min if rate.min > 0 else 1.0
            print(f"{i:2d}. {name:45s} [{rate.min:.2e}, {rate.max:.2e}] ({range_factor:.1f}x)")
    print(f"    ... and {len(tunable_rates) - 10} more")

    params_base = {
        'E_field': 50,
        'L_discharge': 0.45,
        'ne': 3.3e9,
        'Te': 1.0,  # Will be optimized
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

    bounds.append((100.0, 400.0))  # E field (V/cm)
    param_names.append('E_field')

    bounds.append((0.7, 8.0))  # Te (eV) ← NEW! Lowered to 0.7 eV
    param_names.append('Te')

    bounds.append((1e8, 5e9))  # Ne (cm⁻³), relaxed lower bound to 1e8
    param_names.append('ne')

    print(f"\n Optimization parameters:")
    print(f"  Tunable rates: {len(param_names) - 3}")
    print(f"  E field: [100, 400] V/cm")
    print(f"  Te: [0.7, 8] eV  ← NEW! Lowered minimum")
    print(f"  Ne: [1e8, 5e9] cm⁻³  ← Relaxed lower bound")
    print(f"  Total parameters: {len(param_names)}")

    print(f"\n Targets:")
    print(f"  H:  {TARGETS['H']:.2e} cm^-3")
    print(f"  CH: {TARGETS['CH']:.2e} cm^-3  (heavily weighted)")
    print(f"  C2: {TARGETS['C2']:.2e} cm^-3")
    print(f"\n  + Charge balance constraint (weight = 10)")

    print("\n" + "=" * 80)
    print(" Running Optimization (15 iterations, pop=6)")
    print("=" * 80)
    print("\nThis will find Te, E, Ne that:")
    print("  1. Produce n_i ≈ Ne (charge balanced)")
    print("  2. Match experimental targets")
    print("  3. Are consistent with sheath width\n")

    start_time = time.time()

    result = differential_evolution(
        objective_function,
        bounds,
        args=(param_names, params_base),
        strategy='best1bin',
        maxiter=15,
        popsize=6,
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

    # Extract final parameters
    Te_opt = result.x[-2]
    Ne_opt = result.x[-1]
    E_opt = result.x[-3]

    print(f"\nOptimized parameters:")
    print(f"  Te: {Te_opt:.2f} eV")
    print(f"  E-field: {E_opt:.1f} V/cm")
    print(f"  Ne: {Ne_opt:.2e} cm⁻³")

    # Calculate E/N for validation
    n_total = 9.66e15  # cm⁻³
    E_over_N_Td = (E_opt / n_total) * 1e17
    print(f"  E/N: {E_over_N_Td:.1f} Td")
    print(f"  Expected Te from E/N: 5-8 eV")
    if 4 <= Te_opt <= 8:
        print(f"  → Te is consistent with E/N! ✓")
    else:
        print(f"  → Te may be outside expected range")

    # Save final results
    final_params = {
        'Te': float(Te_opt),
        'Ne': float(Ne_opt),
        'E_field': float(E_opt),
        'rates': {name: float(result.x[i]) for i, name in enumerate(param_names[:-3])}
    }

    with open('optimization_results_Te_low/FINAL_RESULT.json', 'w') as f:
        json.dump({
            'objective': float(result.fun),
            'parameters': final_params,
            'success': result.success,
            'message': result.message,
            'iterations': int(result.nit),
            'evaluations': int(result.nfev)
        }, f, indent=2)

    print(f"\n✓ Final parameters saved to optimization_results_Te_low/FINAL_RESULT.json")
    print(f"✓ Best result logs saved to optimization_results_Te_low/best_*.json")
    print("\nRe-run with these parameters to validate results!")


if __name__ == '__main__':
    main()
