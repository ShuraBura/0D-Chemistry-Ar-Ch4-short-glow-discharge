#!/usr/bin/env python3
"""
CHARGE-BALANCED BREAKTHROUGH V6 - Missing CH4+ recombination fixed!

CRITICAL DISCOVERY: CH4+ (38.4% of all ions!) has NO recombination in V5!

KEY IMPROVEMENTS V5 → V6:
1. ALL 15 dissociative recombination reactions NOW TUNABLE (V5 had only 1!)
   - CH4+ + e → CH3 + H          (Literature: 3.08e-7 cm³/s)
   - CH4+ + e → CH2 + H2         (Literature: 9.74e-7 cm³/s) ← DOMINANT!
   - CH4+ + e → CH + H2 + H      (Literature: 3.93e-7 cm³/s)
   - CH4+ + e → C + 2H2          (Literature: 3.42e-8 cm³/s)
   - Plus 11 more ion + e reactions

2. CHARGE BALANCE WEIGHT: 75 → 300 (4x increase!)
   - V5 had Ion/e = 0.524x (target 4.5x) - catastrophic 88.4% deficit
   - Electron imbalance: 97.8% (production 45x loss!)
   - Root cause: Missing CH4+ recombination (major electron sink)

3. Literature-informed rates from Thomas et al. (2013)
   - CH4+ total recombination rate: 1.71×10⁻⁶ cm³/s at 300K
   - Database updated with correct branching fractions

Physics: Ion + e → neutrals is THE critical electron sink!
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


# Experimental targets
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

# Charge balance target (ion/electron ratio)
TARGET_ION_E_RATIO = 4.5  # Middle of 3-6x range
CHARGE_BALANCE_WEIGHT = 300.0  # INCREASED from 75 (4x!) to fix charge crisis

# Create results directory
os.makedirs('optimization_results_charge_balanced_v6', exist_ok=True)

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
    V6: CHARGE-BALANCED - ALL dissociative recombination reactions!

    CRITICAL FIX: V5 had only 1/15 dissociative recombination reactions tunable!
    - CH4+ (38.4% of ions) had ZERO recombination → catastrophic charge imbalance
    - V5 Ion/e ratio: 0.524x (target 4.5x) - 88.4% deficit!
    - V5 Electron imbalance: 97.8% (production 45x loss)

    V6 adds ALL 15 ion + e → neutrals reactions:
    - 4 CH4+ recombination channels (literature-updated rates)
    - 3 CH5+ recombination channels
    - 8 other ion + e reactions (C2Hx+, ArPlus, etc.)

    Strategy: Comprehensive tuning PLUS explicit dissociative recombination!
    """
    db = get_complete_rate_database()

    # Start with empty set
    selected = set()

    # 1. ALL rates with range > 1.5x (same as V5)
    tunable_range = {
        name for name, rate in db.items()
        if (rate.max / rate.min > 1.5) and (rate.min > 0)
    }
    selected |= tunable_range

    # 2. ALL wall sticking rates
    wall_sticking = {name for name in db.keys() if name.startswith('stick_')}
    selected |= wall_sticking

    # 3. ALL volume loss rates
    volume_loss = {name for name in db.keys() if name.startswith('loss_')}
    selected |= volume_loss

    # 4. ALL electron-impact EXCEPT ionization (same as V5)
    ionization_keywords = ['CH4Plus', 'CH3Plus', 'ArPlus', 'CH5Plus',
                          'CHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                          'C2H3Plus', 'C2HPlus', 'H3Plus']

    electron_impact = {
        name for name in db.keys()
        if name.startswith('e_') and
        not any(ion in name for ion in ionization_keywords)
    }
    selected |= electron_impact

    # 5. **NEW IN V6**: ALL DISSOCIATIVE RECOMBINATION REACTIONS (Ion+ + e → neutrals)
    #    This is THE FIX for the charge balance crisis!
    dissoc_recomb = [
        # CH4+ recombination (CRITICAL - 38.4% of all ions!)
        'e_CH4Plus_CH3_H_cm3_6_4',       # CH4+ + e → CH3 + H (18% branch, 3.08e-7)
        'e_CH4Plus_CH2_H2_cm3_6_9',      # CH4+ + e → CH2 + H2 (51%+6% branch, 9.74e-7) ← DOMINANT!
        'e_CH4Plus_CH_H2_H_cm3_6_11',    # CH4+ + e → CH + H2 + H (23% branch, 3.93e-7)
        'e_CH4Plus_C_2H2_cm3_6_13',      # CH4+ + e → C + 2H2 (2% branch, 3.42e-8)

        # CH5+ recombination (26.6% of all ions)
        'CH5Plus_e_CH4_H_cm3_6_3',       # CH5+ + e → CH4 + H
        'CH5Plus_e_CH3_H2_cm3_6_8',      # CH5+ + e → CH3 + H2
        'CH5Plus_e_CH2_H2_H_cm3_6_10',   # CH5+ + e → CH2 + H2 + H
        'CH5Plus_e_CH3_2H_cm3_6_12',     # CH5+ + e → CH3 + 2H

        # Other ion + e recombination
        'ArPlus_e_Ar_cm3_6_1',           # Ar+ + e → Ar (11% of ions)
        'CH3Plus_e_CH3_cm3_6_2',         # CH3+ + e → CH3 (6.6% of ions)
        'C2H5Plus_e_C2H4_H_cm3_6_14',    # C2H5+ + e → C2H4 + H (5.5% of ions)
        'C2H4Plus_e_C2H2_H2_cm3_6_15',   # C2H4+ + e → C2H2 + H2
        'C2H3Plus_e_C2H2_H_cm3_6_16',    # C2H3+ + e → C2H2 + H
        'C2HPlus_e_C2_H_cm3_6_18',       # C2H+ + e → C2 + H (was only one in V5!)
        'C2H5Plus_e_C2H4_H_cm3_6_28',    # C2H5+ + e → C2H4 + H (duplicate?)
    ]
    selected |= set(dissoc_recomb)

    # 6. CRITICAL SECONDARY PATHWAYS (from V5 chemistry analysis)
    critical_secondary = [
        # H chemistry
        'H_CH4_CH3_H2_cm3_7_25',       # H+CH4→CH3+H2 (75.4% of H loss!)
        'stick_H_9_1',                  # H wall sticking (23.4% of H loss)
        'CH3_H_CH2_H2_cm3_7_36',       # CH3+H→CH2+H2

        # C2 chemistry
        'stick_C2_9_9',                 # C2 wall sticking (63.4% of C2 loss!)
        'C2_H_CH_C_cm3_7_6',           # C2+H→CH+C (18.5% of C2 loss)
        'loss_C2_11_3',                 # C2 volume loss (18.1% of C2 loss!)
        'C2H2_H_C2_H2_H_cm3_7_50',     # C2H2+H→C2 (95.5% of C2 production)
        'CH2_CH2_C2_H2_H2_cm3_7_58',   # CH2+CH2→C2

        # CH chemistry
        'CH2_H_CH_H2_cm3_7_1',          # CH2+H→CH (54.2% of CH production)
        'C_H_CH_cm3_7_12',              # C+H→CH
        'CH_CH4_C2H4_H_cm3_7_20',      # CH+CH4→C2H4
        'CH_CH4_CH2_CH3_cm3_7_39',     # CH+CH4→CH2+CH3
        'loss_CH_11_9',                 # CH volume loss
        'stick_CH_9_3',                 # CH wall sticking
    ]
    selected |= set(critical_secondary)

    # 7. Target-specific rates
    h_rates = set(get_tunable_rates_for_target('H').keys())
    ch_rates = set(get_tunable_rates_for_target('CH').keys())
    c2_rates = set(get_tunable_rates_for_target('C2').keys())
    selected |= (h_rates | ch_rates | c2_rates)

    # 8. Flagged rates (uncertain in literature)
    flagged_rates = {name for name, rate in db.items() if rate.flag}
    selected |= flagged_rates

    # Filter to only rates that exist in database
    selected = {name for name in selected if name in db}

    # Count dissociative recombination reactions
    dissoc_count = sum(1 for name in selected if name in dissoc_recomb)

    # Sort by range (largest first) for reporting
    selected_with_range = [
        (name, db[name].max / db[name].min if db[name].min > 0 else 1.0)
        for name in selected
    ]
    selected_with_range.sort(key=lambda x: x[1], reverse=True)

    selected_names = [name for name, _ in selected_with_range]

    print(f"  V6 CHARGE-BALANCED selection: {len(selected_names)} rates")
    print(f"    - Tunable range (>1.5x): {len(tunable_range)}")
    print(f"    - Wall sticking: {len(wall_sticking)}")
    print(f"    - Volume loss: {len(volume_loss)}")
    print(f"    - Electron-impact (no ionization): {len(electron_impact)}")
    print(f"    - ★ DISSOCIATIVE RECOMBINATION: {dissoc_count} (V5 had only 1!)")
    print(f"    - Critical secondary pathways: {len(critical_secondary)}")

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
                        'CH4Plus', 'HMinus', 'CH3Minus']  # Added CH4Plus!

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


def run_simulation_with_logging(rate_values, E_field, ne, params_base, log_file=None):
    """Run simulation and optionally log detailed results."""

    try:
        params = params_base.copy()
        params['E_field'] = E_field
        params['ne'] = ne

        k = define_rates(params)
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

        # Get all species densities for charge balance calculation
        all_densities = {species[i]: float(y_final[i]) for i in range(ns)}

        results = {
            'H': get_density('H'),
            'CH': get_density('CH'),
            'C2': get_density('C2'),
            'all_densities': all_densities  # Include all for charge balance
        }

        # If logging requested, save detailed analysis
        if log_file:
            chemistry = analyze_chemistry(y_final, species, params)

            log_data = {
                'Ne': ne,
                'E_field': E_field,
                'rate_values': {k: float(v) for k, v in rate_values.items()},
                'target_densities': {
                    'H': results['H'],
                    'CH': results['CH'],
                    'C2': results['C2']
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
    """Objective function with ENHANCED charge balance penalty."""
    global iteration_counter, best_result

    if not hasattr(objective_function, 'counter'):
        objective_function.counter = 0
        objective_function.start_time = time.time()

    objective_function.counter += 1

    if objective_function.counter % 10 == 0:
        elapsed = time.time() - objective_function.start_time
        print(f"  [{objective_function.counter} evaluations, {elapsed/60:.1f} min elapsed]")

    rate_values = {name: val for name, val in zip(param_names[:-2], x[:-2])}
    E_field = x[-2]
    ne = x[-1]

    params_updated = params_base.copy()
    params_updated['ne'] = ne

    results = run_simulation_with_logging(rate_values, E_field, ne, params_updated)

    if results is None:
        return 1e10

    # Extract target species errors
    weights = {
        'H': 3.0,    # H is 0.40x, need to increase
        'CH': 30.0,  # CH is critical, balance carefully
        'C2': 5.0    # C2 is 0.37x, need to increase
    }

    error = 0.0
    for species in ['H', 'CH', 'C2']:
        rel_error = (results[species] - TARGETS[species]) / TARGETS[species]
        error += weights[species] * rel_error ** 2

    # ADD CHARGE BALANCE PENALTY (4x stronger than V5!)
    all_dens = results['all_densities']

    # Calculate total positive ions
    positive_ions = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                     'H3Plus', 'CHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus',
                     'C2H3Plus', 'C2HPlus']
    total_positive = sum(all_dens.get(ion, 0) for ion in positive_ions)

    # Ion/electron ratio (target: 4.5x, range 3-6x)
    ion_e_ratio = total_positive / ne if ne > 0 else 0

    # Charge balance error (penalize deviation from target ratio)
    charge_error = abs(ion_e_ratio - TARGET_ION_E_RATIO) / TARGET_ION_E_RATIO

    # Add weighted charge penalty to objective (300x vs V5's 75x!)
    error += CHARGE_BALANCE_WEIGHT * charge_error

    # Save best result with detailed logging
    if error < best_result['objective']:
        best_result['objective'] = error
        best_result['params'] = {
            'Ne': ne,
            'E_field': E_field,
            'rates': dict(rate_values)
        }
        best_result['densities'] = results

        # Save detailed log for best result
        log_file = f'optimization_results_charge_balanced_v6/best_iteration_{iteration_counter:04d}_f{error:.1f}.json'
        run_simulation_with_logging(rate_values, E_field, ne, params_updated, log_file)

        print(f"\n  *** NEW BEST: f(x) = {error:.2f} at evaluation {objective_function.counter}")
        print(f"      H: {results['H']:.2e} (target: {TARGETS['H']:.2e}, ratio: {results['H']/TARGETS['H']:.2f}x)")
        print(f"      CH: {results['CH']:.2e} (target: {TARGETS['CH']:.2e}, ratio: {results['CH']/TARGETS['CH']:.2f}x)")
        print(f"      C2: {results['C2']:.2e} (target: {TARGETS['C2']:.2e}, ratio: {results['C2']/TARGETS['C2']:.2f}x)")
        print(f"      Ion/e ratio: {ion_e_ratio:.2f}x (target: {TARGET_ION_E_RATIO:.1f}x)")
        print(f"      Saved to: {log_file}\n")

    return error


def create_warm_start_population(param_names, bounds, db, popsize=6):
    """
    Create initial population biased for charge balance AND species targets.

    V6 strategy: Start near V5 best but allow CH4+ recombination rates to vary widely.
    """
    print("\n  Creating warm-start population for CHARGE BALANCE breakthrough...")

    population = np.zeros((popsize, len(bounds)))

    # CH4+ recombination reactions - allow FULL range (critical for charge balance!)
    ch4_recomb = {
        'e_CH4Plus_CH3_H_cm3_6_4',
        'e_CH4Plus_CH2_H2_cm3_6_9',
        'e_CH4Plus_CH_H2_H_cm3_6_11',
        'e_CH4Plus_C_2H2_cm3_6_13',
    }

    # CH production reactions - bias toward MINIMUM (same as V5)
    ch_production = {
        'CH2_H_CH_H2_cm3_7_1',
        'C2_H_CH_C_cm3_7_6',
        'C_H_CH_cm3_7_12',
    }

    # CH loss reactions - bias toward MAXIMUM (same as V5)
    ch_loss = {
        'stick_CH_9_3',
        'CH_CH4_C2H4_H_cm3_7_20',
        'CH_CH4_CH2_CH3_cm3_7_39',
        'CH_H_C_H2_cm3_7_3',
        'loss_CH_11_9',
    }

    for i in range(popsize):
        for j, (name, (low, high)) in enumerate(zip(param_names, bounds)):
            if name == 'E_field':
                # Bias toward 50-100 V/cm from V5
                population[i, j] = np.random.uniform(max(50, low), min(100, high))
            elif name == 'ne':
                # Bias toward V5 best: Ne ~ 1.4e9
                population[i, j] = np.random.uniform(max(1.2e9, low), 1.6e9)
            elif name in ch4_recomb:
                # CH4+ recombination - FULL RANGE exploration (this is the breakthrough!)
                population[i, j] = np.random.uniform(low, high)
            elif name in ch_production:
                # Minimize CH production
                range_width = high - low
                population[i, j] = np.random.uniform(low, low + 0.2 * range_width)
            elif name in ch_loss:
                # Maximize CH loss
                range_width = high - low
                population[i, j] = np.random.uniform(high - 0.3 * range_width, high)
            else:
                # Random within bounds for other rates
                population[i, j] = np.random.uniform(low, high)

    print("  ✓ Population biased: E~50-100 V/cm, Ne~1.4e9, CH4+ recomb FULL RANGE")
    return population


def main():
    global iteration_counter

    print("=" * 80)
    print(" V6 OPTIMIZATION - CHARGE BALANCE BREAKTHROUGH!")
    print("=" * 80)
    print("\nCRITICAL DISCOVERY from V5 analysis:")
    print("  ✗ CH4+ (38.4% of all ions!) has NO recombination being tuned!")
    print("  ✗ Only 1/15 dissociative recombination reactions tunable")
    print("  ✗ Ion/electron ratio: 0.524x (target 4.5x) - 88.4% deficit!")
    print("  ✗ Electron imbalance: 97.8% (production 45x loss)")
    print("\nV6 STRATEGY - THE FIX:")
    print("  ✓ ALL 15 dissociative recombination reactions NOW TUNABLE")
    print("  ✓ CH4+ + e rates updated to literature values (Thomas et al. 2013)")
    print("  ✓ Charge balance weight: 75 → 300 (4x increase!)")
    print("  ✓ Expected Ion/e improvement: 0.5x → 2-3x (major breakthrough!)")
    print("  ✓ Expected electron balance improvement: 97.8% → ~40-50%")
    print("\nV5 result: H=0.40x, CH=5.49x, C2=0.37x, Ion/e=0.524x")
    print("V6 goal: Fix charge balance AND maintain species accuracy!")

    print("\nSelecting tunable rates...")
    tunable_rates = select_tunable_rates()
    print(f"  Selected {len(tunable_rates)} key rates")

    db = get_complete_rate_database()

    print("\nKey dissociative recombination rates (NEW IN V6):")
    print("-" * 80)
    ch4_recomb_names = ['e_CH4Plus_CH3_H_cm3_6_4', 'e_CH4Plus_CH2_H2_cm3_6_9',
                        'e_CH4Plus_CH_H2_H_cm3_6_11', 'e_CH4Plus_C_2H2_cm3_6_13']
    for name in ch4_recomb_names:
        if name in db:
            rate = db[name]
            print(f"  {name:45s} [{rate.min:.2e}, {rate.max:.2e}] - {rate.source}")

    params_base = {
        'E_field': 50,
        'L_discharge': 0.45,
        'ne': 3.3e9,
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                    'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                    'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                    'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'C2H2Star'],
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6, 'CH5Plus': 4761.6,
            'ArHPlus': 2969.6, 'CH2Plus': 4949.6, 'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6,
            'C2H3Plus': 4949.6, 'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'CH3Minus': 3000, 'HMinus': 3000
        }
    }

    bounds = []
    param_names = []

    for name in tunable_rates:
        if name in db:
            bounds.append((db[name].min, db[name].max))
            param_names.append(name)

    bounds.append((5.0, 300.0))  # E field
    param_names.append('E_field')

    bounds.append((1.0e9, 5.0e9))  # Ne
    param_names.append('ne')

    print(f"\nOptimization parameters:")
    print(f"  Tunable rates: {len(param_names) - 2}")
    print(f"  E field: [5, 300] V/cm")
    print(f"  Ne: [1.0e9, 5.0e9] cm^-3")
    print(f"  Total parameters: {len(param_names)}")

    print(f"\nTargets:")
    print(f"  H:  {TARGETS['H']:.2e} cm^-3")
    print(f"  CH: {TARGETS['CH']:.2e} cm^-3")
    print(f"  C2: {TARGETS['C2']:.2e} cm^-3")
    print(f"  Ion/e ratio: {TARGET_ION_E_RATIO:.1f}x (HEAVILY weighted: {CHARGE_BALANCE_WEIGHT}x)")

    print("\n" + "=" * 80)
    print(" Running Charge-Balanced Optimization (20 iterations, pop=6)")
    print("=" * 80)
    print("\nDetailed logs will be saved to optimization_results_charge_balanced_v6/")
    print("This is the BREAKTHROUGH run - fixing the charge crisis!\n")

    # Create warm-start population
    init_pop = create_warm_start_population(param_names, bounds, db, popsize=6)

    start_time = time.time()

    result = differential_evolution(
        objective_function,
        bounds,
        args=(param_names, params_base),
        strategy='best1bin',
        maxiter=20,
        popsize=6,
        tol=0.01,
        atol=0.0,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=42,
        init=init_pop,
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

    # Save final results
    final_params = {
        'Ne': result.x[-1],
        'E_field': result.x[-2],
        'rates': {name: result.x[i] for i, name in enumerate(param_names[:-2])}
    }

    with open('optimization_results_charge_balanced_v6/FINAL_RESULT.json', 'w') as f:
        json.dump({
            'objective': result.fun,
            'parameters': final_params,
            'success': result.success,
            'message': result.message,
            'iterations': result.nit,
            'evaluations': result.nfev
        }, f, indent=2)

    print(f"\n✓ Final parameters saved to optimization_results_charge_balanced_v6/FINAL_RESULT.json")
    print(f"✓ Best result logs saved to optimization_results_charge_balanced_v6/best_iteration_*.json")
    print("\nThis should be the BREAKTHROUGH - proper charge balance at last!")


if __name__ == '__main__':
    main()
