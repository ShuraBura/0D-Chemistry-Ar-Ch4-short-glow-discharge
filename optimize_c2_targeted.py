#!/usr/bin/env python3
"""
NEW OPTIMIZATION: Target C2 problem directly
Based on pathway analysis findings:
  - C2 consumption: 59% wall sticking (at min), 28% H+C2→CH+C, 13% loss_C2
  - Add tuning for H+C2→CH+C and loss_C2 reactions
"""

import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from odefun import PlasmaODE

# Load starting point
with open('checkpoint_f3407.json', 'r') as f:
    start_data_checkpoint = json.load(f)

# Also load original for initial densities
with open('optimization_results_C2H2_boost/best_f259.0_Te1.09.json', 'r') as f:
    start_data = json.load(f)

print("="*80)
print("TARGETED C2 OPTIMIZATION - NEW STRATEGY")
print("="*80)
print()
print("Starting from: checkpoint_f3407.json (f(x) = 3407.17)")
print(f"  H:  2.75e+13 (0.53×)")
print(f"  C2: 2.12e+10 (0.16×) ← PRIMARY TARGET")
print(f"  CH: 5.35e+09 (5.35×)")
print()

print("NEW STRATEGY based on pathway analysis:")
print("  C2 consumption breakdown:")
print("    59% - wall sticking (stick_C2 at literature MIN, can't reduce)")
print("    28% - H + C2 → CH + C (NEW: tune this DOWN)")
print("    13% - loss_C2 (NEW: tune this DOWN)")
print()
print("  CH consumption:")
print("    81% - CH + CH4 → C2H4 (chemical, hard to tune)")
print("     2% - wall sticking (stick_CH already ABOVE literature MAX)")
print()

# Targets
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11
}

# Weights (increase C2 weight since it's the worst)
WEIGHTS = {
    'H': 1.0,
    'CH': 3.0,
    'C2': 5.0  # Highest priority
}

# Track best result
best_result = {
    'objective': float('inf'),
    'params': None,
    'densities': None,
    'x': None
}

def run_simulation(params_dict):
    """Run single simulation and return target densities."""
    try:
        db = get_complete_rate_database()
        k = define_rates_tunable(params_dict)

        # Apply rate modifications from params_dict
        for name, val in params_dict.get('rate_values', {}).items():
            if name in k and name in db:
                k[name] = np.clip(val, db[name].min, db[name].max)

        params = params_dict.copy()
        params['k'] = k
        params['R'], params['tags'] = build_reactions(params)

        species = params['species']
        y0 = np.ones(len(species)) * 1e3

        def set_density(name, value):
            try:
                y0[species.index(name)] = value
            except ValueError:
                pass

        set_density('e', params['ne'])
        set_density('Ar', 0.85 * 9.66e15)
        set_density('CH4', 0.15 * 9.66e15)
        set_density('ArPlus', 1e7)
        set_density('CH4Plus', 1e5)
        set_density('H2', 1e12)
        set_density('H', 1e11)
        set_density('C2', 5e7)
        set_density('CH', 5e4)
        set_density('CH3', 5e7)

        sol = solve_ivp(PlasmaODE(params), (0, 100), y0, method='BDF', rtol=1e-5, atol=1e-6, max_step=10.0)

        if not sol.success:
            return None

        y_final = sol.y[:, -1]

        def get_density(name):
            try:
                return y_final[species.index(name)]
            except ValueError:
                return 0.0

        ion_species = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus']
        n_i = sum(get_density(sp) for sp in ion_species)
        n_e = get_density('e')
        charge_imbalance = abs(n_i - n_e) / n_e * 100 if n_e > 0 else 999

        return {
            'H': get_density('H'),
            'CH': get_density('CH'),
            'C2': get_density('C2'),
            'C2H2': get_density('C2H2'),
            'charge_imbalance': charge_imbalance,
        }
    except:
        return None


def objective_function(x):
    if not hasattr(objective_function, 'counter'):
        objective_function.counter = 0
        objective_function.last_print = 0
    objective_function.counter += 1

    if objective_function.counter - objective_function.last_print >= 10:
        print(f"  [{objective_function.counter} evaluations]")
        objective_function.last_print = objective_function.counter

    # Parse parameters (8 total)
    C2_H_CH_C_factor = x[0]      # NEW: H + C2 → CH + C rate
    loss_C2_factor = x[1]         # NEW: C2 wall loss rate
    stick_C2H2_factor = x[2]      # C2H2 wall sticking
    loss_C2H2_factor = x[3]       # C2H2 loss rate
    stick_CH_factor = x[4]        # CH wall sticking (locked at 1.0)
    Te = x[5]
    ne = x[6]
    H_drift_gain = x[7]

    # Build rate values from checkpoint
    rate_values = start_data_checkpoint['params']['rate_values'].copy()

    # Apply NEW C2 reaction tuning
    if 'C2_H_CH_C_cm3_7_6' in rate_values:
        rate_values['C2_H_CH_C_cm3_7_6'] *= C2_H_CH_C_factor

    # Apply NEW C2 loss tuning
    if 'loss_C2_11_3' in rate_values:
        rate_values['loss_C2_11_3'] *= loss_C2_factor

    # Apply C2H2 tuning
    if 'stick_C2H2_9_11' in rate_values:
        rate_values['stick_C2H2_9_11'] *= stick_C2H2_factor

    if 'loss_C2H2_11_19' in rate_values:
        rate_values['loss_C2H2_11_19'] *= loss_C2H2_factor

    # Apply CH tuning (locked at baseline since already above max)
    if 'stick_CH_9_3' in rate_values:
        rate_values['stick_CH_9_3'] *= stick_CH_factor

    params_dict = {
        'Te': Te,
        'ne': ne,
        'E_field': start_data['E_field'],
        'H_drift_gain': H_drift_gain,
        'rate_values': rate_values,
        'species': ['e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                    'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2',
                    'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H',
                    'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',
                    'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus', 'H2Plus', 'C2H2Star'],
        'L_discharge': 0.45,
        'P': 0.4,
        'Tgas': 300,
        'mobilities': {
            'ArPlus': 3057.28, 'CH4Plus': 6432, 'CH3Plus': 4949.6,
            'CH5Plus': 4761.6, 'ArHPlus': 2969.6, 'CH2Plus': 4949.6,
            'C2H5Plus': 4949.6, 'C2H4Plus': 4949.6, 'C2H3Plus': 4949.6,
            'C2HPlus': 5000, 'H3Plus': 5000, 'CHPlus': 5000,
            'H2Plus': 5000, 'CH3Minus': 3000, 'HMinus': 3000
        }
    }

    results = run_simulation(params_dict)

    if results is None:
        return 1e10

    total_error = 0.0
    for species in ['H', 'CH', 'C2']:
        target = TARGETS[species]
        measured = results[species]
        weight = WEIGHTS[species]
        if measured < 1e-10:
            return 1e10
        rel_error = abs(measured - target) / target
        total_error += weight * rel_error ** 2

    charge_penalty = 10.0 * (results['charge_imbalance'] / 100.0) ** 2
    total_error += charge_penalty

    global best_result
    if total_error < best_result['objective']:
        best_result['objective'] = total_error
        best_result['params'] = params_dict
        best_result['densities'] = results
        best_result['x'] = x.tolist()

        print(f"\n  *** NEW BEST: f(x) = {total_error:.2f}")
        print(f"      H:    {results['H']:.2e} ({results['H']/TARGETS['H']:.2f}×)")
        print(f"      CH:   {results['CH']:.2e} ({results['CH']/TARGETS['CH']:.2f}×)")
        print(f"      C2:   {results['C2']:.2e} ({results['C2']/TARGETS['C2']:.2f}×)")
        print(f"      C2H2: {results['C2H2']:.2e}")
        print(f"      C2_H_CH_C factor: {C2_H_CH_C_factor:.2f}×")
        print(f"      loss_C2 factor: {loss_C2_factor:.2e}×")

        # Save checkpoint
        checkpoint_file = f"checkpoint_c2opt_f{total_error:.0f}.json"
        with open(checkpoint_file, 'w') as f:
            json.dump(best_result, f, indent=2)
        print(f"      Saved: {checkpoint_file}")

    return total_error


# Optimization bounds (8 parameters)
bounds = [
    (0.83, 1.25),    # C2_H_CH_C factor (covers literature min/max)
    (0.0001, 2.5),   # loss_C2 factor (HUGE range available!)
    (0.2, 1.0),      # stick_C2H2 factor (reduce to boost C2H2)
    (0.3, 1.0),      # loss_C2H2 factor (reduce to boost C2H2)
    (1.0, 1.0),      # stick_CH factor (LOCKED - already above lit max)
    (0.8, 2.0),      # Te (eV)
    (1e8, 5e9),      # ne
    (1e16, 5e18),    # H_drift_gain
]

print("Optimization parameters (8 total):")
print("  1. C2_H_CH_C factor:   [0.83, 1.25]× (NEW: reduce C2 loss)")
print("  2. loss_C2 factor:     [0.0001, 2.5]× (NEW: HUGE range!)")
print("  3. stick_C2H2 factor:  [0.2, 1.0]× (reduce to boost C2H2)")
print("  4. loss_C2H2 factor:   [0.3, 1.0]× (reduce to boost C2H2)")
print("  5. stick_CH factor:    [1.0, 1.0]× (LOCKED at baseline)")
print("  6. Te:                 [0.8, 2.0] eV")
print("  7. ne:                 [1e8, 5e9] cm^-3")
print("  8. H_drift_gain:       [1e16, 5e18] cm^-3/s")
print()

print("="*80)
print("RUNNING OPTIMIZATION (15 iterations, pop=10)")
print("="*80)
print()

result = differential_evolution(
    objective_function,
    bounds,
    maxiter=15,
    popsize=10,
    seed=43,  # Different seed from previous
    disp=True,
    workers=1,
    updating='deferred',
    polish=True
)

print()
print("="*80)
print("OPTIMIZATION COMPLETE")
print("="*80)
print()
print(f"Best objective: {result.fun:.2f}")
print(f"Total evaluations: {objective_function.counter}")
print()

if best_result['densities']:
    dens = best_result['densities']
    print("Final solution:")
    print(f"  H:    {dens['H']:.2e} ({dens['H']/TARGETS['H']:.2f}×)")
    print(f"  CH:   {dens['CH']:.2e} ({dens['CH']/TARGETS['CH']:.2f}×)")
    print(f"  C2:   {dens['C2']:.2e} ({dens['C2']/TARGETS['C2']:.2f}×)")
    print(f"  C2H2: {dens['C2H2']:.2e}")

    h_ok = 0.6 <= dens['H']/TARGETS['H'] <= 1.4
    ch_ok = 0.6 <= dens['CH']/TARGETS['CH'] <= 1.4
    c2_ok = 0.6 <= dens['C2']/TARGETS['C2'] <= 1.4

    print()
    if h_ok and ch_ok and c2_ok:
        print("★★★ SUCCESS! ALL THREE IN RANGE [0.6-1.4×]! ★★★")
    else:
        print("Status:")
        print(f"  H:  {'✓' if h_ok else '✗'}")
        print(f"  CH: {'✓' if ch_ok else '✗'}")
        print(f"  C2: {'✓' if c2_ok else '✗'}")

    print()
    x = best_result['x']
    print("Optimized parameters:")
    print(f"  C2_H_CH_C factor: {x[0]:.4f}×")
    print(f"  loss_C2 factor:   {x[1]:.6f}×")
    print(f"  stick_C2H2:       {x[2]:.4f}×")
    print(f"  loss_C2H2:        {x[3]:.4f}×")
    print(f"  stick_CH:         {x[4]:.4f}× (locked)")
    print(f"  Te:               {x[5]:.3f} eV")
    print(f"  ne:               {x[6]:.3e} cm^-3")
    print(f"  H_drift_gain:     {x[7]:.3e} cm^-3/s")

print()
print("="*80)
