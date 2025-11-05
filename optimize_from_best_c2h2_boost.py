#!/usr/bin/env python3
"""
TARGETED OPTIMIZATION FROM BEST H~0.49×, C2~0.46×, CH~8× RESULT

USER'S INSIGHT: Boost C2H2 to simultaneously increase H and C2

STRATEGY:
Start from best_f259.0_Te1.09.json and:
1. REDUCE stick_C2H2 (allow more C2H2 to accumulate)
2. REDUCE loss_C2H2 (allow more C2H2 to accumulate)
3. INCREASE stick_CH (suppress CH from C2 + H → CH coupling)
4. Fine-tune other parameters

Expected: H and C2 boost via C2H2 + H → C2 + H2 + H
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import time
import os
import json

from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions

TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

WEIGHTS = {
    'H': 50.0,
    'CH': 50.0,
    'C2': 50.0,
}

os.makedirs('optimization_results_c2h2_boost_targeted', exist_ok=True)

best_result = {'objective': 1e10, 'params': None, 'densities': None}

# Load starting point
with open('optimization_results_C2H2_boost/best_f259.0_Te1.09.json', 'r') as f:
    start_data = json.load(f)

print("="*80)
print("TARGETED C2H2 BOOST FROM BEST RESULT")
print("="*80)
print()
print("Starting from: optimization_results_C2H2_boost/best_f259.0_Te1.09.json")
print("  H:  2.53e+13 (0.49×)")
print("  C2: 6.04e+10 (0.46×)")
print("  CH: 8.03e+09 (8.03×)")
print()
print("STRATEGY:")
print("  1. Reduce stick_C2H2 (↑ C2H2)")
print("  2. Reduce loss_C2H2 (↑ C2H2)")
print("  3. Increase stick_CH (↓ CH despite C2 + H → CH coupling)")
print("  4. Result: C2H2 + H → C2 + H2 + H boosts BOTH H and C2!")
print()


class PlasmaODE:
    def __init__(self, params):
        self.species = params['species']
        self.ns = len(self.species)
        self.R = params['R']
        self.nr = len(self.R)  # FIX: was missing!
        self.k = params['k']
        self.tags = params['tags']
        self.H_drift_gain = params.get('H_drift_gain', 3.2e17)
        self.e_idx = self.species.index('e')
        self.Ar_idx = self.species.index('Ar')
        self.CH4_idx = self.species.index('CH4')
        self.H_idx = self.species.index('H')
        self._build_stoich()
    
    def _build_stoich(self):
        from scipy import sparse
        rows_prod, cols_prod, vals_prod = [], [], []
        rows_react, cols_react, vals_react = [], [], []
        for rxn_idx, reaction in enumerate(self.R):
            prod_species = np.where(reaction.products > 0)[0]
            for sp_idx in prod_species:
                rows_prod.append(sp_idx)
                cols_prod.append(rxn_idx)
                vals_prod.append(reaction.products[sp_idx])
            react_species = np.where(reaction.reactants > 0)[0]
            for sp_idx in react_species:
                rows_react.append(sp_idx)
                cols_react.append(rxn_idx)
                vals_react.append(reaction.reactants[sp_idx])
        self.prod_matrix = sparse.csr_matrix((vals_prod, (rows_prod, cols_prod)), shape=(self.ns, self.nr))
        self.react_matrix = sparse.csr_matrix((vals_react, (rows_react, cols_react)), shape=(self.ns, self.nr))
    
    def __call__(self, t, y):
        # Prevent negative densities (higher floor = more stable)
        y = np.maximum(y, 1e-6)

        rates = np.zeros(len(self.R))
        for rxn_idx, reaction in enumerate(self.R):
            rate_constant = self.k[self.tags[rxn_idx]]  # FIX: use tags, not reaction.rate!
            if 'drift' in self.tags[rxn_idx]:
                rates[rxn_idx] = rate_constant * y[self.e_idx]
            elif 'loss' in self.tags[rxn_idx] or 'stick' in self.tags[rxn_idx]:
                react_species = np.where(reaction.reactants > 0)[0]
                if len(react_species) > 0:
                    rates[rxn_idx] = rate_constant * y[react_species[0]]
            else:
                rate = rate_constant
                react_species = np.where(reaction.reactants > 0)[0]
                for sp_idx in react_species:
                    rate *= y[sp_idx] ** reaction.reactants[sp_idx]
                rates[rxn_idx] = rate

        dydt = self.prod_matrix.dot(rates) - self.react_matrix.dot(rates)
        dydt[self.H_idx] += self.H_drift_gain

        # Fix constant species (prevent drift)
        dydt[self.e_idx] = 0
        dydt[self.Ar_idx] = 0
        dydt[self.CH4_idx] = 0

        return dydt


def run_simulation(params_dict):
    try:
        params = {'species': [], 'P': 0.4, 'Tgas': 300}
        params.update(params_dict)
        
        k = define_rates_tunable(params)
        db = get_complete_rate_database()
        
        for name, val in params_dict.get('rate_values', {}).items():
            if name in k and name in db:
                k[name] = np.clip(val, db[name].min, db[name].max)
        
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
    
    # Parse parameters
    stick_C2H2_factor = x[0]  # Multiply baseline by this
    loss_C2H2_factor = x[1]
    stick_CH_factor = x[2]
    Te = x[3]
    ne = x[4]
    H_drift_gain = x[5]
    
    # Build rate values starting from baseline
    rate_values = start_data['rate_values'].copy()
    
    # Adjust C2H2 losses (reduce to boost C2H2)
    if 'stick_C2H2_9_11' in rate_values:
        rate_values['stick_C2H2_9_11'] *= stick_C2H2_factor
    if 'loss_C2H2_11_19' in rate_values:
        rate_values['loss_C2H2_11_19'] *= loss_C2H2_factor
    
    # Adjust CH sticking (increase to suppress CH)
    if 'stick_CH_9_3' in rate_values:
        rate_values['stick_CH_9_3'] *= stick_CH_factor
    
    params_dict = {
        'Te': Te,
        'ne': ne,
        'E_field': start_data['E_field'],  # Keep baseline
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
            'H2Plus': 5000,
            'CH3Minus': 3000,
            'HMinus': 3000
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
        
        print(f"\n  *** NEW BEST: f(x) = {total_error:.2f}")
        print(f"      H:    {results['H']:.2e} ({results['H']/TARGETS['H']:.2f}×)")
        print(f"      CH:   {results['CH']:.2e} ({results['CH']/TARGETS['CH']:.2f}×)")
        print(f"      C2:   {results['C2']:.2e} ({results['C2']/TARGETS['C2']:.2f}×)")
        print(f"      C2H2: {results['C2H2']:.2e}")
        print(f"      stick_C2H2 factor: {stick_C2H2_factor:.2f}×")
        print(f"      stick_CH factor: {stick_CH_factor:.2f}×")
    
    return total_error


# Optimization bounds
bounds = [
    (0.2, 1.0),      # stick_C2H2_factor (reduce 20-100% of baseline)
    (0.3, 1.0),      # loss_C2H2_factor (reduce 30-100% of baseline)
    (1.0, 2.0),      # stick_CH_factor (increase 100-200% of baseline)
    (0.8, 2.0),      # Te (eV) - narrow range around baseline 1.09
    (1e8, 5e9),      # ne
    (1e16, 5e18),    # H_drift_gain
]

print("Optimization bounds:")
print("  stick_C2H2 factor: [0.2, 1.0]× (reduce to boost C2H2)")
print("  loss_C2H2 factor:  [0.3, 1.0]× (reduce to boost C2H2)")
print("  stick_CH factor:   [1.0, 2.0]× (increase to suppress CH)")
print("  Te:  [0.8, 2.0] eV")
print("  ne:  [1e8, 5e9] cm^-3")
print("  H_drift_gain: [1e16, 5e18] cm^-3/s")
print()

print("="*80)
print("RUNNING OPTIMIZATION (20 iterations, pop=8)")
print("="*80)
print()

result = differential_evolution(
    objective_function,
    bounds,
    maxiter=20,
    popsize=8,
    seed=42,
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
