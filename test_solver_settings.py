#!/usr/bin/env python3
"""
TEST DIFFERENT ODE SOLVER SETTINGS

Goal: Find solver settings that are more robust to parameter changes
"""

import numpy as np
from scipy.integrate import solve_ivp
import json

from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions

# Load the best known solution
with open('optimization_results_C2H2_boost/best_f259.0_Te1.09.json', 'r') as f:
    start_data = json.load(f)


class PlasmaODE:
    def __init__(self, params):
        self.species = params['species']
        self.ns = len(self.species)
        self.R = params['R']
        self.k = params['k']
        self.tags = params['tags']
        self.H_drift_gain = params.get('H_drift_gain', 3.2e17)
        self.e_idx = self.species.index('e')
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
        self.nr = len(self.R)
        self.prod_matrix = sparse.csr_matrix((vals_prod, (rows_prod, cols_prod)), shape=(self.ns, self.nr))
        self.react_matrix = sparse.csr_matrix((vals_react, (rows_react, cols_react)), shape=(self.ns, self.nr))

    def __call__(self, t, y):
        y = np.maximum(y, 1e-10)
        rates = np.zeros(len(self.R))
        for rxn_idx, reaction in enumerate(self.R):
            rate_constant = self.k.get(reaction.rate, 0.0)
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
        return dydt


def run_simulation(params_dict, solver_settings):
    """Run simulation with specified solver settings"""
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

        sol = solve_ivp(PlasmaODE(params), (0, 100), y0, **solver_settings)

        return sol.success, sol.message if not sol.success else "Success", sol.nfev if hasattr(sol, 'nfev') else 0

    except Exception as e:
        return False, str(e), 0


print("="*80)
print("TESTING DIFFERENT SOLVER SETTINGS")
print("="*80)
print()
print("Testing baseline solution with various solver configurations")
print()

# Baseline parameters
baseline_params = {
    'Te': start_data['Te'],
    'ne': start_data['Ne'],  # Note: capital N in JSON
    'E_field': start_data['E_field'],
    'H_drift_gain': start_data.get('H_drift_gain', 3.2e17),
    'rate_values': start_data['rate_values'],
    'L_discharge': 0.45,  # cm (from documentation)
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

# Also test with slightly modified parameters (10% reduction in stick_C2H2)
modified_rates = start_data['rate_values'].copy()
if 'stick_C2H2_9_11' in modified_rates:
    modified_rates['stick_C2H2_9_11'] *= 0.9

modified_params = {
    'Te': start_data['Te'],
    'ne': start_data['Ne'],  # Note: capital N in JSON
    'E_field': start_data['E_field'],
    'H_drift_gain': start_data.get('H_drift_gain', 3.2e17),
    'rate_values': modified_rates,
    'L_discharge': 0.45,  # cm (from documentation)
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

# Test configurations
solver_configs = [
    {
        'name': 'Current (BDF, rtol=1e-5, atol=1e-6)',
        'settings': {
            'method': 'BDF',
            'rtol': 1e-5,
            'atol': 1e-6,
            'max_step': 10.0
        }
    },
    {
        'name': 'Tighter tolerances (BDF, rtol=1e-6, atol=1e-7)',
        'settings': {
            'method': 'BDF',
            'rtol': 1e-6,
            'atol': 1e-7,
            'max_step': 10.0
        }
    },
    {
        'name': 'Looser tolerances (BDF, rtol=1e-4, atol=1e-5)',
        'settings': {
            'method': 'BDF',
            'rtol': 1e-4,
            'atol': 1e-5,
            'max_step': 10.0
        }
    },
    {
        'name': 'Smaller max_step (BDF, max_step=1.0)',
        'settings': {
            'method': 'BDF',
            'rtol': 1e-5,
            'atol': 1e-6,
            'max_step': 1.0
        }
    },
    {
        'name': 'Larger max_step (BDF, max_step=50.0)',
        'settings': {
            'method': 'BDF',
            'rtol': 1e-5,
            'atol': 1e-6,
            'max_step': 50.0
        }
    },
    {
        'name': 'Radau method (implicit)',
        'settings': {
            'method': 'Radau',
            'rtol': 1e-5,
            'atol': 1e-6,
            'max_step': 10.0
        }
    },
    {
        'name': 'LSODA (adaptive method)',
        'settings': {
            'method': 'LSODA',
            'rtol': 1e-5,
            'atol': 1e-6,
            'max_step': 10.0
        }
    },
    {
        'name': 'Very tight + small steps (BDF, rtol=1e-7, atol=1e-8, max_step=0.1)',
        'settings': {
            'method': 'BDF',
            'rtol': 1e-7,
            'atol': 1e-8,
            'max_step': 0.1
        }
    },
]

print("="*80)
print("TEST 1: BASELINE PARAMETERS (known to work)")
print("="*80)
print()

for config in solver_configs:
    success, message, nfev = run_simulation(baseline_params, config['settings'])
    status = "✓" if success else "✗"
    print(f"{status} {config['name']:60s} nfev={nfev:5d}  {message[:50]}")

print()
print("="*80)
print("TEST 2: MODIFIED PARAMETERS (90% stick_C2H2)")
print("="*80)
print()

for config in solver_configs:
    success, message, nfev = run_simulation(modified_params, config['settings'])
    status = "✓" if success else "✗"
    print(f"{status} {config['name']:60s} nfev={nfev:5d}  {message[:50]}")

print()
print("="*80)
print("SUMMARY")
print("="*80)
print()
print("Recommendations:")
print("1. If a solver method works for modified parameters, use it in optimization")
print("2. Looser tolerances may improve robustness at cost of accuracy")
print("3. Smaller max_step may help with stiff systems")
print("4. Alternative methods (Radau, LSODA) may be more robust")
