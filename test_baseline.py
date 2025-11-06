#!/usr/bin/env python3
"""Quick test to see if baseline simulation works"""
import json
import numpy as np
from scipy.integrate import solve_ivp
from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions
from scipy import sparse

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
        self.prod_matrix = sparse.csr_matrix((vals_prod, (rows_prod, cols_prod)), shape=(self.ns, len(self.R)))
        self.react_matrix = sparse.csr_matrix((vals_react, (rows_react, cols_react)), shape=(self.ns, len(self.R)))

    def __call__(self, t, y):
        y = np.maximum(y, 1e-10)
        rates = np.zeros(len(self.R))
        for rxn_idx, reaction in enumerate(self.R):
            rate_constant = self.k[self.tags[rxn_idx]]
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

# Load start data
with open('optimization_results_C2H2_boost/best_f259.0_Te1.09.json', 'r') as f:
    start_data = json.load(f)

params_dict = {
    'Te': start_data['Te'],
    'ne': start_data['Ne'],
    'E_field': start_data['E_field'],
    'H_drift_gain': start_data.get('H_drift_gain', 3.2e17),
    'rate_values': start_data['rate_values'],
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

print(f'Calling solve_ivp...')
sol = solve_ivp(PlasmaODE(params), (0, 100), y0, method='BDF', rtol=1e-5, atol=1e-6, max_step=10.0)

print(f'Success: {sol.success}')
print(f'Message: {sol.message}')
print(f'nfev: {sol.nfev}')

if sol.success:
    y_final = sol.y[:, -1]
    H = y_final[species.index('H')]
    CH = y_final[species.index('CH')]
    C2 = y_final[species.index('C2')]
    print(f'H:  {H:.2e}')
    print(f'CH: {CH:.2e}')
    print(f'C2: {C2:.2e}')
