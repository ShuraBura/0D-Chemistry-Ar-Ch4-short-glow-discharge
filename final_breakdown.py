#!/usr/bin/env python3
"""
Extract complete species breakdown from best checkpoint using working simulation code
"""

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
        self.nr = len(self.R)
        self.k = params['k']
        self.tags = params['tags']
        self.H_drift_gain = params.get('H_drift_gain', 3.2e17)
        self.e_idx = self.species.index('e')
        self.Ar_idx = self.species.index('Ar')
        self.CH4_idx = self.species.index('CH4')
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
        y = np.maximum(y, 1e-6)
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
        dydt[self.e_idx] = 0
        dydt[self.Ar_idx] = 0
        dydt[self.CH4_idx] = 0
        return dydt


print("="*80)
print("COMPLETE SPECIES BREAKDOWN FROM CHECKPOINT_F3407.JSON")
print("="*80)
print()

# Load checkpoint
with open('checkpoint_f3407.json', 'r') as f:
    best_result = json.load(f)

params_dict = best_result['params']
species = params_dict['species']

print(f"Best objective: f(x) = {best_result['objective']:.2f}")
print()

# Build reactions exactly as in optimization
db = get_complete_rate_database()
k = define_rates_tunable(params_dict)

for name, val in params_dict.get('rate_values', {}).items():
    if name in k and name in db:
        k[name] = np.clip(val, db[name].min, db[name].max)

params = params_dict.copy()
params['k'] = k
params['R'], params['tags'] = build_reactions(params)

# Initial conditions - same as optimization
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

print("Running simulation to steady state...")
print()

# Solve exactly as in optimization
sol = solve_ivp(PlasmaODE(params), (0, 100), y0, method='BDF', rtol=1e-5, atol=1e-6, max_step=10.0)

if not sol.success:
    print(f"ERROR: Simulation failed - {sol.message}")
    exit(1)

y_final = sol.y[:, -1]

def get_density(name):
    try:
        return y_final[species.index(name)]
    except ValueError:
        return 0.0

# Calculate ion density and charge balance
ion_species = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
               'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus',
               'H3Plus', 'CHPlus', 'H2Plus']
n_i = sum(get_density(sp) for sp in ion_species)
n_e = get_density('e')

print("="*80)
print("ALL SPECIES DENSITIES (cm^-3)")
print("="*80)
print()

for sp in species:
    dens = get_density(sp)
    marker = ""
    if sp in ['H', 'CH', 'C2', 'C2H2']:
        marker = " ← TARGET SPECIES"
    elif sp == 'e':
        marker = " ← Ne (electron density)"
    elif sp in ion_species:
        marker = " (ion)"

    print(f"  {sp:12s}: {dens:.3e}{marker}")

print()
print("="*80)
print("SUMMARY")
print("="*80)
print()
print(f"Total ion density (Ni):  {n_i:.3e} cm^-3")
print(f"Electron density (Ne):   {n_e:.3e} cm^-3")
print(f"Electric field (E):      {params_dict['E_field']:.2f} V/cm")
print(f"Electron temp (Te):      {params_dict['Te']:.3f} eV")
print(f"H drift gain:            {params_dict['H_drift_gain']:.3e} cm^-3/s")
print(f"Pressure (P):            {params_dict['P']} Torr")
print(f"Gas temperature (Tgas):  {params_dict['Tgas']} K")
print()

# Target comparison
TARGETS = {'H': 5.18e13, 'CH': 1.0e9, 'C2': 1.3e11}
print("TARGET SPECIES COMPARISON:")
print(f"  H:    {get_density('H'):.3e} cm^-3  ({get_density('H')/TARGETS['H']:.2f}× target, goal: 0.6-1.4×)")
print(f"  CH:   {get_density('CH'):.3e} cm^-3  ({get_density('CH')/TARGETS['CH']:.2f}× target, goal: 0.6-1.4×)")
print(f"  C2:   {get_density('C2'):.3e} cm^-3  ({get_density('C2')/TARGETS['C2']:.2f}× target, goal: 0.6-1.4×)")
print(f"  C2H2: {get_density('C2H2'):.3e} cm^-3")
print()

# Status
h_ok = 0.6 <= get_density('H')/TARGETS['H'] <= 1.4
ch_ok = 0.6 <= get_density('CH')/TARGETS['CH'] <= 1.4
c2_ok = 0.6 <= get_density('C2')/TARGETS['C2'] <= 1.4

print("STATUS:")
print(f"  H:  {'✓ IN RANGE' if h_ok else '✗ OUT OF RANGE'}")
print(f"  CH: {'✓ IN RANGE' if ch_ok else '✗ OUT OF RANGE'}")
print(f"  C2: {'✓ IN RANGE' if c2_ok else '✗ OUT OF RANGE'}")
print()
print("="*80)
