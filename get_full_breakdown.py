#!/usr/bin/env python3
"""
Extract complete species breakdown from best checkpoint
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

        # Fix constant species
        dydt[self.e_idx] = 0
        dydt[self.Ar_idx] = 0
        dydt[self.CH4_idx] = 0

        return dydt


# Load best checkpoint
with open('checkpoint_f3407.json', 'r') as f:
    best_result = json.load(f)

with open('optimization_results_C2H2_boost/best_f259.0_Te1.09.json', 'r') as f:
    start_data = json.load(f)

print("="*80)
print("COMPLETE SPECIES BREAKDOWN FROM BEST CHECKPOINT")
print("="*80)
print()
print(f"Best objective: f(x) = {best_result['objective']:.2f}")
print()

params_final = best_result['params']
species = params_final['species']

# Run final simulation
rate_db = get_complete_rate_database()
k = define_rates_tunable(params_final)

# Update rate values
for name, val in params_final.get('rate_values', {}).items():
    if name in k and name in rate_db:
        k[name] = np.clip(val, rate_db[name].min, rate_db[name].max)

params_final['k'] = k
R, tags = build_reactions(params_final)

params_ode = {
    'species': species,
    'R': R,
    'k': k,
    'tags': tags,
    'H_drift_gain': params_final['H_drift_gain']
}

ode = PlasmaODE(params_ode)

# Initial conditions
y0 = np.zeros(len(species))
for sp in species:
    if sp in start_data['all_densities']:
        y0[species.index(sp)] = start_data['all_densities'][sp]

print("Running final simulation...")
# Solve
sol = solve_ivp(
    ode,
    (0, 1e-3),
    y0,
    method='BDF',
    dense_output=True,
    rtol=1e-6,
    atol=1e-10
)

if sol.success:
    y_final = sol.y[:, -1]

    print("SUCCESS!")
    print()
    print("All species densities (cm^-3):")
    print("-" * 80)

    total_ions = 0.0
    ion_species = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                   'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus',
                   'H3Plus', 'CHPlus', 'H2Plus']

    for i, sp in enumerate(species):
        density = y_final[i]
        marker = ""
        if sp in ['H', 'CH', 'C2', 'C2H2']:
            marker = " ← TARGET"
        elif sp == 'e':
            marker = " ← Ne"
        elif sp in ion_species:
            total_ions += density
            marker = " (ion)"

        print(f"  {sp:12s}: {density:.3e}{marker}")

    print("-" * 80)
    print()
    print("SUMMARY:")
    print(f"  Total ion density (Ni): {total_ions:.3e} cm^-3")
    print(f"  Electron density (Ne):  {y_final[species.index('e')]:.3e} cm^-3")
    print(f"  Electric field (E):     {params_final['E_field']:.2f} V/cm")
    print(f"  Electron temp (Te):     {params_final['Te']:.3f} eV")
    print(f"  H drift gain:           {params_final['H_drift_gain']:.3e} cm^-3/s")
    print()
    print("TARGET SPECIES:")
    print(f"  H:    {y_final[species.index('H')]:.3e} cm^-3 ({y_final[species.index('H')]/5.18e13:.2f}× target)")
    print(f"  CH:   {y_final[species.index('CH')]:.3e} cm^-3 ({y_final[species.index('CH')]/1.0e9:.2f}× target)")
    print(f"  C2:   {y_final[species.index('C2')]:.3e} cm^-3 ({y_final[species.index('C2')]/1.3e11:.2f}× target)")
    print(f"  C2H2: {y_final[species.index('C2H2')]:.3e} cm^-3")
    print()
else:
    print("ERROR: Final simulation failed to converge")
    print(f"Message: {sol.message}")
    print()
