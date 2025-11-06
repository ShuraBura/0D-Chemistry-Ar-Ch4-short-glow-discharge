#!/usr/bin/env python3
"""
Extract best result from optimization log and run final simulation
to get complete species breakdown.
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


print("="*80)
print("FINAL OPTIMIZATION RESULTS")
print("="*80)
print()

# Load starting data
with open('optimization_results_C2H2_boost/best_f259.0_Te1.09.json', 'r') as f:
    start_data = json.load(f)

print("Best result from log (step 19, evaluation ~930):")
print("  f(x) = 1260.63")
print("  H:    2.63e+13 (0.51×)")
print("  CH:   2.77e+09 (2.77×)")
print("  C2:   1.84e+10 (0.14×)")
print("  C2H2: 6.17e+11")
print("  stick_C2H2 factor: 0.87×")
print("  stick_CH factor: 1.08×")
print()

# Best parameters from step 19 log output
# The log shows stick_C2H2=0.87 and stick_CH=1.08
# But we don't have the other 4 parameters (loss_C2H2, Te, ne, H_drift_gain)
# These were optimized but not printed in the "NEW BEST" output

# The script optimizes 6 parameters:
# x[0] = stick_C2H2_factor
# x[1] = loss_C2H2_factor
# x[2] = stick_CH_factor
# x[3] = Te
# x[4] = ne
# x[5] = H_drift_gain

# Since we only have 2 of the 6 parameters from the log,
# we can't reconstruct the exact simulation.

print("NOTE: The log only prints stick_C2H2 and stick_CH factors.")
print("      The optimization also tuned: loss_C2H2, Te, ne, H_drift_gain")
print("      These parameters are not shown in the 'NEW BEST' output.")
print()
print("To get the full species breakdown, the optimization would need to")
print("save the complete parameter set and run a final simulation.")
print()

# Check if there's a saved result file
import os
result_files = [f for f in os.listdir('.') if f.startswith('best_') and f.endswith('.json')]
if result_files:
    print("Found result files:")
    for f in result_files:
        print(f"  {f}")
else:
    print("No saved result files found.")
    print()
    print("The optimization stopped at step 19 before completing and saving results.")

print()
print("="*80)
