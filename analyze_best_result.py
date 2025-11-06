#!/usr/bin/env python3
"""
Analyze the best result from optimization
Based on step 15, evaluation 760:
  f(x) = 3407.17
  H: 2.75e+13 (0.53×)
  CH: 5.35e+09 (5.35×)
  C2: 2.12e+10 (0.16×)
  C2H2: 6.85e+11
  stick_C2H2 factor: 0.87×
  stick_CH factor: 1.20×
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

# Load starting data
with open('optimization_results_C2H2_boost/best_f259.0_Te1.09.json', 'r') as f:
    start_data = json.load(f)

# Best parameters from step 15:
# From log output, these are the best parameters found
# We need to reconstruct the exact parameters
# stick_C2H2 factor: 0.87×
# stick_CH factor: 1.20×

# Read log to extract best parameters
print("="*80)
print("ANALYZING BEST RESULT FROM OPTIMIZATION")
print("="*80)
print()
print("Best result found at step 15, evaluation 760:")
print("  f(x) = 3407.17")
print("  H:    2.75e+13 (0.53×)")
print("  CH:   5.35e+09 (5.35×)")
print("  C2:   2.12e+10 (0.16×)")
print("  C2H2: 6.85e+11")
print("  stick_C2H2 factor: 0.87×")
print("  stick_CH factor: 1.20×")
print()

# Since we don't have the exact parameters from the optimization,
# let me parse the log file to get them
import re

with open('final_optimization_run.log', 'r') as f:
    log_content = f.read()

# Find the last "NEW BEST" entry with f(x) = 3407
pattern = r'\*\*\* NEW BEST: f\(x\) = 3407\.17[\s\S]*?stick_C2H2 factor: ([\d.]+)×[\s\S]*?stick_CH factor: ([\d.]+)×'
matches = list(re.finditer(pattern, log_content))

if matches:
    last_match = matches[-1]
    stick_C2H2_factor = float(last_match.group(1))
    stick_CH_factor = float(last_match.group(2))
    print(f"Extracted from log:")
    print(f"  stick_C2H2 factor: {stick_C2H2_factor}×")
    print(f"  stick_CH factor: {stick_CH_factor}×")
    print()
else:
    print("Could not extract parameters from log, using values from log output")
    stick_C2H2_factor = 0.87
    stick_CH_factor = 1.20

# We still need the other parameters (Te, ne, H_drift_gain, loss_C2H2)
# The log doesn't show these, so I need to search more carefully
# Let me look for the full parameter set

print("Attempting to reconstruct full parameter set...")
print()
print("NOTE: Log format doesn't include all parameters (Te, ne, H_drift_gain, loss_C2H2)")
print("      These would be in the optimization but not printed in the best result output.")
print()
print("To get the complete breakdown, we need the optimization to finish and save results.")
print()
print("="*80)
