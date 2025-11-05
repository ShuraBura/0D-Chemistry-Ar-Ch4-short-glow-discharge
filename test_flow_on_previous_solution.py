#!/usr/bin/env python3
"""
Test: Take the best previous solution and add flow to see if ODE solver can handle it
"""

import numpy as np
from scipy.integrate import solve_ivp
import json

from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions

# Load previous successful solution
with open('optimization_results_CH_suppression/best_f72.7_Hdrift1.05e+17.json', 'r') as f:
    prev_solution = json.load(f)

# Flow parameters
FLOW_RATE_CM3_S = 2307.4  # 138,445 cm³/min
REACTOR_VOLUME_CM3 = 400.0
RESIDENCE_TIME_S = REACTOR_VOLUME_CM3 / FLOW_RATE_CM3_S  # 0.173 seconds
FLOW_LOSS_RATE = 1.0 / RESIDENCE_TIME_S  # 5.77 /s

print(f"Flow loss rate: {FLOW_LOSS_RATE:.2f} /s")
print(f"Residence time: {RESIDENCE_TIME_S:.3f} s")
print()

# Build chemistry with previous solution's parameters
Te = prev_solution['Te']
Ne = prev_solution['Ne']
E_field = prev_solution['E_field']
H_drift_gain = prev_solution['H_drift_gain']

print(f"Using previous solution parameters:")
print(f"  Te: {Te:.2f} eV")
print(f"  Ne: {Ne:.2e} cm^-3")
print(f"  E_field: {E_field:.1f} V/cm")
print(f"  H_drift_gain: {H_drift_gain:.2e} cm^-3/s")
print()

# Get reactions
all_rates_tunable = define_rates_tunable()
params = build_reactions(all_rates_tunable, Te=Te, Ne=Ne, E_field=E_field)

# Apply rate values from previous solution
for rate_name, rate_value in prev_solution['rate_values'].items():
    if rate_name in params['k']:
        params['k'][rate_name] = rate_value

params['H_drift_gain'] = H_drift_gain

# Test 1: Run WITHOUT flow (should work)
print("="*80)
print("TEST 1: Running WITHOUT flow (baseline)")
print("="*80)

class PlasmaODE_NoFlow:
    def __init__(self, params):
        self.species = params['species']
        self.ns = len(self.species)
        self.R = params['R']
        self.nr = len(self.R)
        self.k = params['k']
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
        self.prod_matrix = sparse.csr_matrix((vals_prod, (rows_prod, cols_prod)), shape=(self.ns, self.nr))
        self.react_matrix = sparse.csr_matrix((vals_react, (rows_react, cols_react)), shape=(self.ns, self.nr))
    
    def __call__(self, t, y):
        y = np.maximum(y, 1e-10)
        rates = np.zeros(self.nr)
        for rxn_idx, reaction in enumerate(self.R):
            rate_constant = self.k[self.tags[rxn_idx]]
            if 'drift' in reaction.tags:
                rates[rxn_idx] = rate_constant * y[self.e_idx]
            elif 'loss' in reaction.tags or 'stick' in reaction.tags:
                react_species = np.where(reaction.reactants > 0)[0]
                if len(react_species) > 0:
                    rates[rxn_idx] = rate_constant * y[react_species[0]]
            else:
                rate = rate_constant
                react_species = np.where(reaction.reactants > 0)[0]
                for sp_idx in react_species:
                    rate *= y[sp_idx] ** reaction.reactants[sp_idx]
                rates[rxn_idx] = rate
        production = self.prod_matrix.dot(rates)
        consumption = self.react_matrix.dot(rates)
        dydt = production - consumption
        dydt[self.H_idx] += self.H_drift_gain
        return dydt

ode_func = PlasmaODE_NoFlow(params)
y0 = np.array([prev_solution['all_densities'].get(sp, 1e8) for sp in params['species']])

try:
    sol = solve_ivp(ode_func, [0, 1e-3], y0, method='BDF', rtol=1e-6, atol=1e-8)
    if sol.success:
        print(f"✓ SUCCESS! Reached t={sol.t[-1]:.2e} s")
        final = sol.y[:, -1]
        H_idx = params['species'].index('H')
        CH_idx = params['species'].index('CH')
        C2_idx = params['species'].index('C2')
        C2H2_idx = params['species'].index('C2H2')
        print(f"  H:    {final[H_idx]:.2e} ({final[H_idx]/5.18e13:.2f}×)")
        print(f"  CH:   {final[CH_idx]:.2e} ({final[CH_idx]/1.0e9:.2f}×)")
        print(f"  C2:   {final[C2_idx]:.2e} ({final[C2_idx]/1.3e11:.2f}×)")
        print(f"  C2H2: {final[C2H2_idx]:.2e} ({final[C2H2_idx]/4.0e12:.2f}×)")
    else:
        print(f"✗ FAILED: {sol.message}")
except Exception as e:
    print(f"✗ EXCEPTION: {e}")

print()

# Test 2: Run WITH flow
print("="*80)
print("TEST 2: Running WITH flow")
print("="*80)

class PlasmaODE_WithFlow:
    def __init__(self, params, flow_loss_rate):
        self.species = params['species']
        self.ns = len(self.species)
        self.R = params['R']
        self.nr = len(self.R)
        self.k = params['k']
        self.H_drift_gain = params.get('H_drift_gain', 3.2e17)
        self.flow_loss_rate = flow_loss_rate
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
        self.prod_matrix = sparse.csr_matrix((vals_prod, (rows_prod, cols_prod)), shape=(self.ns, self.nr))
        self.react_matrix = sparse.csr_matrix((vals_react, (rows_react, cols_react)), shape=(self.ns, self.nr))
    
    def __call__(self, t, y):
        y = np.maximum(y, 1e-10)
        rates = np.zeros(self.nr)
        for rxn_idx, reaction in enumerate(self.R):
            rate_constant = self.k[self.tags[rxn_idx]]
            if 'drift' in reaction.tags:
                rates[rxn_idx] = rate_constant * y[self.e_idx]
            elif 'loss' in reaction.tags or 'stick' in reaction.tags:
                react_species = np.where(reaction.reactants > 0)[0]
                if len(react_species) > 0:
                    rates[rxn_idx] = rate_constant * y[react_species[0]]
            else:
                rate = rate_constant
                react_species = np.where(reaction.reactants > 0)[0]
                for sp_idx in react_species:
                    rate *= y[sp_idx] ** reaction.reactants[sp_idx]
                rates[rxn_idx] = rate
        production = self.prod_matrix.dot(rates)
        consumption = self.react_matrix.dot(rates)
        dydt = production - consumption
        dydt[self.H_idx] += self.H_drift_gain
        
        # ADD FLOW LOSSES (only neutrals, not feedstock or electrons)
        for i, sp in enumerate(self.species):
            if sp not in ['Ar', 'CH4', 'e']:
                dydt[i] -= self.flow_loss_rate * y[i]
        
        return dydt

ode_func_flow = PlasmaODE_WithFlow(params, FLOW_LOSS_RATE)

try:
    sol = solve_ivp(ode_func_flow, [0, 1e-3], y0, method='BDF', rtol=1e-6, atol=1e-8)
    if sol.success:
        print(f"✓ SUCCESS! Reached t={sol.t[-1]:.2e} s")
        final = sol.y[:, -1]
        H_idx = params['species'].index('H')
        CH_idx = params['species'].index('CH')
        C2_idx = params['species'].index('C2')
        C2H2_idx = params['species'].index('C2H2')
        print(f"  H:    {final[H_idx]:.2e} ({final[H_idx]/5.18e13:.2f}×)")
        print(f"  CH:   {final[CH_idx]:.2e} ({final[CH_idx]/1.0e9:.2f}×)")
        print(f"  C2:   {final[C2_idx]:.2e} ({final[C2_idx]/1.3e11:.2f}×)")
        print(f"  C2H2: {final[C2H2_idx]:.2e} ({final[C2H2_idx]/4.0e12:.2f}×)")
        
        print()
        print("Comparing C2H2 losses:")
        stick_rate = prev_solution['rate_values'].get('stick_C2H2_9_11', 0)
        print(f"  Wall sticking rate: {stick_rate:.1f} /s")
        print(f"  Flow loss rate: {FLOW_LOSS_RATE:.2f} /s")
        print(f"  Ratio (wall/flow): {stick_rate/FLOW_LOSS_RATE:.0f}×")
        
    else:
        print(f"✗ FAILED: {sol.message}")
except Exception as e:
    print(f"✗ EXCEPTION: {e}")

