#!/usr/bin/env python3
"""
WALL STICKING FOCUSED OPTIMIZATION

USER CLARIFICATION:
- Only H, CH, C2 were measured (NOT C2H2!)
- Discharge is continuous (steady-state model OK)
- Wall sticking can be tuned within literature ranges

PROBLEM:
Out of 380 previous results, ZERO had all three (H, CH, C2) in 0.6-1.4× range

STRATEGY:
Treat wall sticking coefficients as primary tunable parameters
Goal: Find wall sticking values that achieve H, CH, C2 all in range [0.6-1.4×]
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import time
import signal
from contextlib import contextmanager
import os
import json

from define_rates_tunable import define_rates_tunable
from rate_database_complete import get_complete_rate_database
from build_reactions import build_reactions

# REAL targets (C2H2 NOT measured!)
TARGETS = {
    'H': 5.18e13,
    'CH': 1.0e9,
    'C2': 1.3e11,
}

# Equal high weights - all three equally important
WEIGHTS = {
    'H': 50.0,
    'CH': 50.0,
    'C2': 50.0,
}

# Create results directory
os.makedirs('optimization_results_wall_sticking', exist_ok=True)

# Global counter
best_result = {'objective': 1e10, 'params': None, 'densities': None}


class TimeoutException(Exception):
    pass


@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


class PlasmaODE:
    """Standard ODE with H_drift_gain"""
    
    def __init__(self, params):
        self.species = params['species']
        self.ns = len(self.species)
        self.R = params['R']
        self.nr = len(self.R)
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
        
        self.prod_matrix = sparse.csr_matrix(
            (vals_prod, (rows_prod, cols_prod)),
            shape=(self.ns, self.nr)
        )
        self.react_matrix = sparse.csr_matrix(
            (vals_react, (rows_react, cols_react)),
            shape=(self.ns, self.nr)
        )
    
    def __call__(self, t, y):
        y = np.maximum(y, 1e-10)

        rates = np.zeros(self.nr)
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
        
        production = self.prod_matrix.dot(rates)
        consumption = self.react_matrix.dot(rates)
        dydt = production - consumption
        
        dydt[self.H_idx] += self.H_drift_gain
        
        return dydt


def get_all_sticking_rates():
    """Get all wall sticking rate names from database"""
    db = get_complete_rate_database()
    sticking_rates = [name for name in db.keys() if 'stick_' in name]
    return sticking_rates


def run_simulation(rate_values, E_field, ne, Te, H_drift_gain, params_base):
    """Run plasma simulation"""
    try:
        params = params_base.copy()
        params['E_field'] = E_field
        params['ne'] = ne
        params['Te'] = Te
        params['H_drift_gain'] = H_drift_gain
        
        k = define_rates_tunable(params)
        db = get_complete_rate_database()
        
        # Apply rate values
        for name, val in rate_values.items():
            if name in k and name in db:
                k[name] = np.clip(val, db[name].min, db[name].max)
        
        params['k'] = k
        params['R'], params['tags'] = build_reactions(params)
        
        species = params['species']
        ns = len(species)
        
        # Initial conditions
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
        set_density('H2', 1e12)
        set_density('H', 1e11)
        set_density('C2', 5e7)
        set_density('CH', 5e4)
        set_density('CH3', 5e7)
        
        ode_func = PlasmaODE(params)
        
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
        
        # Calculate charge balance
        ion_species = ['ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus',
                      'CHPlus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'C2HPlus',
                      'H2Plus', 'H3Plus', 'CH2Plus']
        n_i = sum(get_density(sp) for sp in ion_species)
        n_e = get_density('e')
        charge_imbalance = abs(n_i - n_e) / n_e * 100 if n_e > 0 else 999
        
        results = {
            'H': get_density('H'),
            'CH': get_density('CH'),
            'C2': get_density('C2'),
            'C2H2': get_density('C2H2'),
            'n_i': n_i,
            'n_e': n_e,
            'charge_imbalance': charge_imbalance,
        }
        
        return results
        
    except Exception as e:
        return None


def objective_function(x, param_names, params_base):
    """Objective: minimize error for H, CH, C2 (ignore C2H2!)"""
    
    # Track evaluations
    if not hasattr(objective_function, 'counter'):
        objective_function.counter = 0
        objective_function.last_print = 0
    
    objective_function.counter += 1
    
    # Print progress every 10 evaluations
    if objective_function.counter - objective_function.last_print >= 10:
        elapsed = (time.time() - objective_function.start_time) / 60
        print(f"  [{objective_function.counter} evaluations, {elapsed:.1f} min elapsed]")
        objective_function.last_print = objective_function.counter
    
    # Parse parameters
    rate_values = {}
    for i, name in enumerate(param_names[:-4]):
        rate_values[name] = x[i]
    
    E_field = x[-4]
    ne = x[-3]
    Te = x[-2]
    H_drift_gain = x[-1]
    
    # Run simulation
    results = run_simulation(rate_values, E_field, ne, Te, H_drift_gain, params_base)
    
    if results is None:
        return 1e10
    
    # Calculate errors for ONLY H, CH, C2 (NOT C2H2!)
    errors = {}
    total_error = 0.0
    
    for species in ['H', 'CH', 'C2']:
        target = TARGETS[species]
        measured = results[species]
        weight = WEIGHTS[species]
        
        if measured < 1e-10:
            error = 1e10
        else:
            rel_error = abs(measured - target) / target
            error = weight * rel_error ** 2
        
        errors[species] = error
        total_error += error
    
    # Add charge balance penalty
    charge_penalty = 10.0 * (results['charge_imbalance'] / 100.0) ** 2
    total_error += charge_penalty
    
    # Track best result
    global best_result
    if total_error < best_result['objective']:
        best_result['objective'] = total_error
        best_result['params'] = {
            'E_field': E_field,
            'ne': ne,
            'Te': Te,
            'H_drift_gain': H_drift_gain,
            'rates': dict(rate_values)
        }
        best_result['densities'] = results
        
        print(f"\n  *** NEW BEST: f(x) = {total_error:.2f} at evaluation {objective_function.counter}")
        print(f"      H:  {results['H']:.2e} ({results['H']/TARGETS['H']:.2f}×)")
        print(f"      CH: {results['CH']:.2e} ({results['CH']/TARGETS['CH']:.2f}×)")
        print(f"      C2: {results['C2']:.2e} ({results['C2']/TARGETS['C2']:.2f}×)")
        print(f"      Charge: {results['charge_imbalance']:.1f}%")
    
    return total_error


if __name__ == '__main__':
    print("=" * 80)
    print(" WALL STICKING FOCUSED OPTIMIZATION")
    print("=" * 80)
    print()
    print("GOAL: Find wall sticking coefficients that achieve:")
    print("  H, CH, C2 all in range [0.6-1.4×] simultaneously")
    print()
    print("CLARIFICATION: C2H2 was NOT measured experimentally - ignore it!")
    print()
    
    # Get all sticking rates
    all_sticking = get_all_sticking_rates()
    print(f"Found {len(all_sticking)} wall sticking rates")
    
    # Select key sticking rates to tune
    key_sticking = [
        'stick_H_9_1',      # H sticking
        'stick_CH_9_3',     # CH sticking - KEY!
        'stick_C2_9_9',     # C2 sticking
        'stick_C2H2_9_11',  # C2H2 sticking
        'stick_CH3_9_2',    # CH3 sticking
        'stick_C2H4_9_12',  # C2H4 sticking
        'stick_C2H6_9_14',  # C2H6 sticking
        'stick_C_9_10',     # C sticking
        'stick_CH2_9_13',   # CH2 sticking
        'stick_C2H5_9_17',  # C2H5 sticking
    ]
    
    # Build parameter list
    param_names = key_sticking + ['E_field', 'ne', 'Te', 'H_drift_gain']
    
    # Get bounds from database
    db = get_complete_rate_database()
    bounds = []
    
    print()
    print("Tunable parameters:")
    print("-" * 80)
    
    for i, name in enumerate(param_names[:-4]):
        if name in db:
            min_val = db[name].min
            max_val = db[name].max
            bounds.append((min_val, max_val))
            print(f"{i+1:2d}. {name:40s} [{min_val:.2e}, {max_val:.2e}]")
    
    # Plasma parameters
    bounds.append((100, 400))      # E_field (V/cm)
    bounds.append((1e8, 5e9))      # ne (cm^-3)
    bounds.append((0.7, 8.0))      # Te (eV)
    bounds.append((1e16, 5e18))    # H_drift_gain (cm^-3/s)
    
    print(f"{len(param_names)-3:2d}. E_field                                  [100, 400] V/cm")
    print(f"{len(param_names)-2:2d}. ne                                       [1e8, 5e9] cm^-3")
    print(f"{len(param_names)-1:2d}. Te                                       [0.7, 8] eV")
    print(f"{len(param_names):2d}. H_drift_gain                             [1e16, 5e18] cm^-3/s")
    
    print()
    print(f"Total parameters: {len(param_names)}")
    print()
    
    print(" Targets:")
    for species, target in TARGETS.items():
        weight = WEIGHTS[species]
        print(f"  {species:4s}: {target:.2e} cm⁻³  (weight: {weight})")
    print("  + Charge balance penalty: 10.0×")
    print()
    
    # Base parameters
    params_base = {
        'species': [],
        'P': 0.4,  # Torr
        'Tgas': 300,  # K
    }
    
    print("=" * 80)
    print(" Running Optimization (30 iterations, pop=10)")
    print("=" * 80)
    print()
    
    objective_function.start_time = time.time()
    objective_function.counter = 0
    objective_function.last_print = 0
    
    # Run optimization
    result = differential_evolution(
        lambda x: objective_function(x, param_names, params_base),
        bounds,
        maxiter=30,
        popsize=10,
        seed=42,
        disp=True,
        workers=1,
        updating='deferred',
        polish=True
    )
    
    print()
    print("=" * 80)
    print(" OPTIMIZATION COMPLETE")
    print("=" * 80)
    print()
    print(f"Best objective value: {result.fun:.2f}")
    print(f"Total evaluations: {objective_function.counter}")
    print()
    
    if best_result['densities']:
        print("Final solution:")
        dens = best_result['densities']
        print(f"  H:  {dens['H']:.2e} ({dens['H']/TARGETS['H']:.2f}× target)")
        print(f"  CH: {dens['CH']:.2e} ({dens['CH']/TARGETS['CH']:.2f}× target)")
        print(f"  C2: {dens['C2']:.2e} ({dens['C2']/TARGETS['C2']:.2f}× target)")
        print(f"  C2H2: {dens['C2H2']:.2e} (unmeasured)")
        print(f"  Charge imbalance: {dens['charge_imbalance']:.1f}%")
        
        # Check if all three in range
        h_ok = 0.6 <= dens['H']/TARGETS['H'] <= 1.4
        ch_ok = 0.6 <= dens['CH']/TARGETS['CH'] <= 1.4
        c2_ok = 0.6 <= dens['C2']/TARGETS['C2'] <= 1.4
        
        if h_ok and ch_ok and c2_ok:
            print()
            print("★★★ SUCCESS! ALL THREE TARGETS IN RANGE [0.6-1.4×]! ★★★")
        else:
            print()
            print("Status:")
            print(f"  H:  {'✓' if h_ok else '✗'}")
            print(f"  CH: {'✓' if ch_ok else '✗'}")
            print(f"  C2: {'✓' if c2_ok else '✗'}")

