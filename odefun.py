"""
odefun.py - ODE function for Ar/CH4 plasma simulation
Converted and optimized from MATLAB odefun_diagnostic.m
"""

import numpy as np

class PlasmaODE:
    """ODE function class for plasma simulation with diagnostics."""

    def __init__(self, params):
        self.species = params['species']
        self.ns = len(self.species)
        self.R = params['R']
        self.nr = len(self.R)
        self.k = params['k']
        self.tags = params['tags']
        self.H_drift_gain = 3.2e17

        # Cache species indices
        self.e_idx = self.species.index('e')
        self.Ar_idx = self.species.index('Ar')
        self.CH4_idx = self.species.index('CH4')

        # Pre-allocate arrays for speed
        self.dydt = np.zeros(self.ns)
        self.reaction_rates = np.zeros(self.nr)

        # Pre-compute reaction stoichiometry for faster evaluation
        self._precompute_stoichiometry()

    def _precompute_stoichiometry(self):
        """Pre-compute stoichiometry matrices for faster ODE evaluation."""
        # Create sparse-like representation using lists of indices
        self.reactant_indices = []
        self.reactant_coeffs = []
        self.product_indices = []
        self.product_coeffs = []

        for i, reaction in enumerate(self.R):
            # Find non-zero reactants
            r_idx = np.where(reaction.reactants > 0)[0]
            r_coeff = reaction.reactants[r_idx]
            self.reactant_indices.append(r_idx)
            self.reactant_coeffs.append(r_coeff)

            # Find non-zero products
            p_idx = np.where(reaction.products > 0)[0]
            p_coeff = reaction.products[p_idx]
            self.product_indices.append(p_idx)
            self.product_coeffs.append(p_coeff)

    def __call__(self, t, y):
        """
        ODE function for scipy.integrate.solve_ivp

        Parameters:
        -----------
        t : float
            Current time
        y : ndarray
            Current species densities

        Returns:
        --------
        dydt : ndarray
            Time derivatives of species densities
        """
        # Prevent negative densities
        y = np.maximum(y, 1e-6)

        # Compute reaction rates efficiently
        for i in range(self.nr):
            rate = self.k[self.tags[i]]
            r_idx = self.reactant_indices[i]
            r_coeff = self.reactant_coeffs[i]

            # Multiply by reactant concentrations
            for idx, coeff in zip(r_idx, r_coeff):
                rate *= y[idx]**coeff

            self.reaction_rates[i] = rate

        # Compute dydt using vectorized operations
        self.dydt.fill(0.0)
        for i in range(self.nr):
            rate = self.reaction_rates[i]

            # Subtract reactants
            for idx, coeff in zip(self.reactant_indices[i], self.reactant_coeffs[i]):
                self.dydt[idx] -= coeff * rate

            # Add products
            for idx, coeff in zip(self.product_indices[i], self.product_coeffs[i]):
                self.dydt[idx] += coeff * rate

        # Add H drift source term
        H_idx = self.species.index('H')
        self.dydt[H_idx] += self.H_drift_gain

        # Fix constant species
        self.dydt[self.e_idx] = 0
        self.dydt[self.Ar_idx] = 0
        self.dydt[self.CH4_idx] = 0

        return self.dydt
