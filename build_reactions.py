"""
build_reactions.py - Build reaction network for Ar/CH4 plasma simulation
Converted from MATLAB build_reactions.m
"""

import numpy as np

class Reaction:
    """Class to represent a single chemical reaction."""
    def __init__(self, reactants, products, rate, tag):
        self.reactants = reactants
        self.products = products
        self.rate = rate
        self.tag = tag

def build_reactions(params):
    """
    Build the reaction network for the plasma simulation.

    Parameters:
    -----------
    params : dict
        Dictionary containing species list and rate constants

    Returns:
    --------
    R : list of Reaction objects
    tags : list of str
        Reaction tags
    """
    species = params['species']
    ns = len(species)
    k = params['k']

    R = []
    tags = []

    def sto(*args):
        """Create stoichiometry vector for given species and coefficients."""
        vec = np.zeros(ns)
        for i in range(0, len(args), 2):
            species_name = args[i]
            coeff = args[i+1]
            try:
                idx = species.index(species_name)
                vec[idx] = coeff
            except ValueError:
                raise ValueError(f'Species {species_name} not found')
        return vec

    def push(reactants, products, rate, tag):
        """Add a reaction to the list."""
        R.append(Reaction(reactants, products, rate, tag))
        tags.append(tag)

    # Group 1: Electron-Impact Reactions (Neutral Products)
    push(sto('e', 1, 'CH4', 1), sto('e', 1, 'CH3', 1, 'H', 1), k['e_CH4_CH3_H_cm3_1_1'], 'e_CH4_CH3_H_cm3_1_1')
    push(sto('e', 1, 'CH4', 1), sto('e', 1, 'CH2', 1, 'H2', 1), k['e_CH4_CH2_H2_cm3_1_2'], 'e_CH4_CH2_H2_cm3_1_2')
    push(sto('e', 1, 'CH4', 1), sto('e', 1, 'CH', 1, 'H2', 1, 'H', 1), k['e_CH4_CH_H2_H_vib_cm3_1_3'], 'e_CH4_CH_H2_H_vib_cm3_1_3')
    push(sto('e', 1, 'H2', 1), sto('e', 1, 'H', 2), k['e_H2_H_H_cm3_1_4'], 'e_H2_H_H_cm3_1_4')
    push(sto('e', 1, 'CH3', 1), sto('e', 1, 'CH2', 1, 'H', 1), k['e_CH3_CH2_H_cm3_1_5'], 'e_CH3_CH2_H_cm3_1_5')
    push(sto('e', 1, 'C2H4', 1), sto('e', 1, 'C2H2', 1, 'H2', 1), k['e_C2H4_C2H2_H2_cm3_1_6'], 'e_C2H4_C2H2_H2_cm3_1_6')
    push(sto('e', 1, 'Ar', 1), sto('e', 1, 'ArStar', 1), k['e_Ar_ArStar_cm3_1_7'], 'e_Ar_ArStar_cm3_1_7')
    push(sto('e', 1, 'C2H6', 1), sto('e', 1, 'C2H4', 1, 'H2', 1), k['e_C2H6_C2H4_H2_cm3_1_8'], 'e_C2H6_C2H4_H2_cm3_1_8')
    push(sto('e', 1, 'C2H6', 1), sto('e', 1, 'C2H4', 1, 'H2', 1), k['e_C2H6_C2H4_H2_e_cm3_1_9'], 'e_C2H6_C2H4_H2_e_cm3_1_9')
    push(sto('e', 1, 'CH4', 1), sto('CH3Minus', 1, 'H', 1), k['e_CH4_CH3Minus_H_cm3_1_10'], 'e_CH4_CH3Minus_H_cm3_1_10')
    push(sto('e', 1, 'CH4', 1), sto('e', 1, 'CH', 1, 'H', 1, 'H2', 1), k['e_CH4_CH_H_H2_cm3_1_11'], 'e_CH4_CH_H_H2_cm3_1_11')
    push(sto('e', 1, 'CH', 1), sto('e', 1, 'C', 1, 'H', 1), k['e_CH_CH_C_H_e_cm3_1_12'], 'e_CH_CH_C_H_e_cm3_1_12')
    push(sto('e', 1, 'H2', 1), sto('HMinus', 1, 'H', 1), k['e_H2_HMinus_H_cm3_1_13'], 'e_H2_HMinus_H_cm3_1_13')
    push(sto('e', 1, 'CH3', 1), sto('CH3Minus', 1), k['e_CH3_CH3Minus_cm3_1_14'], 'e_CH3_CH3Minus_cm3_1_14')
    push(sto('e', 1, 'C2H4', 1), sto('C2H2', 1, 'H2', 1), k['e_C2H4_C2H2_H2_cm3_1_15'], 'e_C2H4_C2H2_H2_cm3_1_15')
    push(sto('e', 1, 'C2H2', 1), sto('C2', 1, 'H2', 1), k['e_C2H2_C2_H2_cm3_1_16'], 'e_C2H2_C2_H2_cm3_1_16')
    push(sto('e', 1, 'C2H4', 1), sto('C2H2', 1, 'H', 2), k['e_C2H4_C2H2_H_H_cm3_1_17'], 'e_C2H4_C2H2_H_H_cm3_1_17')
    push(sto('e', 1, 'C2H6', 1), sto('C2H2', 1, 'H2', 2), k['e_C2H6_C2H2_2H2_cm3_1_18'], 'e_C2H6_C2H2_2H2_cm3_1_18')
    # New reactions from audit
    push(sto('e', 1, 'C2H2', 1), sto('e', 1, 'C2H', 1, 'H', 1), k['e_C2H2_C2H_H_cm3_1_19'], 'e_C2H2_C2H_H_cm3_1_19')
    push(sto('e', 1, 'C2H4', 1), sto('e', 1, 'C2H3', 1, 'H', 1), k['e_C2H4_C2H3_H_cm3_1_20'], 'e_C2H4_C2H3_H_cm3_1_20')
    push(sto('e', 1, 'C2H6', 1), sto('e', 1, 'C2H5', 1, 'H', 1), k['e_C2H6_C2H5_H_cm3_1_21'], 'e_C2H6_C2H5_H_cm3_1_21')

    # Group 2: Electron-Impact Ionization
    push(sto('e', 1, 'CH4', 1), sto('CH3Plus', 1, 'H', 1, 'e', 2), k['e_CH4_CH3Plus_H_cm3_2_1'], 'e_CH4_CH3Plus_H_cm3_2_1')
    push(sto('e', 1, 'CH4', 1), sto('e', 2, 'CH4Plus', 1), k['e_CH4_CH4Plus_cm3_2_2'], 'e_CH4_CH4Plus_cm3_2_2')
    push(sto('e', 1, 'Ar', 1), sto('e', 2, 'ArPlus', 1), k['e_Ar_ArPlus_cm3_2_3'], 'e_Ar_ArPlus_cm3_2_3')
    push(sto('e', 1, 'ArStar', 1), sto('e', 2, 'ArPlus', 1), k['e_ArStar_ArPlus_cm3_2_4'], 'e_ArStar_ArPlus_cm3_2_4')
    push(sto('e', 1, 'C2H6', 1), sto('e', 2, 'C2H5Plus', 1, 'H', 1), k['e_C2H6_C2H5Plus_H_2e_cm3_2_5'], 'e_C2H6_C2H5Plus_H_2e_cm3_2_5')
    push(sto('e', 1, 'C2H4', 1), sto('e', 2, 'C2H4Plus', 1), k['e_C2H4_C2H4Plus_2e_cm3_2_6'], 'e_C2H4_C2H4Plus_2e_cm3_2_6')
    push(sto('e', 1, 'C2H4', 1), sto('e', 2, 'C2H3Plus', 1, 'H', 1), k['e_C2H4_C2H3Plus_H_2e_cm3_2_7'], 'e_C2H4_C2H3Plus_H_2e_cm3_2_7')
    push(sto('e', 1, 'C2H2', 1), sto('e', 2, 'C2HPlus', 1, 'H', 1), k['e_C2H2_C2HPlus_2e_cm3_2_8'], 'e_C2H2_C2HPlus_2e_cm3_2_8')
    push(sto('e', 1, 'H2', 1), sto('e', 2, 'H2Plus', 1), k['e_H2_H2Plus_2e_cm3_2_9'], 'e_H2_H2Plus_2e_cm3_2_9')
    push(sto('e', 1, 'CH', 1), sto('e', 2, 'CHPlus', 1), k['e_CH_CHPlus_2e_cm3_2_10'], 'e_CH_CHPlus_2e_cm3_2_10')  # CH ionization (NEW!)

    # Group 3: ArStar Reactions
    push(sto('ArStar', 1, 'CH4', 1), sto('Ar', 1, 'CH3', 1, 'H', 1), k['ArStar_CH4_CH3_H_cm3_3_1'], 'ArStar_CH4_CH3_H_cm3_3_1')
    push(sto('ArStar', 1, 'H2', 1), sto('Ar', 1, 'H', 2), k['ArStar_H2_H_H_cm3_3_2'], 'ArStar_H2_H_H_cm3_3_2')
    push(sto('ArStar', 1, 'H2', 1), sto('ArHPlus', 1, 'H', 1, 'e', 1), k['ArStar_H2_ArHPlus_H_cm3_3_3'], 'ArStar_H2_ArHPlus_H_cm3_3_3')
    push(sto('ArStar', 1, 'CH2', 1), sto('Ar', 1, 'CH', 1, 'H', 1), k['ArStar_CH2_CH_H_cm3_3_4'], 'ArStar_CH2_CH_H_cm3_3_4')
    push(sto('ArStar', 1, 'C', 1), sto('Ar', 1, 'C', 1), k['ArStar_C_Ar_CStar_cm3_3_5'], 'ArStar_C_Ar_CStar_cm3_3_5')
    push(sto('ArStar', 1, 'H', 1), sto('ArHPlus', 1, 'e', 1), k['ArStar_H_ArHPlus_e_cm3_3_6'], 'ArStar_H_ArHPlus_e_cm3_3_6')
    push(sto('ArStar', 1, 'CH3', 1), sto('Ar', 1, 'CH2', 1, 'H', 1), k['ArStar_CH3_Ar_CH2_H_cm3_3_7'], 'ArStar_CH3_Ar_CH2_H_cm3_3_7')
    push(sto('ArStar', 1, 'C2', 1), sto('Ar', 1, 'C2', 1), k['ArStar_C2_Ar_C2_cm3_3_8'], 'ArStar_C2_Ar_C2_cm3_3_8')
    push(sto('ArStar', 1, 'C2H4', 1), sto('Ar', 1, 'C2H4', 1), k['ArStar_C2H4_Ar_C2H4_cm3_3_9'], 'ArStar_C2H4_Ar_C2H4_cm3_3_9')
    push(sto('ArStar', 1, 'CH2', 1), sto('Ar', 1, 'CH2', 1), k['ArStar_CH2_Ar_CH2_cm3_3_10'], 'ArStar_CH2_Ar_CH2_cm3_3_10')
    push(sto('ArStar', 1, 'CH3', 1), sto('Ar', 1, 'CH3', 1), k['ArStar_CH3_Ar_CH3_cm3_3_11'], 'ArStar_CH3_Ar_CH3_cm3_3_11')
    push(sto('ArStar', 1), sto('Ar', 1), k['ArStar_M_Ar_3_12'], 'ArStar_M_Ar_3_12')
    push(sto('ArStar', 1, 'Ar', 1), sto('Ar', 2), k['ArStar_Ar_Ar_Ar_cm3_3_13'], 'ArStar_Ar_Ar_Ar_cm3_3_13')
    push(sto('CH', 1, 'ArStar', 1), sto('C', 1, 'H', 1, 'Ar', 1), k['CH_ArStar_C_H_Ar_cm3_3_14'], 'CH_ArStar_C_H_Ar_cm3_3_14')
    push(sto('CH', 1, 'ArStar', 1), sto('CH', 1, 'Ar', 1), k['CH_ArStar_CH_Ar_cm3_3_15'], 'CH_ArStar_CH_Ar_cm3_3_15')
    push(sto('ArStar', 1, 'CH3', 1), sto('Ar', 1, 'CH2', 1, 'H', 1), k['ArStar_CH3_CH2_H_Ar_cm3_3_16'], 'ArStar_CH3_CH2_H_Ar_cm3_3_16')
    push(sto('ArStar', 1, 'e', 1), sto('Ar', 1, 'e', 1), k['ArStar_e_Ar_e_cm3_3_17'], 'ArStar_e_Ar_e_cm3_3_17')
    push(sto('ArStar', 1, 'H2', 1), sto('Ar', 1, 'H2', 1), k['ArStar_H2_Ar_H2Star_cm3_3_18'], 'ArStar_H2_Ar_H2Star_cm3_3_18')
    push(sto('ArStar', 1, 'CH4', 1), sto('Ar', 1, 'CH4', 1), k['ArStar_CH4_Ar_CH4Star_cm3_3_19'], 'ArStar_CH4_Ar_CH4Star_cm3_3_19')
    push(sto('ArStar', 1, 'CH3', 1), sto('Ar', 1, 'CH3', 1), k['ArStar_CH3_Ar_CH3Star_cm3_3_20'], 'ArStar_CH3_Ar_CH3Star_cm3_3_20')
    push(sto('ArStar', 1, 'H2', 1), sto('Ar', 1, 'H', 2), k['ArStar_H2_Ar_H_H_cm3_3_21'], 'ArStar_H2_Ar_H_H_cm3_3_21')
    push(sto('ArStar', 1, 'C2H2', 1), sto('Ar', 1, 'C2H2Star', 1), k['ArStar_C2H2_Ar_C2H2Star_cm3_3_22'], 'ArStar_C2H2_Ar_C2H2Star_cm3_3_22')
    push(sto('ArStar', 1, 'CH3', 1), sto('Ar', 1, 'CH2', 1, 'H', 1), k['ArStar_CH3_Ar_CH2_H_cm3_3_23'], 'ArStar_CH3_Ar_CH2_H_cm3_3_23')
    push(sto('ArStar', 1, 'C2H4', 1), sto('Ar', 1, 'C2H4', 1), k['ArStar_C2H4_Ar_C2H4Star_cm3_3_24'], 'ArStar_C2H4_Ar_C2H4Star_cm3_3_24')

    # Group 4: Penning Ionization
    push(sto('ArStar', 1, 'CH4', 1), sto('Ar', 1, 'CH4Plus', 1, 'e', 1), k['ArStar_CH4_CH4Plus_cm3_4_1'], 'ArStar_CH4_CH4Plus_cm3_4_1')
    push(sto('ArStar', 1, 'CH4', 1), sto('Ar', 1, 'CH3Plus', 1, 'H', 1, 'e', 1), k['ArStar_CH4_CH3Plus_H_cm3_4_2'], 'ArStar_CH4_CH3Plus_H_cm3_4_2')
    push(sto('ArStar', 1, 'Ar', 1), sto('ArPlus', 1, 'Ar', 1, 'e', 1), k['ArStar_Ar_ArPlus_cm3_4_3'], 'ArStar_Ar_ArPlus_cm3_4_3')
    push(sto('ArStar', 1, 'CH4', 1), sto('ArPlus', 1, 'CH3', 1, 'H', 1, 'e', 1), k['ArStar_CH4_ArPlus_CH3_H_e_cm3_4_4'], 'ArStar_CH4_ArPlus_CH3_H_e_cm3_4_4')
    push(sto('ArStar', 1, 'CH3', 1), sto('ArPlus', 1, 'CH2', 1, 'H', 1, 'e', 1), k['ArStar_CH3_ArPlus_CH2_H_e_cm3_4_5'], 'ArStar_CH3_ArPlus_CH2_H_e_cm3_4_5')
    push(sto('ArStar', 1, 'H', 1), sto('ArPlus', 1, 'HMinus', 1, 'e', 1), k['ArStar_H_ArPlus_HMinus_cm3_4_6'], 'ArStar_H_ArPlus_HMinus_cm3_4_6')
    push(sto('ArStar', 1, 'CH2', 1), sto('ArPlus', 1, 'CH', 1, 'H', 1, 'e', 1), k['ArStar_CH2_ArPlus_CH_H_e_cm3_4_7'], 'ArStar_CH2_ArPlus_CH_H_e_cm3_4_7')
    push(sto('ArStar', 1, 'H2', 1), sto('ArPlus', 1, 'H2', 1, 'e', 1), k['ArStar_H2_ArPlus_H2_e_cm3_4_8'], 'ArStar_H2_ArPlus_H2_e_cm3_4_8')
    push(sto('ArStar', 1, 'C2H2', 1), sto('ArPlus', 1, 'C2H2', 1, 'e', 1), k['ArStar_C2H2_ArPlus_C2H2_e_cm3_4_9'], 'ArStar_C2H2_ArPlus_C2H2_e_cm3_4_9')
    push(sto('ArStar', 1, 'C2H5', 1), sto('ArPlus', 1, 'C2H5', 1, 'e', 1), k['ArStar_C2H5_ArPlus_C2H5_e_cm3_4_10'], 'ArStar_C2H5_ArPlus_C2H5_e_cm3_4_10')
    push(sto('ArStar', 1, 'H', 1), sto('ArPlus', 1, 'H', 1, 'e', 1), k['ArStar_H_ArPlus_H_e_cm3_4_11'], 'ArStar_H_ArPlus_H_e_cm3_4_11')
    push(sto('ArStar', 1, 'CH4', 1), sto('ArPlus', 1, 'CH4', 1, 'e', 1), k['ArStar_CH4_ArPlus_CH4_e_cm3_4_12'], 'ArStar_CH4_ArPlus_CH4_e_cm3_4_12')
    push(sto('ArStar', 1, 'CH3', 1), sto('ArPlus', 1, 'CH3', 1, 'e', 1), k['ArStar_CH3_ArPlus_CH3_e_cm3_4_13'], 'ArStar_CH3_ArPlus_CH3_e_cm3_4_13')
    push(sto('ArStar', 1, 'C2H4', 1), sto('ArPlus', 1, 'C2H4', 1, 'e', 1), k['ArStar_C2H4_ArPlus_C2H4_e_cm3_4_14'], 'ArStar_C2H4_ArPlus_C2H4_e_cm3_4_14')
    push(sto('ArStar', 1, 'C2H5', 1), sto('ArPlus', 1, 'C2H5', 1, 'e', 1), k['ArStar_C2H5_ArPlus_C2H5_e_cm3_4_15'], 'ArStar_C2H5_ArPlus_C2H5_e_cm3_4_15')
    push(sto('ArStar', 1, 'CH2', 1), sto('ArPlus', 1, 'CH2', 1, 'e', 1), k['ArStar_CH2_ArPlus_CH2_e_cm3_4_16'], 'ArStar_CH2_ArPlus_CH2_e_cm3_4_16')
    push(sto('ArStar', 1, 'C2H6', 1), sto('ArPlus', 1, 'C2H6', 1, 'e', 1), k['ArStar_C2H6_ArPlus_C2H6_e_cm3_4_17'], 'ArStar_C2H6_ArPlus_C2H6_e_cm3_4_17')

    # Group 5: Ion-Neutral Reactions
    push(sto('ArPlus', 1, 'CH4', 1), sto('CH3Plus', 1, 'H', 1, 'Ar', 1), k['ArPlus_CH4_CH3Plus_H_cm3_5_1'], 'ArPlus_CH4_CH3Plus_H_cm3_5_1')
    push(sto('CH3Plus', 1, 'CH4', 1), sto('CH5Plus', 1, 'CH2', 1), k['CH3Plus_CH4_CH5Plus_CH2_cm3_5_2'], 'CH3Plus_CH4_CH5Plus_CH2_cm3_5_2')
    push(sto('ArPlus', 1, 'CH3', 1), sto('CH3Plus', 1, 'Ar', 1), k['ArPlus_CH3_CH3Plus_cm3_5_3'], 'ArPlus_CH3_CH3Plus_cm3_5_3')
    push(sto('CH', 2), sto('C2', 1, 'H2', 1), k['CH_CH_C2_H2_cm3_5_4'], 'CH_CH_C2_H2_cm3_5_4')
    push(sto('ArPlus', 1, 'CH4', 1), sto('Ar', 1, 'CH4Plus', 1), k['ArPlus_CH4_Ar_CH4Plus_cm3_5_5'], 'ArPlus_CH4_Ar_CH4Plus_cm3_5_5')
    push(sto('CH4Plus', 1, 'H2', 1), sto('CH5Plus', 1, 'H', 1), k['CH4Plus_H2_CH5Plus_H_cm3_5_6'], 'CH4Plus_H2_CH5Plus_H_cm3_5_6')
    push(sto('ArPlus', 1, 'H2', 1), sto('ArHPlus', 1, 'H', 1), k['ArPlus_H2_ArHPlus_H_cm3_5_7'], 'ArPlus_H2_ArHPlus_H_cm3_5_7')
    push(sto('ArHPlus', 1, 'CH4', 1), sto('Ar', 1, 'CH5Plus', 1), k['ArHPlus_CH4_Ar_CH5Plus_cm3_5_8'], 'ArHPlus_CH4_Ar_CH5Plus_cm3_5_8')
    push(sto('ArPlus', 1, 'CH4', 1), sto('ArHPlus', 1, 'CH3', 1), k['ArPlus_CH4_ArHPlus_CH3_cm3_5_9'], 'ArPlus_CH4_ArHPlus_CH3_cm3_5_9')
    push(sto('CH2', 1, 'CH3Plus', 1), sto('CH3', 1, 'CH2Plus', 1), k['CH2_CH3Plus_CH3_CH2Plus_cm3_5_10'], 'CH2_CH3Plus_CH3_CH2Plus_cm3_5_10')
    push(sto('CH3Plus', 1, 'CH4', 1), sto('C2H5Plus', 1, 'H2', 1), k['CH3Plus_CH4_C2H5Plus_H2_cm3_5_11'], 'CH3Plus_CH4_C2H5Plus_H2_cm3_5_11')
    push(sto('CH5Plus', 1, 'C2H4', 1), sto('C2H5Plus', 1, 'CH4', 1), k['CH5Plus_C2H4_C2H5Plus_CH4_cm3_5_12'], 'CH5Plus_C2H4_C2H5Plus_CH4_cm3_5_12')
    # H3+ chemistry
    push(sto('H2Plus', 1, 'H2', 1), sto('H3Plus', 1, 'H', 1), k['H2Plus_H2_H3Plus_H_cm3_5_13'], 'H2Plus_H2_H3Plus_H_cm3_5_13')
    push(sto('H3Plus', 1, 'CH4', 1), sto('CH5Plus', 1, 'H2', 1), k['H3Plus_CH4_CH5Plus_H2_cm3_5_14'], 'H3Plus_CH4_CH5Plus_H2_cm3_5_14')
    push(sto('H3Plus', 1, 'H2', 1), sto('H2Plus', 1, 'H2', 1), k['H3Plus_H2_H2Plus_H2_cm3_5_15'], 'H3Plus_H2_H2Plus_H2_cm3_5_15')
    # ADDED: CH + Ar⁺ → CHPlus + Ar (fast ion-neutral charge transfer)
    if 'ArPlus_CH_CHPlus_Ar_cm3_5_16' in k:
        push(sto('ArPlus', 1, 'CH', 1), sto('CHPlus', 1, 'Ar', 1), k['ArPlus_CH_CHPlus_Ar_cm3_5_16'], 'ArPlus_CH_CHPlus_Ar_cm3_5_16')

    # Group 6: Dissociative Recombination
    push(sto('ArPlus', 1, 'e', 1), sto('Ar', 1), k['ArPlus_e_Ar_cm3_6_1'], 'ArPlus_e_Ar_cm3_6_1')
    push(sto('CH3Plus', 1, 'e', 1), sto('CH3', 1), k['CH3Plus_e_CH3_cm3_6_2'], 'CH3Plus_e_CH3_cm3_6_2')
    push(sto('CH5Plus', 1, 'e', 1), sto('CH4', 1, 'H', 1), k['CH5Plus_e_CH4_H_cm3_6_3'], 'CH5Plus_e_CH4_H_cm3_6_3')
    push(sto('e', 1, 'CH4Plus', 1), sto('CH3', 1, 'H', 1), k['e_CH4Plus_CH3_H_cm3_6_4'], 'e_CH4Plus_CH3_H_cm3_6_4')
    push(sto('CH3Minus', 1, 'ArPlus', 1), sto('CH3', 1, 'Ar', 1), k['CH3Minus_ArPlus_CH3_Ar_cm3_6_5'], 'CH3Minus_ArPlus_CH3_Ar_cm3_6_5')
    push(sto('CH3Minus', 1, 'CH4Plus', 1), sto('CH4', 1, 'CH3', 1), k['CH3Minus_CH4Plus_CH4_CH3_cm3_6_6'], 'CH3Minus_CH4Plus_CH4_CH3_cm3_6_6')
    push(sto('CH3Minus', 1, 'CH3Plus', 1), sto('CH4', 1, 'CH2', 1), k['CH3Minus_CH3Plus_CH4_CH2_cm3_6_7'], 'CH3Minus_CH3Plus_CH4_CH2_cm3_6_7')
    push(sto('CH5Plus', 1, 'e', 1), sto('CH3', 1, 'H2', 1), k['CH5Plus_e_CH3_H2_cm3_6_8'], 'CH5Plus_e_CH3_H2_cm3_6_8')
    push(sto('CH4Plus', 1, 'e', 1), sto('CH2', 1, 'H2', 1), k['e_CH4Plus_CH2_H2_cm3_6_9'], 'e_CH4Plus_CH2_H2_cm3_6_9')
    push(sto('CH5Plus', 1, 'e', 1), sto('CH2', 1, 'H2', 1, 'H', 1), k['CH5Plus_e_CH2_H2_H_cm3_6_10'], 'CH5Plus_e_CH2_H2_H_cm3_6_10')
    push(sto('CH4Plus', 1, 'e', 1), sto('CH', 1, 'H2', 1, 'H', 1), k['e_CH4Plus_CH_H2_H_cm3_6_11'], 'e_CH4Plus_CH_H2_H_cm3_6_11')
    push(sto('CH5Plus', 1, 'e', 1), sto('CH3', 1, 'H', 2), k['CH5Plus_e_CH3_2H_cm3_6_12'], 'CH5Plus_e_CH3_2H_cm3_6_12')
    push(sto('CH4Plus', 1, 'e', 1), sto('C', 1, 'H2', 2), k['e_CH4Plus_C_2H2_cm3_6_13'], 'e_CH4Plus_C_2H2_cm3_6_13')
    push(sto('C2H5Plus', 1, 'e', 1), sto('C2H4', 1, 'H', 1), k['C2H5Plus_e_C2H4_H_cm3_6_14'], 'C2H5Plus_e_C2H4_H_cm3_6_14')
    push(sto('C2H4Plus', 1, 'e', 1), sto('C2H2', 1, 'H2', 1), k['C2H4Plus_e_C2H2_H2_cm3_6_15'], 'C2H4Plus_e_C2H2_H2_cm3_6_15')
    push(sto('C2H3Plus', 1, 'e', 1), sto('C2H2', 1, 'H', 1), k['C2H3Plus_e_C2H2_H_cm3_6_16'], 'C2H3Plus_e_C2H2_H_cm3_6_16')
    push(sto('HMinus', 1, 'ArPlus', 1), sto('H', 1, 'Ar', 1), k['HMinus_ArPlus_H_Ar_cm3_6_17'], 'HMinus_ArPlus_H_Ar_cm3_6_17')
    push(sto('C2HPlus', 1, 'e', 1), sto('C2', 1, 'H', 1), k['C2HPlus_e_C2_H_cm3_6_18'], 'C2HPlus_e_C2_H_cm3_6_18')
    push(sto('HMinus', 1, 'CH5Plus', 1), sto('CH4', 1, 'H2', 1, 'H', 1), k['HMinus_CH5Plus_CH4_H2_H_cm3_6_19'], 'HMinus_CH5Plus_CH4_H2_H_cm3_6_19')
    push(sto('CH4Plus', 1, 'HMinus', 1), sto('CH4', 1, 'H', 1), k['CH4Plus_HMinus_CH4_H_cm3_6_20'], 'CH4Plus_HMinus_CH4_H_cm3_6_20')
    push(sto('CH3Plus', 1, 'HMinus', 1), sto('CH4', 1, 'H2', 1), k['CH3Plus_HMinus_CH4_H2_cm3_6_21'], 'CH3Plus_HMinus_CH4_H2_cm3_6_21')
    push(sto('C2H5Plus', 1, 'HMinus', 1), sto('C2H6', 1, 'H', 1), k['C2H5Plus_HMinus_C2H6_H_cm3_6_22'], 'C2H5Plus_HMinus_C2H6_H_cm3_6_22')
    push(sto('ArHPlus', 1, 'HMinus', 1), sto('Ar', 1, 'H2', 1, 'H', 1), k['ArHPlus_HMinus_Ar_H2_H_cm3_6_23'], 'ArHPlus_HMinus_Ar_H2_H_cm3_6_23')
    push(sto('CH5Plus', 1, 'CH3Minus', 1), sto('CH4', 2, 'H', 1), k['CH5Plus_CH3Minus_CH4_CH4_H_cm3_6_24'], 'CH5Plus_CH3Minus_CH4_CH4_H_cm3_6_24')
    push(sto('CH4Plus', 1, 'CH3Minus', 1), sto('CH4', 1, 'CH3', 1, 'H', 1), k['CH4Plus_CH3Minus_CH4_CH3_H_cm3_6_25'], 'CH4Plus_CH3Minus_CH4_CH3_H_cm3_6_25')
    push(sto('CH3Plus', 1, 'CH3Minus', 1), sto('CH4', 1, 'CH2', 1, 'H', 1), k['CH3Plus_CH3Minus_CH4_CH2_H_cm3_6_26'], 'CH3Plus_CH3Minus_CH4_CH2_H_cm3_6_26')
    push(sto('C2H5Plus', 1, 'CH3Minus', 1), sto('C2H6', 1, 'H', 1), k['C2H5Plus_CH3Minus_C2H6_H_cm3_6_27'], 'C2H5Plus_CH3Minus_C2H6_H_cm3_6_27')
    push(sto('C2H5Plus', 1, 'e', 1), sto('C2H4', 1, 'H', 1), k['C2H5Plus_e_C2H4_H_cm3_6_28'], 'C2H5Plus_e_C2H4_H_cm3_6_28')
    # H2+ and H3+ recombination
    push(sto('H2Plus', 1, 'e', 1), sto('H', 2), k['H2Plus_e_H_H_cm3_6_29'], 'H2Plus_e_H_H_cm3_6_29')
    push(sto('H3Plus', 1, 'e', 1), sto('H2', 1, 'H', 1), k['H3Plus_e_H2_H_cm3_6_30'], 'H3Plus_e_H2_H_cm3_6_30')
    push(sto('H3Plus', 1, 'e', 1), sto('H', 3), k['H3Plus_e_H_H_H_cm3_6_31'], 'H3Plus_e_H_H_H_cm3_6_31')
    push(sto('CHPlus', 1, 'e', 1), sto('CH', 1), k['CHPlus_e_CH_cm3_6_32'], 'CHPlus_e_CH_cm3_6_32')  # CHPlus recombination (NEW!)

    # NEW IN V7: ArH+ dissociative recombination
    push(sto('ArHPlus', 1, 'e', 1), sto('Ar', 1, 'H', 1), k['ArHPlus_e_Ar_H_cm3_6_29'], 'ArHPlus_e_Ar_H_cm3_6_29')

    # Group 8: Attachment Reactions
    # NEW IN V7: Dissociative attachment
    push(sto('e', 1, 'CH4', 1), sto('CH3', 1, 'HMinus', 1), k['e_CH4_CH3_HMinus_cm3_8_1'], 'e_CH4_CH3_HMinus_cm3_8_1')

    # Group 7: Neutral-Neutral Reactions (continuing...)
    push(sto('CH2', 1, 'H', 1), sto('CH', 1, 'H2', 1), k['CH2_H_CH_H2_cm3_7_1'], 'CH2_H_CH_H2_cm3_7_1')
    push(sto('CH2', 1, 'H', 1), sto('C', 1, 'H2', 1, 'H', 1), k['CH2_H_C_H2_H_cm3_7_2'], 'CH2_H_C_H2_H_cm3_7_2')
    push(sto('CH', 1, 'H', 1), sto('C', 1, 'H2', 1), k['CH_H_C_H2_cm3_7_3'], 'CH_H_C_H2_cm3_7_3')
    push(sto('C', 1, 'CH', 1), sto('C2', 1, 'H', 1), k['C_CH_C2_H_cm3_7_4'], 'C_CH_C2_H_cm3_7_4')
    push(sto('CH', 1, 'CH3', 1), sto('C2H4', 1), k['CH_CH3_C2H4_cm3_7_5'], 'CH_CH3_C2H4_cm3_7_5')
    push(sto('C2', 1, 'H', 1), sto('CH', 1, 'C', 1), k['C2_H_CH_C_cm3_7_6'], 'C2_H_CH_C_cm3_7_6')
    push(sto('CH', 1, 'CH2', 1), sto('C2H2', 1, 'H', 1), k['CH_CH2_C2H2_H_cm3_7_7'], 'CH_CH2_C2H2_H_cm3_7_7')
    push(sto('C', 1, 'CH3', 1), sto('C2H2', 1, 'H', 1), k['C_CH3_C2_H2_H_cm3_7_8'], 'C_CH3_C2_H2_H_cm3_7_8')
    push(sto('CH', 1, 'C', 1), sto('C2', 1, 'H', 1), k['CH_C_C2_H_cm3_7_9'], 'CH_C_C2_H_cm3_7_9')
    push(sto('CH', 1, 'CH3', 1), sto('C2H3', 1, 'H', 1), k['CH_CH3_C2H3_H_cm3_7_10'], 'CH_CH3_C2H3_H_cm3_7_10')
    push(sto('CH', 1, 'Ar', 1), sto('Ar', 1, 'C', 1, 'H', 1), k['CH_Ar_Ar_C_H_cm3_7_11'], 'CH_Ar_Ar_C_H_cm3_7_11')
    push(sto('C', 1, 'H', 1), sto('CH', 1), k['C_H_CH_cm3_7_12'], 'C_H_CH_cm3_7_12')
    push(sto('CH2', 2), sto('C2H4', 1), k['CH2_CH2_C2H4_cm3_7_13'], 'CH2_CH2_C2H4_cm3_7_13')
    push(sto('CH3', 1, 'CH2', 1), sto('C2H5', 1), k['CH3_CH2_C2H5_cm3_7_14'], 'CH3_CH2_C2H5_cm3_7_14')
    push(sto('CH2', 2), sto('C2H2', 1, 'H2', 1), k['CH2_CH2_C2H2_H2_cm3_7_15'], 'CH2_CH2_C2H2_H2_cm3_7_15')
    push(sto('CH3', 1, 'CH', 1), sto('C2H2', 1, 'H2', 1), k['CH3_CH_C2H2_H2_cm3_7_16'], 'CH3_CH_C2H2_H2_cm3_7_16')
    push(sto('CH2', 1, 'C', 1), sto('C2H2', 1), k['CH2_C_C2H2_cm3_7_17'], 'CH2_C_C2H2_cm3_7_17')
    push(sto('CH', 1, 'C2H4', 1), sto('C2H2', 1, 'CH3', 1), k['CH_C2H4_C2H2_CH3_cm3_7_18'], 'CH_C2H4_C2H2_CH3_cm3_7_18')
    push(sto('C2H2', 1, 'C', 1), sto('C2', 1, 'CH2', 1), k['C2H2_C_C2_CH2_cm3_7_19'], 'C2H2_C_C2_CH2_cm3_7_19')
    push(sto('CH', 1, 'CH4', 1), sto('C2H4', 1, 'H', 1), k['CH_CH4_C2H4_H_cm3_7_20'], 'CH_CH4_C2H4_H_cm3_7_20')
    push(sto('CH', 1, 'H', 1), sto('CH2', 1), k['CH_H_CH2_cm3_7_21'], 'CH_H_CH2_cm3_7_21')
    push(sto('CH', 1, 'C2H2', 1), sto('C3H2', 1, 'H', 1), k['CH_C2H2_C3H2_H_cm3_7_22'], 'CH_C2H2_C3H2_H_cm3_7_22')
    push(sto('CH', 1, 'CH3', 1), sto('C2H2', 1, 'H2', 1), k['CH_CH3_C2H2_H2_cm3_7_23'], 'CH_CH3_C2H2_H2_cm3_7_23')
    push(sto('CH', 1, 'C', 1), sto('C2', 1, 'H2', 1), k['CH_C_C2_H2_cm3_7_24'], 'CH_C_C2_H2_cm3_7_24')
    push(sto('H', 1, 'CH4', 1), sto('CH3', 1, 'H2', 1), k['H_CH4_CH3_H2_cm3_7_25'], 'H_CH4_CH3_H2_cm3_7_25')
    push(sto('CH2', 1, 'CH', 1), sto('C2', 1, 'H2', 1, 'H', 1), k['CH2_CH_C2_H2_H_cm3_7_26'], 'CH2_CH_C2_H2_H_cm3_7_26')
    push(sto('CH', 1, 'C2H2', 1), sto('C3H', 1, 'H2', 1), k['CH_C2H2_C3H_H2_cm3_7_27'], 'CH_C2H2_C3H_H2_cm3_7_27')
    push(sto('CH', 1, 'C3H', 1), sto('C4H2', 1, 'H', 1), k['CH_C3H_C4H2_H_cm3_7_28'], 'CH_C3H_C4H2_H_cm3_7_28')
    push(sto('CH', 1, 'C2H2', 1), sto('C2H', 1, 'CH2', 1), k['CH_C2H2_C2H_CH2_cm3_7_29'], 'CH_C2H2_C2H_CH2_cm3_7_29')
    push(sto('CH', 1, 'H2', 1), sto('CH2', 1, 'H', 1), k['CH_H2_CH2_H_cm3_7_30'], 'CH_H2_CH2_H_cm3_7_30')
    push(sto('CH', 1, 'C2H3', 1), sto('C3H3', 1, 'H', 1), k['CH_C2H3_C3H3_H_cm3_7_31'], 'CH_C2H3_C3H3_H_cm3_7_31')
    push(sto('CH', 1, 'C2H4', 1), sto('C3H4', 1, 'H', 1), k['CH_C2H4_C3H4_H_cm3_7_32'], 'CH_C2H4_C3H4_H_cm3_7_32')
    push(sto('CH', 1, 'C2', 1), sto('C3', 1, 'H', 1), k['CH_C2_C3_H_cm3_7_33'], 'CH_C2_C3_H_cm3_7_33')
    push(sto('CH', 1, 'C2H5', 1), sto('C3H5', 1, 'H', 1), k['CH_C2H5_C3H5_H_cm3_7_34'], 'CH_C2H5_C3H5_H_cm3_7_34')
    push(sto('CH', 1, 'C3H2', 1), sto('C4H', 1, 'H2', 1), k['CH_C3H2_C4H_H2_cm3_7_35'], 'CH_C3H2_C4H_H2_cm3_7_35')
    push(sto('CH3', 1, 'H', 1), sto('CH2', 1, 'H2', 1), k['CH3_H_CH2_H2_cm3_7_36'], 'CH3_H_CH2_H2_cm3_7_36')
    push(sto('CH', 1, 'C2H6', 1), sto('C3H6', 1, 'H', 1), k['CH_C2H6_C3H6_H_cm3_7_37'], 'CH_C2H6_C3H6_H_cm3_7_37')
    push(sto('CH3', 2), sto('CH2', 1, 'CH4', 1), k['CH3_CH3_CH2_CH4_cm3_7_38'], 'CH3_CH3_CH2_CH4_cm3_7_38')
    push(sto('CH', 1, 'CH4', 1), sto('CH2', 1, 'CH3', 1), k['CH_CH4_CH2_CH3_cm3_7_39'], 'CH_CH4_CH2_CH3_cm3_7_39')
    push(sto('CH3', 2), sto('C2H6', 1), k['CH3_CH3_C2H6_cm3_7_40'], 'CH3_CH3_C2H6_cm3_7_40')
    push(sto('CH', 1, 'C2H5', 1), sto('C3H6', 1), k['CH_C2H5_C3H6_cm3_7_41'], 'CH_C2H5_C3H6_cm3_7_41')
    push(sto('CH2', 2), sto('C2H2', 1, 'H2', 1), k['CH2_CH2_C2H2_H2_cm3_7_42'], 'CH2_CH2_C2H2_H2_cm3_7_42')
    push(sto('CH', 1, 'C', 1), sto('C2', 1, 'H', 1), k['CH_C_C2_H_cm3_7_43'], 'CH_C_C2_H_cm3_7_43')
    push(sto('CH', 2), sto('C2', 1, 'H2', 1), k['CH_CH_C2_H2_cm3_7_44'], 'CH_CH_C2_H2_cm3_7_44')
    push(sto('CH', 1, 'C2H6', 1), sto('C3H6', 1, 'H', 1), k['CH_C2H6_C3H6_H_cm3_7_45'], 'CH_C2H6_C3H6_H_cm3_7_45')
    push(sto('CH', 1, 'C2H4', 1), sto('C2H2', 1, 'CH3', 1), k['CH_C2H4_C2H2_CH3_cm3_7_46'], 'CH_C2H4_C2H2_CH3_cm3_7_46')
    push(sto('C2H', 1, 'H', 1), sto('C2', 1, 'H2', 1), k['C2H_H_C2_H2_cm3_7_47'], 'C2H_H_C2_H2_cm3_7_47')
    push(sto('CH', 1, 'CH2', 1), sto('C2H2', 1, 'H', 1), k['CH_CH2_C2H2_H_cm3_7_48'], 'CH_CH2_C2H2_H_cm3_7_48')
    push(sto('CH3', 2), sto('C2H2', 1, 'H2', 2), k['CH3_CH3_C2H2_H2_H2_cm3_7_49'], 'CH3_CH3_C2H2_H2_H2_cm3_7_49')
    push(sto('C2H2', 1, 'H', 1), sto('C2', 1, 'H2', 1, 'H', 1), k['C2H2_H_C2_H2_H_cm3_7_50'], 'C2H2_H_C2_H2_H_cm3_7_50')
    push(sto('CH', 1, 'H2', 1), sto('CH2', 1, 'H', 1), k['CH_H2_CH2_H_cm3_7_51'], 'CH_H2_CH2_H_cm3_7_51')
    push(sto('C2', 1, 'CH', 1), sto('C3', 1, 'H', 1), k['C2_CH_C3_H_cm3_7_52'], 'C2_CH_C3_H_cm3_7_52')
    push(sto('CH2', 1, 'C2H3', 1), sto('C2H2', 1, 'CH3', 1), k['CH2_C2H3_C2H2_CH3_cm3_7_53'], 'CH2_C2H3_C2H2_CH3_cm3_7_53')
    push(sto('C', 1, 'C2H3', 1), sto('C2', 1, 'CH3', 1), k['C_C2H3_C2_CH3_cm3_7_54'], 'C_C2H3_C2_CH3_cm3_7_54')
    push(sto('CH', 1, 'C2H5', 1), sto('C3H6', 1), k['CH_C2H5_C3H6_cm3_7_55'], 'CH_C2H5_C3H6_cm3_7_55')
    push(sto('C2H2', 1, 'CH', 1), sto('C3', 1, 'H2', 1, 'H', 1), k['C2H2_CH_C3_H2_H_cm3_7_56'], 'C2H2_CH_C3_H2_H_cm3_7_56')
    push(sto('CH', 1, 'C2H6', 1), sto('C2H2', 1, 'CH3', 1, 'H', 1), k['CH_C2H6_C2H2_CH3_H_cm3_7_57'], 'CH_C2H6_C2H2_CH3_H_cm3_7_57')
    push(sto('CH2', 2), sto('C2', 1, 'H2', 2), k['CH2_CH2_C2_H2_H2_cm3_7_58'], 'CH2_CH2_C2_H2_H2_cm3_7_58')
    push(sto('CH', 1, 'CH4', 1), sto('C2H4', 1, 'H', 1), k['CH_CH4_C2H4_H_cm3_7_59'], 'CH_CH4_C2H4_H_cm3_7_59')
    push(sto('CH', 1, 'C2H4', 1), sto('C3H4', 1, 'H', 1), k['CH_C2H4_C3H4_H_cm3_7_60'], 'CH_C2H4_C3H4_H_cm3_7_60')
    push(sto('CH3', 1, 'C2H5', 1), sto('C2H2', 1, 'CH3', 1, 'H2', 1), k['CH3_C2H5_C2H2_CH3_H2_cm3_7_61'], 'CH3_C2H5_C2H2_CH3_H2_cm3_7_61')
    push(sto('CH2', 1, 'CH3', 1), sto('C2H2', 1, 'H', 1, 'H2', 1), k['CH2_CH3_C2H2_H_H2_cm3_7_62'], 'CH2_CH3_C2H2_H_H2_cm3_7_62')
    push(sto('CH2', 1, 'C2H5', 1), sto('C2H2', 1, 'CH3', 1, 'H', 1), k['CH2_C2H5_C2H2_CH3_H_cm3_7_63'], 'CH2_C2H5_C2H2_CH3_H_cm3_7_63')
    # New reactions from audit
    push(sto('C', 2), sto('C2', 1), k['C_C_M_C2_M_cm6_7_64'], 'C_C_M_C2_M_cm6_7_64')
    push(sto('H', 1, 'C2H4', 1), sto('C2H3', 1, 'H2', 1), k['H_C2H4_C2H3_H2_cm3_7_65'], 'H_C2H4_C2H3_H2_cm3_7_65')

    # MISSING CH3 PRODUCTION PATHWAYS (added for physical realism)
    push(sto('e', 1, 'C2H4', 1), sto('CH3', 1, 'CH', 1), k['e_C2H4_CH3_CH_cm3_7_66'], 'e_C2H4_CH3_CH_cm3_7_66')
    push(sto('e', 1, 'C2H6', 1), sto('CH3', 2), k['e_C2H6_CH3_CH3_cm3_7_67'], 'e_C2H6_CH3_CH3_cm3_7_67')
    push(sto('e', 1, 'C2H5', 1), sto('CH3', 1, 'CH2', 1), k['e_C2H5_CH3_CH2_cm3_7_68'], 'e_C2H5_CH3_CH2_cm3_7_68')
    push(sto('ArStar', 1, 'C2H4', 1), sto('Ar', 1, 'CH3', 1, 'CH', 1), k['ArStar_C2H4_CH3_CH_cm3_7_69'], 'ArStar_C2H4_CH3_CH_cm3_7_69')
    push(sto('ArStar', 1, 'C2H6', 1), sto('Ar', 1, 'CH3', 2), k['ArStar_C2H6_CH3_CH3_cm3_7_70'], 'ArStar_C2H6_CH3_CH3_cm3_7_70')
    push(sto('H', 1, 'C2H5', 1), sto('CH3', 1, 'CH2', 1), k['H_C2H5_CH3_CH2_cm3_7_71'], 'H_C2H5_CH3_CH2_cm3_7_71')
    push(sto('CH2', 2), sto('CH3', 1, 'CH', 1), k['CH2_CH2_CH3_CH_cm3_7_72'], 'CH2_CH2_CH3_CH_cm3_7_72')
    push(sto('C2H5Plus', 1, 'e', 1), sto('CH3', 1, 'CH2', 1), k['C2H5Plus_e_CH3_CH2_cm3_7_73'], 'C2H5Plus_e_CH3_CH2_cm3_7_73')
    push(sto('CH2', 1, 'H', 1), sto('CH3', 1), k['CH2_H_M_CH3_M_cm6_7_74'], 'CH2_H_M_CH3_M_cm6_7_74')
    push(sto('ArPlus', 1, 'CH4', 1), sto('CH3Plus', 1, 'Ar', 1, 'H', 1), k['ArPlus_CH4_CH3Plus_ArH_cm3_7_75'], 'ArPlus_CH4_CH3Plus_ArH_cm3_7_75')

    # Group 8: Termolecular Recombination
    push(sto('H', 2), sto('H2', 1), k['H_H_M_H2_M_cm6_8_1'], 'H_H_M_H2_M_cm6_8_1')
    push(sto('CH3', 2), sto('C2H6', 1), k['CH3_CH3_M_C2H6_M_cm6_8_2'], 'CH3_CH3_M_C2H6_M_cm6_8_2')
    push(sto('CH3', 1, 'H', 1), sto('CH4', 1), k['CH3_H_M_CH4_M_cm6_8_3'], 'CH3_H_M_CH4_M_cm6_8_3')

    # Three-body electron-ion recombination (CRITICAL for high-density stabilization)
    push(sto('e', 1, 'ArPlus', 1), sto('Ar', 1), k['e_ArPlus_M_Ar_M_cm6_8_4'], 'e_ArPlus_M_Ar_M_cm6_8_4')
    push(sto('e', 1, 'CH4Plus', 1), sto('CH4', 1), k['e_CH4Plus_M_CH4_M_cm6_8_5'], 'e_CH4Plus_M_CH4_M_cm6_8_5')
    push(sto('e', 1, 'CH3Plus', 1), sto('CH3', 1), k['e_CH3Plus_M_CH3_M_cm6_8_6'], 'e_CH3Plus_M_CH3_M_cm6_8_6')
    push(sto('e', 1, 'CH5Plus', 1), sto('CH4', 1, 'H', 1), k['e_CH5Plus_M_CH5_M_cm6_8_7'], 'e_CH5Plus_M_CH5_M_cm6_8_7')  # CH5 → CH4 + H
    push(sto('e', 1, 'ArHPlus', 1), sto('Ar', 1, 'H', 1), k['e_ArHPlus_M_ArH_M_cm6_8_8'], 'e_ArHPlus_M_ArH_M_cm6_8_8')  # ArH → Ar + H
    push(sto('e', 1, 'C2H5Plus', 1), sto('C2H4', 1, 'H', 1), k['e_C2H5Plus_M_C2H5_M_cm6_8_9'], 'e_C2H5Plus_M_C2H5_M_cm6_8_9')  # C2H5 → C2H4 + H

    # Group 9: Stick Reactions
    push(sto('H', 1), sto(), k['stick_H_9_1'], 'stick_H_9_1')
    push(sto('CH3', 1), sto(), k['stick_CH3_9_2'], 'stick_CH3_9_2')
    push(sto('CH', 1), sto(), k['stick_CH_9_3'], 'stick_CH_9_3')
    push(sto('ArPlus', 1), sto(), k['stick_ArPlus_9_4'], 'stick_ArPlus_9_4')
    push(sto('ArStar', 1), sto(), k['stick_ArStar_9_5'], 'stick_ArStar_9_5')
    push(sto('CH3Plus', 1), sto(), k['stick_CH3Plus_9_6'], 'stick_CH3Plus_9_6')
    push(sto('CH5Plus', 1), sto(), k['stick_CH5Plus_9_7'], 'stick_CH5Plus_9_7')
    push(sto('ArHPlus', 1), sto(), k['stick_ArHPlus_9_8'], 'stick_ArHPlus_9_8')
    push(sto('C2', 1), sto(), k['stick_C2_9_9'], 'stick_C2_9_9')
    push(sto('C', 1), sto(), k['stick_C_9_10'], 'stick_C_9_10')
    push(sto('C2H2', 1), sto(), k['stick_C2H2_9_11'], 'stick_C2H2_9_11')
    push(sto('C2H4', 1), sto(), k['stick_C2H4_9_12'], 'stick_C2H4_9_12')
    push(sto('CH2', 1), sto(), k['stick_CH2_9_13'], 'stick_CH2_9_13')
    push(sto('C2H6', 1), sto(), k['stick_C2H6_9_14'], 'stick_C2H6_9_14')
    push(sto('CH3Minus', 1), sto(), k['stick_CH3Minus_9_15'], 'stick_CH3Minus_9_15')
    push(sto('H2', 1), sto(), k['stick_H2_9_16'], 'stick_H2_9_16')
    push(sto('C2H5', 1), sto(), k['stick_C2H5_9_17'], 'stick_C2H5_9_17')
    push(sto('HMinus', 1), sto(), k['stick_HMinus_9_18'], 'stick_HMinus_9_18')
    push(sto('C3H', 1), sto(), k['stick_C3H_9_19'], 'stick_C3H_9_19')
    push(sto('C4H2', 1), sto(), k['stick_C4H2_9_20'], 'stick_C4H2_9_20')
    push(sto('C3H3', 1), sto(), k['stick_C3H3_9_21'], 'stick_C3H3_9_21')
    push(sto('C3H4', 1), sto(), k['stick_C3H4_9_22'], 'stick_C3H4_9_22')
    push(sto('C3', 1), sto(), k['stick_C3_9_23'], 'stick_C3_9_23')
    push(sto('C4H', 1), sto(), k['stick_C4H_9_24'], 'stick_C4H_9_24')
    push(sto('C3H6', 1), sto(), k['stick_C3H6_9_25'], 'stick_C3H6_9_25')
    push(sto('H3Plus', 1), sto(), k['stick_H3Plus_9_26'], 'stick_H3Plus_9_26')
    push(sto('CHPlus', 1), sto(), k['stick_CHPlus_9_27'], 'stick_CHPlus_9_27')
    push(sto('C2HPlus', 1), sto(), k['stick_C2HPlus_9_28'], 'stick_C2HPlus_9_28')
    push(sto('H2Plus', 1), sto(), k['stick_H2Plus_9_29'], 'stick_H2Plus_9_29')

    # Group 10: Drift Losses
    push(sto('ArPlus', 1), sto(), k['drift_ArPlus_10_1'], 'drift_ArPlus_10_1')
    push(sto('CH4Plus', 1), sto(), k['drift_CH4Plus_10_2'], 'drift_CH4Plus_10_2')
    push(sto('CH3Plus', 1), sto(), k['drift_CH3Plus_10_3'], 'drift_CH3Plus_10_3')
    push(sto('CH5Plus', 1), sto(), k['drift_CH5Plus_10_4'], 'drift_CH5Plus_10_4')
    push(sto('ArHPlus', 1), sto(), k['drift_ArHPlus_10_5'], 'drift_ArHPlus_10_5')
    push(sto('CH2Plus', 1), sto(), k['drift_CH2Plus_10_6'], 'drift_CH2Plus_10_6')
    push(sto('C2H5Plus', 1), sto(), k['drift_C2H5Plus_10_7'], 'drift_C2H5Plus_10_7')
    push(sto('C2H4Plus', 1), sto(), k['drift_C2H4Plus_10_8'], 'drift_C2H4Plus_10_8')
    push(sto('C2H3Plus', 1), sto(), k['drift_C2H3Plus_10_9'], 'drift_C2H3Plus_10_9')
    push(sto('H3Plus', 1), sto(), k['drift_H3Plus_10_10'], 'drift_H3Plus_10_10')
    push(sto('CHPlus', 1), sto(), k['drift_CHPlus_10_11'], 'drift_CHPlus_10_11')
    push(sto('CH3Minus', 1), sto(), k['drift_CH3Minus_10_12'], 'drift_CH3Minus_10_12')
    push(sto('C2HPlus', 1), sto(), k['drift_C2HPlus_10_13'], 'drift_C2HPlus_10_13')
    push(sto('H2Plus', 1), sto(), k['drift_H2Plus_10_14'], 'drift_H2Plus_10_14')

    # Group 11: Loss Reactions
    push(sto('CH2', 1), sto(), k['loss_CH2_11_1'], 'loss_CH2_11_1')
    push(sto('H2', 1), sto(), k['loss_H2_11_2'], 'loss_H2_11_2')
    push(sto('C2', 1), sto(), k['loss_C2_11_3'], 'loss_C2_11_3')
    push(sto('e', 1), sto(), k['loss_e_11_4'], 'loss_e_11_4')
    push(sto('C2H6', 1), sto(), k['loss_C2H6_11_5'], 'loss_C2H6_11_5')
    push(sto('CH4', 1), sto(), k['loss_CH4_11_6'], 'loss_CH4_11_6')
    push(sto('Ar', 1), sto(), k['loss_Ar_11_7'], 'loss_Ar_11_7')
    push(sto('C', 1), sto(), k['loss_C_11_8'], 'loss_C_11_8')
    push(sto('CH', 1), sto(), k['loss_CH_11_9'], 'loss_CH_11_9')
    push(sto('C3H', 1), sto(), k['loss_C3H_11_10'], 'loss_C3H_11_10')
    push(sto('C4H2', 1), sto(), k['loss_C4H2_11_11'], 'loss_C4H2_11_11')
    push(sto('C2H', 1), sto(), k['loss_C2H_11_12'], 'loss_C2H_11_12')
    push(sto('C3H3', 1), sto(), k['loss_C3H3_11_13'], 'loss_C3H3_11_13')
    push(sto('C3', 1), sto(), k['loss_C3_11_14'], 'loss_C3_11_14')
    push(sto('C4H', 1), sto(), k['loss_C4H_11_15'], 'loss_C4H_11_15')
    push(sto('C3H6', 1), sto(), k['loss_C3H6_11_16'], 'loss_C3H6_11_16')
    push(sto('C2H5', 1), sto(), k['loss_C2H5_11_17'], 'loss_C2H5_11_17')
    push(sto('C3H4', 1), sto(), k['loss_C3H4_11_18'], 'loss_C3H4_11_18')
    push(sto('C2H2', 1), sto(), k['loss_C2H2_11_19'], 'loss_C2H2_11_19')
    push(sto('C2H4', 1), sto(), k['loss_C2H4_11_20'], 'loss_C2H4_11_20')
    push(sto('CH3', 1), sto(), k['loss_CH3_11_21'], 'loss_CH3_11_21')
    push(sto('C2H3', 1), sto(), k['loss_C2H3_11_22'], 'loss_C2H3_11_22')
    push(sto('C3H2', 1), sto(), k['loss_C3H2_11_23'], 'loss_C3H2_11_23')
    push(sto('C3H5', 1), sto(), k['loss_C3H5_11_24'], 'loss_C3H5_11_24')
    push(sto('C2H2Star', 1), sto('C2H2', 1), k['loss_C2H2Star_11_25'], 'loss_C2H2Star_11_25')

    # DUST/NANOPARTICLE LOSS REACTIONS (if enabled)
    # Loss to growing dust particles via sticking
    for sp_name in ['CH', 'CH2', 'CH3', 'C', 'C2', 'H']:
        dust_key = f'dust_loss_{sp_name}_12'
        if dust_key in k:
            push(sto(sp_name, 1), sto(), k[dust_key], dust_key)

    return R, tags
