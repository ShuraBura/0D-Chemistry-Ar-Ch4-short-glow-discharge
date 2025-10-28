"""
define_rates_tunable.py - Tunable rate coefficients for CG region
Based on experimental cathode fall: 253 V / 4 mm = 632.5 V/cm
Based on diffusion length: 0-1 mm
"""

import numpy as np


def define_rates_tunable(params):
    """
    Define all reaction rate constants with tunability.

    Key tunable parameters in params:
    - Te: electron temperature (eV) - range: 1-7 eV
    - ne: electron density (cm⁻³) - range: 1e7-1e11 cm⁻³
    - E_field: electric field (V/cm) - range: 600-1000 V/cm
    - Tgas: gas temperature (K) - range: 300-500 K
    - L_diff: diffusion length (cm) - range: 0.05-0.3 cm
    - scale_e_impact: scaling factor for e-impact rates - range: 0.5-2.0
    """
    k = {}

    # Extract parameters
    E_field = params.get('E_field', 800)  # V/cm
    L_discharge = params['L_discharge']
    mobilities = params['mobilities']
    Tgas = params.get('Tgas', params.get('Tg', 400))  # K
    L_diff = params.get('L_diff', 0.1)  # cm (CRITICAL - from experiment!)
    scale_e = params.get('scale_e_impact', 1.0)  # Rate scaling

    # Diffusion coefficients at reference conditions (400 K, 0.4 Torr)
    # Scale with temperature: D ∝ T^1.75 / P
    T_ref = 400  # K
    P = params.get('P', 0.4)  # Torr
    T_scale = (Tgas / T_ref)**1.75
    P_scale = (0.4 / P)
    D_scale = T_scale * P_scale

    # Reference diffusion coefficients (cm²/s at 400 K, 0.4 Torr)
    D_ref = {
        'H': 300, 'CH': 150, 'CH2': 120, 'CH3': 100, 'CH4': 80,
        'C': 120, 'C2': 100, 'C2H': 90, 'C2H2': 70, 'C2H3': 80,
        'C2H4': 70, 'C2H5': 70, 'C2H6': 60, 'H2': 180, 'Ar': 50,
        'e': 1e5,  # Electrons much faster (but ambipolar limited)
        'C3H': 60, 'C3H2': 60, 'C3H3': 50, 'C3H4': 50, 'C3H5': 50,
        'C3H6': 40, 'C3': 80, 'C4H': 50, 'C4H2': 50,
    }

    # ===================================================================
    # Group 1: Electron-Impact Reactions (Neutral Products)
    # Updated from Janev-Reiter 2002, with scaling factor
    # ===================================================================
    k['e_CH4_CH3_H_cm3_1_1'] = 4.2e-11 * scale_e
    k['e_CH4_CH2_H2_cm3_1_2'] = 1.1e-11 * scale_e
    k['e_CH4_CH_H2_H_vib_cm3_1_3'] = 0.7e-11 * scale_e
    k['e_H2_H_H_cm3_1_4'] = 6e-12 * scale_e
    k['e_CH3_CH2_H_cm3_1_5'] = 3e-11 * scale_e
    k['e_C2H4_C2H2_H2_cm3_1_6'] = 1.0e-11 * scale_e
    k['e_Ar_ArStar_cm3_1_7'] = 6e-11 * scale_e
    k['e_C2H6_C2H4_H2_cm3_1_8'] = 7.0e-12 * scale_e
    k['e_C2H6_C2H4_H2_e_cm3_1_9'] = 1.2e-11 * scale_e
    k['e_CH4_CH3Minus_H_cm3_1_10'] = 6e-18 * scale_e
    k['e_CH4_CH_H_H2_cm3_1_11'] = 2e-11 * scale_e
    k['e_CH_CH_C_H_e_cm3_1_12'] = 6e-11 * scale_e
    k['e_H2_HMinus_H_cm3_1_13'] = 6e-16 * scale_e
    k['e_CH3_CH3Minus_cm3_1_14'] = 5e-13 * scale_e
    k['e_C2H4_C2H2_H2_cm3_1_15'] = 5e-11 * scale_e
    k['e_C2H2_C2_H2_cm3_1_16'] = 5e-11 * scale_e
    k['e_C2H4_C2H2_H_H_cm3_1_17'] = 2.5e-11 * scale_e
    k['e_C2H6_C2H2_2H2_cm3_1_18'] = 1.5e-11 * scale_e
    # New reactions from audit
    k['e_C2H2_C2H_H_cm3_1_19'] = 1.0e-11 * scale_e
    k['e_C2H4_C2H3_H_cm3_1_20'] = 8.0e-12 * scale_e
    k['e_C2H6_C2H5_H_cm3_1_21'] = 1.2e-11 * scale_e

    # ===================================================================
    # Group 2: Electron-Impact Ionization
    # ===================================================================
    k['e_CH4_CH3Plus_H_cm3_2_1'] = 1e-11 * scale_e
    k['e_CH4_CH4Plus_cm3_2_2'] = 1e-11 * scale_e
    k['e_Ar_ArPlus_cm3_2_3'] = 8e-12 * scale_e
    k['e_ArStar_ArPlus_cm3_2_4'] = 1e-10 * scale_e
    k['e_C2H6_C2H5Plus_H_2e_cm3_2_5'] = 8e-12 * scale_e
    k['e_C2H4_C2H4Plus_2e_cm3_2_6'] = 1.2e-11 * scale_e
    k['e_C2H4_C2H3Plus_H_2e_cm3_2_7'] = 8e-12 * scale_e
    k['e_C2H2_C2HPlus_2e_cm3_2_8'] = 8e-12 * scale_e

    # ===================================================================
    # Group 3: Ar* Reactions
    # ===================================================================
    k['ArStar_CH4_CH3_H_cm3_3_1'] = 5e-10
    k['ArStar_H2_H_H_cm3_3_2'] = 6e-11
    k['ArStar_H2_ArHPlus_H_cm3_3_3'] = 6e-11
    k['ArStar_CH2_CH_H_cm3_3_4'] = 8e-11
    k['ArStar_C_Ar_CStar_cm3_3_5'] = 1e-10
    k['ArStar_H_ArHPlus_e_cm3_3_6'] = 5e-11
    k['ArStar_CH3_Ar_CH2_H_cm3_3_7'] = 1.2e-10
    k['ArStar_C2_Ar_C2_cm3_3_8'] = 1e-10
    k['ArStar_C2H4_Ar_C2H4_cm3_3_9'] = 1.2e-10
    k['ArStar_CH2_Ar_CH2_cm3_3_10'] = 1e-10
    k['ArStar_CH3_Ar_CH3_cm3_3_11'] = 1e-10
    k['ArStar_M_Ar_3_12'] = 1e-13
    k['ArStar_Ar_Ar_Ar_cm3_3_13'] = 1e-10
    k['CH_ArStar_C_H_Ar_cm3_3_14'] = 1e-11
    k['CH_ArStar_CH_Ar_cm3_3_15'] = 5e-11
    k['ArStar_CH3_CH2_H_Ar_cm3_3_16'] = 8e-11
    k['ArStar_e_Ar_e_cm3_3_17'] = 1e-9
    k['ArStar_H2_Ar_H2Star_cm3_3_18'] = 5e-11
    k['ArStar_CH4_Ar_CH4Star_cm3_3_19'] = 5e-11
    k['ArStar_CH3_Ar_CH3Star_cm3_3_20'] = 6e-11
    k['ArStar_H2_Ar_H_H_cm3_3_21'] = 6e-11
    k['ArStar_C2H2_Ar_C2H2Star_cm3_3_22'] = 6e-11
    k['ArStar_CH3_Ar_CH2_H_cm3_3_23'] = 6e-11
    k['ArStar_C2H4_Ar_C2H4Star_cm3_3_24'] = 5e-11

    # ===================================================================
    # Group 4: Penning Ionization
    # ===================================================================
    k['ArStar_CH4_CH4Plus_cm3_4_1'] = 4.8e-11
    k['ArStar_CH4_CH3Plus_H_cm3_4_2'] = 4e-11
    k['ArStar_Ar_ArPlus_cm3_4_3'] = 4e-11
    k['ArStar_CH4_ArPlus_CH3_H_e_cm3_4_4'] = 2.4e-11
    k['ArStar_CH3_ArPlus_CH2_H_e_cm3_4_5'] = 2.4e-11
    k['ArStar_H_ArPlus_HMinus_cm3_4_6'] = 4e-11
    k['ArStar_CH2_ArPlus_CH_H_e_cm3_4_7'] = 2e-11
    k['ArStar_H2_ArPlus_H2_e_cm3_4_8'] = 5e-11
    k['ArStar_C2H2_ArPlus_C2H2_e_cm3_4_9'] = 5e-11
    k['ArStar_C2H5_ArPlus_C2H5_e_cm3_4_10'] = 5e-11
    k['ArStar_H_ArPlus_H_e_cm3_4_11'] = 5e-11
    k['ArStar_CH4_ArPlus_CH4_e_cm3_4_12'] = 6e-11
    k['ArStar_CH3_ArPlus_CH3_e_cm3_4_13'] = 6e-11
    k['ArStar_C2H4_ArPlus_C2H4_e_cm3_4_14'] = 6e-11
    k['ArStar_C2H5_ArPlus_C2H5_e_cm3_4_15'] = 6e-11
    k['ArStar_CH2_ArPlus_CH2_e_cm3_4_16'] = 6e-11
    k['ArStar_C2H6_ArPlus_C2H6_e_cm3_4_17'] = 6e-11

    # ===================================================================
    # Group 5: Ion-Neutral Reactions
    # ===================================================================
    k['ArPlus_CH4_CH3Plus_H_cm3_5_1'] = 9e-10
    k['CH3Plus_CH4_CH5Plus_CH2_cm3_5_2'] = 7.68e-10
    k['ArPlus_CH3_CH3Plus_cm3_5_3'] = 1.2e-9
    k['CH_CH_C2_H2_cm3_5_4'] = 2.16e-10
    k['ArPlus_CH4_Ar_CH4Plus_cm3_5_5'] = 5.6e-10
    k['CH4Plus_H2_CH5Plus_H_cm3_5_6'] = 1.6e-10
    k['ArPlus_H2_ArHPlus_H_cm3_5_7'] = 1.4e-10
    k['ArHPlus_CH4_Ar_CH5Plus_cm3_5_8'] = 1e-9
    k['ArPlus_CH4_ArHPlus_CH3_cm3_5_9'] = 8e-10
    k['CH2_CH3Plus_CH3_CH2Plus_cm3_5_10'] = 1e-9
    k['CH3Plus_CH4_C2H5Plus_H2_cm3_5_11'] = 1e-9
    k['CH5Plus_C2H4_C2H5Plus_CH4_cm3_5_12'] = 8e-10

    # ===================================================================
    # Group 6: Dissociative Recombination
    # ===================================================================
    k['ArPlus_e_Ar_cm3_6_1'] = 1.5e-7
    k['CH3Plus_e_CH3_cm3_6_2'] = 4.5e-7
    k['CH5Plus_e_CH4_H_cm3_6_3'] = 7.5e-7
    k['e_CH4Plus_CH3_H_cm3_6_4'] = 6e-7
    k['CH3Minus_ArPlus_CH3_Ar_cm3_6_5'] = 1.5e-7
    k['CH3Minus_CH4Plus_CH4_CH3_cm3_6_6'] = 1.5e-7
    k['CH3Minus_CH3Plus_CH4_CH2_cm3_6_7'] = 1.5e-7
    k['CH5Plus_e_CH3_H2_cm3_6_8'] = 1.5e-7
    k['e_CH4Plus_CH2_H2_cm3_6_9'] = 6e-7
    k['CH5Plus_e_CH2_H2_H_cm3_6_10'] = 1.5e-7
    k['e_CH4Plus_CH_H2_H_cm3_6_11'] = 6e-7
    k['CH5Plus_e_CH3_2H_cm3_6_12'] = 1.5e-7
    k['e_CH4Plus_C_2H2_cm3_6_13'] = 5e-7
    k['C2H5Plus_e_C2H4_H_cm3_6_14'] = 3e-7
    k['C2H4Plus_e_C2H2_H2_cm3_6_15'] = 3e-7
    k['C2H3Plus_e_C2H2_H_cm3_6_16'] = 3e-7
    k['HMinus_ArPlus_H_Ar_cm3_6_17'] = 1.8e-7
    k['C2HPlus_e_C2_H_cm3_6_18'] = 3.6e-7
    k['HMinus_CH5Plus_CH4_H2_H_cm3_6_19'] = 1.25e-7
    k['CH4Plus_HMinus_CH4_H_cm3_6_20'] = 1.25e-7
    k['CH3Plus_HMinus_CH4_H2_cm3_6_21'] = 1.25e-7
    k['C2H5Plus_HMinus_C2H6_H_cm3_6_22'] = 1.25e-7
    k['ArHPlus_HMinus_Ar_H2_H_cm3_6_23'] = 1.25e-7
    k['CH5Plus_CH3Minus_CH4_CH4_H_cm3_6_24'] = 1.25e-7
    k['CH4Plus_CH3Minus_CH4_CH3_H_cm3_6_25'] = 1.25e-7
    k['CH3Plus_CH3Minus_CH4_CH2_H_cm3_6_26'] = 1.25e-7
    k['C2H5Plus_CH3Minus_C2H6_H_cm3_6_27'] = 1.25e-7
    k['C2H5Plus_e_C2H4_H_cm3_6_28'] = 3.6e-7

    # ===================================================================
    # Group 7: Neutral-Neutral Reactions (Baulch 2005, updated)
    # ===================================================================
    k['CH2_H_CH_H2_cm3_7_1'] = 1.0e-11
    k['CH2_H_C_H2_H_cm3_7_2'] = 1.2e-11
    k['CH_H_C_H2_cm3_7_3'] = 1.2e-10
    k['C_CH_C2_H_cm3_7_4'] = 1.2e-10
    k['CH_CH3_C2H4_cm3_7_5'] = 1.5e-10  # Updated!
    k['C2_H_CH_C_cm3_7_6'] = 9.6e-11
    k['CH_CH2_C2H2_H_cm3_7_7'] = 1.2e-10
    k['C_CH3_C2_H2_H_cm3_7_8'] = 1.2e-10
    k['CH_C_C2_H_cm3_7_9'] = 1.2e-10
    k['CH_CH3_C2H3_H_cm3_7_10'] = 8e-11
    k['CH_Ar_Ar_C_H_cm3_7_11'] = 1e-15
    k['C_H_CH_cm3_7_12'] = 8e-11
    k['CH2_CH2_C2H4_cm3_7_13'] = 1e-10
    k['CH3_CH2_C2H5_cm3_7_14'] = 8e-11
    k['CH2_CH2_C2H2_H2_cm3_7_15'] = 1e-10
    k['CH3_CH_C2H2_H2_cm3_7_16'] = 1.2e-10
    k['CH2_C_C2H2_cm3_7_17'] = 1e-10
    k['CH_C2H4_C2H2_CH3_cm3_7_18'] = 1e-10
    k['C2H2_C_C2_CH2_cm3_7_19'] = 1e-10
    k['CH_CH4_C2H4_H_cm3_7_20'] = 1.5e-10  # Updated!
    k['CH_H_CH2_cm3_7_21'] = 1e-10
    k['CH_C2H2_C3H2_H_cm3_7_22'] = 1e-10
    k['CH_CH3_C2H2_H2_cm3_7_23'] = 1e-10
    k['CH_C_C2_H2_cm3_7_24'] = 1e-10
    k['H_CH4_CH3_H2_cm3_7_25'] = 6e-12
    k['CH2_CH_C2_H2_H_cm3_7_26'] = 1.2e-10
    k['CH_C2H2_C3H_H2_cm3_7_27'] = 1e-10
    k['CH_C3H_C4H2_H_cm3_7_28'] = 1e-10
    k['CH_C2H2_C2H_CH2_cm3_7_29'] = 1e-10
    k['CH_H2_CH2_H_cm3_7_30'] = 1e-11
    k['CH_C2H3_C3H3_H_cm3_7_31'] = 1e-10
    k['CH_C2H4_C3H4_H_cm3_7_32'] = 1e-10
    k['CH_C2_C3_H_cm3_7_33'] = 1e-10
    k['CH_C2H5_C3H5_H_cm3_7_34'] = 1e-10
    k['CH_C3H2_C4H_H2_cm3_7_35'] = 1e-10
    k['CH3_H_CH2_H2_cm3_7_36'] = 6e-12
    k['CH_C2H6_C3H6_H_cm3_7_37'] = 1.2e-10
    k['CH3_CH3_CH2_CH4_cm3_7_38'] = 1.2e-11
    k['CH_CH4_CH2_CH3_cm3_7_39'] = 1e-11
    k['CH3_CH3_C2H6_cm3_7_40'] = 1e-11
    k['CH_C2H5_C3H6_cm3_7_41'] = 1e-10
    k['CH2_CH2_C2H2_H2_cm3_7_42'] = 1e-11
    k['CH_C_C2_H_cm3_7_43'] = 1e-10
    k['CH_CH_C2_H2_cm3_7_44'] = 1.0e-10  # Confirmed
    k['CH_C2H6_C3H6_H_cm3_7_45'] = 1e-11
    k['CH_C2H4_C2H2_CH3_cm3_7_46'] = 1e-10
    k['C2H_H_C2_H2_cm3_7_47'] = 1e-10
    k['CH_CH2_C2H2_H_cm3_7_48'] = 1e-10
    k['CH3_CH3_C2H2_H2_H2_cm3_7_49'] = 1e-11
    k['C2H2_H_C2_H2_H_cm3_7_50'] = 1e-11
    k['CH_H2_CH2_H_cm3_7_51'] = 1e-11
    k['C2_CH_C3_H_cm3_7_52'] = 1e-10
    k['CH2_C2H3_C2H2_CH3_cm3_7_53'] = 1e-10
    k['C_C2H3_C2_CH3_cm3_7_54'] = 1e-10
    k['CH_C2H5_C3H6_cm3_7_55'] = 1e-10
    k['C2H2_CH_C3_H2_H_cm3_7_56'] = 1e-10
    k['CH_C2H6_C2H2_CH3_H_cm3_7_57'] = 1e-11
    k['CH2_CH2_C2_H2_H2_cm3_7_58'] = 1e-11
    k['CH_CH4_C2H4_H_cm3_7_59'] = 6e-12
    k['CH_C2H4_C3H4_H_cm3_7_60'] = 1e-11
    k['CH3_C2H5_C2H2_CH3_H2_cm3_7_61'] = 1e-11
    k['CH2_CH3_C2H2_H_H2_cm3_7_62'] = 1.2e-11
    k['CH2_C2H5_C2H2_CH3_H_cm3_7_63'] = 1.2e-11
    # New reactions from audit
    k['C_C_M_C2_M_cm6_7_64'] = 1.0e-32
    k['H_C2H4_C2H3_H2_cm3_7_65'] = 1.0e-11

    # ===================================================================
    # Group 8: Termolecular Recombination
    # ===================================================================
    k['H_H_M_H2_M_cm6_8_1'] = 1.0e-32
    k['CH3_CH3_M_C2H6_M_cm6_8_2'] = 3.6e-29
    k['CH3_H_M_CH4_M_cm6_8_3'] = 5.0e-31  # New!

    # ===================================================================
    # Group 9: Stick Reactions
    # ===================================================================
    k['stick_H_9_1'] = 3.89e2
    k['stick_CH3_9_2'] = 3.51e3
    k['stick_CH_9_3'] = 6.25e3
    k['stick_ArPlus_9_4'] = 7.14e3
    k['stick_ArStar_9_5'] = 3.57e2
    k['stick_CH3Plus_9_6'] = 1.16e4
    k['stick_CH5Plus_9_7'] = 1.09e4
    k['stick_ArHPlus_9_8'] = 7.14e3
    k['stick_C2_9_9'] = 1.25e3
    k['stick_C_9_10'] = 6.25e3
    k['stick_C2H2_9_11'] = 5e2
    k['stick_C2H4_9_12'] = 5e2
    k['stick_CH2_9_13'] = 2e3
    k['stick_C2H6_9_14'] = 3.2e2
    k['stick_CH3Minus_9_15'] = 1.5e3
    k['stick_H2_9_16'] = 6.25e2
    k['stick_C2H5_9_17'] = 2e3
    k['stick_HMinus_9_18'] = 5e3
    k['stick_C3H_9_19'] = 1.25e3
    k['stick_C4H2_9_20'] = 1.25e3
    k['stick_C3H3_9_21'] = 1.25e3
    k['stick_C3H4_9_22'] = 1.25e3
    k['stick_C3_9_23'] = 1.25e3
    k['stick_C4H_9_24'] = 1.25e3
    k['stick_C3H6_9_25'] = 1e3
    k['stick_H3Plus_9_26'] = 5e3
    k['stick_CHPlus_9_27'] = 5e3
    k['stick_C2HPlus_9_28'] = 5e3

    # ===================================================================
    # Group 10: Drift Losses (HIGH E-FIELD!)
    # ===================================================================
    k['drift_ArPlus_10_1'] = mobilities['ArPlus'] * E_field / L_discharge
    k['drift_CH4Plus_10_2'] = mobilities['CH4Plus'] * E_field / L_discharge
    k['drift_CH3Plus_10_3'] = mobilities['CH3Plus'] * E_field / L_discharge
    k['drift_CH5Plus_10_4'] = mobilities['CH5Plus'] * E_field / L_discharge
    k['drift_ArHPlus_10_5'] = mobilities['ArHPlus'] * E_field / L_discharge
    k['drift_CH2Plus_10_6'] = mobilities['CH2Plus'] * E_field / L_discharge
    k['drift_C2H5Plus_10_7'] = mobilities['C2H5Plus'] * E_field / L_discharge
    k['drift_C2H4Plus_10_8'] = mobilities['C2H4Plus'] * E_field / L_discharge
    k['drift_C2H3Plus_10_9'] = mobilities['C2H3Plus'] * E_field / L_discharge
    k['drift_H3Plus_10_10'] = mobilities['H3Plus'] * E_field / L_discharge
    k['drift_CHPlus_10_11'] = mobilities['CHPlus'] * E_field / L_discharge
    k['drift_CH3Minus_10_12'] = mobilities['CH3Minus'] * E_field / L_discharge
    k['drift_C2HPlus_10_13'] = mobilities['C2HPlus'] * E_field / L_discharge

    # ===================================================================
    # Group 11: Loss Reactions (BASED ON L_diff!)
    # Formula: k_loss = D / L_diff²
    # For H on copper: k_loss_eff = gamma_H × (D / L_diff²)
    # ===================================================================
    # Get wall recombination coefficient for H on copper
    gamma_H = params.get('gamma_H', 0.01)  # Default: 1% sticks on copper, 99% reflects

    # Calculate wall losses from diffusion length
    for species_name, D_400K in D_ref.items():
        D = D_400K * D_scale
        k_loss_diffusion = D / (L_diff ** 2)

        # For H atoms on copper: only fraction gamma_H actually recombines and is lost
        # Fraction (1-gamma_H) reflects back to gas phase
        if species_name == 'H':
            k_loss = gamma_H * k_loss_diffusion
        else:
            k_loss = k_loss_diffusion  # Other species: assume full loss

        # Map to reaction names
        loss_map = {
            'H': None,  # H uses stick reaction instead
            'CH': 'loss_CH_11_9',
            'CH2': 'loss_CH2_11_1',
            'CH3': 'loss_CH3_11_21',
            'CH4': 'loss_CH4_11_6',
            'C': 'loss_C_11_8',
            'C2': 'loss_C2_11_3',
            'C2H': 'loss_C2H_11_12',
            'C2H2': 'loss_C2H2_11_19',
            'C2H3': 'loss_C2H3_11_22',
            'C2H4': 'loss_C2H4_11_20',
            'C2H5': 'loss_C2H5_11_17',
            'C2H6': 'loss_C2H6_11_5',
            'H2': 'loss_H2_11_2',
            'Ar': 'loss_Ar_11_7',
            'e': 'loss_e_11_4',
            'C3H': 'loss_C3H_11_10',
            'C3H2': 'loss_C3H2_11_23',
            'C3H3': 'loss_C3H3_11_13',
            'C3H4': 'loss_C3H4_11_18',
            'C3H5': 'loss_C3H5_11_24',
            'C3H6': 'loss_C3H6_11_16',
            'C3': 'loss_C3_11_14',
            'C4H': 'loss_C4H_11_15',
            'C4H2': 'loss_C4H2_11_11',
        }

        if species_name in loss_map and loss_map[species_name] is not None:
            k[loss_map[species_name]] = k_loss

    # Excited state losses
    k['loss_C2H2Star_11_25'] = D_ref['C2H2'] * D_scale / (L_diff ** 2)

    return k
