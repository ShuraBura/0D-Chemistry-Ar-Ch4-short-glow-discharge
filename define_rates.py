"""
define_rates.py for Ar/CH4 plasma with temperature-dependent rates
Updated to include Te and Tgas dependence (backward compatible with Te=1 eV, Tgas=400 K)
Sources: Morgan (1992), Janev & Reiter (2002), Phelps (1999), Anicich (2003), Baulch et al. (2005), UMIST (2012)
"""

import numpy as np


def define_rates(params):
    """
    Define all reaction rate constants for the plasma simulation.

    NEW: Temperature-dependent rates using Te and Tgas parameters
    - Te: electron temperature (eV), default = 1.0 eV
    - Tgas: gas temperature (K), default = 400 K (or from params['Tg'])

    Temperature scaling:
    - Electron-impact dissociation: k ~ sqrt(Te) * exp(-E_threshold/Te)
    - Ionization: k ~ sqrt(Te) * exp(-E_ion/Te)
    - Recombination: k ~ Te^(-alpha), alpha ~ 0.5-1.0
    - Neutral-neutral: Independent of Te (thermal energies)
    """
    k = {}

    # Extract parameters
    E_field = params['E_field']
    L_discharge = params['L_discharge']
    mobilities = params['mobilities']

    # Temperature parameters (backward compatible defaults)
    Te = params.get('Te', 1.0)  # eV
    Tgas = params.get('Tgas', params.get('Tg', 400))  # K

    # Calculate total gas density from pressure (for three-body reactions)
    P_mTorr = params.get('P', 500.0)  # mTorr
    P_Pa = P_mTorr * 0.133322
    k_B = 1.380649e-23  # J/K
    n_total_m3 = P_Pa / (k_B * Tgas)  # m⁻³
    n_total = n_total_m3 * 1e-6  # Convert to cm⁻³

    # ===================================================================
    # Temperature Scaling Functions
    # ===================================================================
    def scale_electron_impact(k_ref, Te, Te_ref=1.0, E_threshold=None):
        """
        Scale electron-impact rates with Te.
        For threshold processes: k ~ sqrt(Te) * exp(-E_threshold/Te)
        For high-energy processes: k ~ Te^0.7

        Args:
            k_ref: Reference rate at Te_ref (cm³/s)
            Te: Electron temperature (eV)
            Te_ref: Reference electron temperature (eV), default 1.0
            E_threshold: Threshold energy (eV), if applicable
        """
        if E_threshold is not None and E_threshold > 0:
            # Threshold process - exponential + power law
            return k_ref * np.sqrt(Te/Te_ref) * np.exp(-E_threshold * (1/Te - 1/Te_ref))
        else:
            # High-energy process - power law scaling
            return k_ref * (Te/Te_ref)**0.7

    def scale_ionization(k_ref, Te, Te_ref=1.0, E_ion=12.0):
        """
        Scale ionization rates with Te using threshold behavior.
        k ~ sqrt(Te) * exp(-E_ion/Te)

        Args:
            k_ref: Reference rate at Te_ref (cm³/s)
            Te: Electron temperature (eV)
            Te_ref: Reference electron temperature (eV), default 1.0
            E_ion: Ionization threshold (eV)
        """
        return k_ref * np.sqrt(Te/Te_ref) * np.exp(-E_ion * (1/Te - 1/Te_ref))

    def scale_recombination(k_ref, Te, Te_ref=1.0, alpha=0.7):
        """
        Scale dissociative recombination rates with Te.
        k ~ Te^(-alpha) where alpha typically 0.5 to 1.5

        Args:
            k_ref: Reference rate at Te_ref (cm³/s)
            Te: Electron temperature (eV)
            Te_ref: Reference electron temperature (eV), default 1.0
            alpha: Temperature exponent (dimensionless)
        """
        return k_ref * (Te/Te_ref)**(-alpha)

    # ===================================================================
    # Group 1: Electron-Impact Reactions (Neutral Products)
    # NOW TEMPERATURE-DEPENDENT via Te scaling
    # ===================================================================
    # Reference values are for Te = 1 eV
    k['e_CH4_CH3_H_cm3_1_1'] = scale_electron_impact(4.2e-11, Te, E_threshold=8.5)
    k['e_CH4_CH2_H2_cm3_1_2'] = scale_electron_impact(1.1e-11, Te, E_threshold=9.5)
    k['e_CH4_CH_H2_H_vib_cm3_1_3'] = scale_electron_impact(0.7e-11, Te, E_threshold=10.5)
    k['e_H2_H_H_cm3_1_4'] = scale_electron_impact(6e-12, Te, E_threshold=8.8)
    k['e_CH3_CH2_H_cm3_1_5'] = scale_electron_impact(3e-11, Te, E_threshold=7.5)
    k['e_C2H4_C2H2_H2_cm3_1_6'] = scale_electron_impact(1.0e-11, Te, E_threshold=7.0)
    k['e_Ar_ArStar_cm3_1_7'] = scale_electron_impact(6e-11, Te, E_threshold=11.5)
    k['e_C2H6_C2H4_H2_cm3_1_8'] = scale_electron_impact(7.0e-12, Te, E_threshold=7.5)
    k['e_C2H6_C2H4_H2_e_cm3_1_9'] = scale_electron_impact(1.2e-11, Te, E_threshold=7.5)
    k['e_CH4_CH3Minus_H_cm3_1_10'] = scale_electron_impact(6e-18, Te, E_threshold=8.0)  # Attachment
    k['e_CH4_CH_H_H2_cm3_1_11'] = scale_electron_impact(2e-11, Te, E_threshold=11.0)
    k['e_CH_CH_C_H_e_cm3_1_12'] = scale_electron_impact(6e-11, Te, E_threshold=8.0)
    k['e_H2_HMinus_H_cm3_1_13'] = scale_electron_impact(6e-16, Te, E_threshold=3.75)  # Attachment
    k['e_CH3_CH3Minus_cm3_1_14'] = scale_electron_impact(5e-13, Te, E_threshold=1.0)  # Attachment
    k['e_C2H4_C2H2_H2_cm3_1_15'] = scale_electron_impact(5e-11, Te, E_threshold=7.0)
    k['e_C2H2_C2_H2_cm3_1_16'] = scale_electron_impact(5e-11, Te, E_threshold=9.0)
    k['e_C2H4_C2H2_H_H_cm3_1_17'] = scale_electron_impact(2.5e-11, Te, E_threshold=8.0)
    k['e_C2H6_C2H2_2H2_cm3_1_18'] = scale_electron_impact(1.5e-11, Te, E_threshold=9.0)
    k['e_CH4_CH3_HMinus_cm3_8_1'] = scale_electron_impact(1.0e-13, Te, E_threshold=3.0)  # NEW IN V7: Dissociative attachment
    # New reactions from audit
    k['e_C2H2_C2H_H_cm3_1_19'] = scale_electron_impact(1.0e-11, Te, E_threshold=8.0)
    k['e_C2H4_C2H3_H_cm3_1_20'] = scale_electron_impact(8.0e-12, Te, E_threshold=7.5)
    k['e_C2H6_C2H5_H_cm3_1_21'] = scale_electron_impact(1.2e-11, Te, E_threshold=7.5)

    # ===================================================================
    # Group 2: Electron-Impact Ionization
    # NOW TEMPERATURE-DEPENDENT via Te scaling
    # ===================================================================
    k['e_CH4_CH3Plus_H_cm3_2_1'] = scale_ionization(1e-11, Te, E_ion=12.6)
    k['e_CH4_CH4Plus_cm3_2_2'] = scale_ionization(1e-11, Te, E_ion=12.6)
    k['e_Ar_ArPlus_cm3_2_3'] = scale_ionization(8e-12, Te, E_ion=15.76)
    k['e_ArStar_ArPlus_cm3_2_4'] = scale_ionization(1e-10, Te, E_ion=4.2)  # Ar* already excited
    k['e_C2H6_C2H5Plus_H_2e_cm3_2_5'] = scale_ionization(8e-12, Te, E_ion=11.5)
    k['e_C2H4_C2H4Plus_2e_cm3_2_6'] = scale_ionization(1.2e-11, Te, E_ion=10.5)
    k['e_C2H4_C2H3Plus_H_2e_cm3_2_7'] = scale_ionization(8e-12, Te, E_ion=10.5)
    k['e_C2H2_C2HPlus_2e_cm3_2_8'] = scale_ionization(8e-12, Te, E_ion=11.4)
    k['e_H2_H2Plus_2e_cm3_2_9'] = scale_ionization(3e-12, Te, E_ion=15.43)  # H2 ionization
    k['e_CH_CHPlus_2e_cm3_2_10'] = scale_ionization(5e-12, Te, E_ion=10.64)  # CH ionization (NEW!)

    # ===================================================================
    # Group 3: Ar* Reactions (Temperature-independent - thermal)
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
    # Group 4: Penning Ionization (Temperature-independent - excitation energy driven)
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
    # Group 5: Ion-Neutral Reactions (Temperature-independent - thermal)
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
    # H3+ chemistry (CRITICAL!)
    k['H2Plus_H2_H3Plus_H_cm3_5_13'] = 2.0e-9  # Very fast H3+ formation
    k['H3Plus_CH4_CH5Plus_H2_cm3_5_14'] = 1.5e-9  # Proton transfer
    k['H3Plus_H2_H2Plus_H2_cm3_5_15'] = 6.4e-10  # Reverse formation (minor)

    # ADDED: CH + Ar⁺ → products (fast ion-neutral reaction)
    k['ArPlus_CH_CHPlus_Ar_cm3_5_16'] = 5.0e-9  # Charge transfer / ionization (typical 1e-9 to 1e-8)

    # ===================================================================
    # Group 6: Dissociative Recombination
    # NOW TEMPERATURE-DEPENDENT via Te scaling
    # ===================================================================
    # Rates scale as Te^(-alpha) where alpha ~ 0.5-1.0 for molecular ions
    k['ArPlus_e_Ar_cm3_6_1'] = scale_recombination(1.5e-7, Te, alpha=0.7)
    k['CH3Plus_e_CH3_cm3_6_2'] = scale_recombination(4.5e-7, Te, alpha=0.7)
    k['CH5Plus_e_CH4_H_cm3_6_3'] = scale_recombination(7.5e-7, Te, alpha=0.7)
    k['e_CH4Plus_CH3_H_cm3_6_4'] = scale_recombination(6e-7, Te, alpha=0.7)
    k['CH3Minus_ArPlus_CH3_Ar_cm3_6_5'] = 1.5e-7  # Ion-ion recombination - no Te dependence
    k['CH3Minus_CH4Plus_CH4_CH3_cm3_6_6'] = 1.5e-7  # Ion-ion recombination
    k['CH3Minus_CH3Plus_CH4_CH2_cm3_6_7'] = 1.5e-7  # Ion-ion recombination
    k['CH5Plus_e_CH3_H2_cm3_6_8'] = scale_recombination(1.5e-7, Te, alpha=0.7)
    k['e_CH4Plus_CH2_H2_cm3_6_9'] = scale_recombination(6e-7, Te, alpha=0.7)
    k['CH5Plus_e_CH2_H2_H_cm3_6_10'] = scale_recombination(1.5e-7, Te, alpha=0.7)
    k['e_CH4Plus_CH_H2_H_cm3_6_11'] = scale_recombination(6e-7, Te, alpha=0.7)
    k['CH5Plus_e_CH3_2H_cm3_6_12'] = scale_recombination(1.5e-7, Te, alpha=0.7)
    k['e_CH4Plus_C_2H2_cm3_6_13'] = scale_recombination(5e-7, Te, alpha=0.7)
    k['C2H5Plus_e_C2H4_H_cm3_6_14'] = scale_recombination(3e-7, Te, alpha=0.75)
    k['C2H4Plus_e_C2H2_H2_cm3_6_15'] = scale_recombination(3e-7, Te, alpha=0.75)
    k['C2H3Plus_e_C2H2_H_cm3_6_16'] = scale_recombination(3e-7, Te, alpha=0.75)
    k['HMinus_ArPlus_H_Ar_cm3_6_17'] = 1.8e-7  # Ion-ion recombination
    k['C2HPlus_e_C2_H_cm3_6_18'] = scale_recombination(3.6e-7, Te, alpha=0.75)
    k['HMinus_CH5Plus_CH4_H2_H_cm3_6_19'] = 1.25e-7  # Ion-ion recombination
    k['CH4Plus_HMinus_CH4_H_cm3_6_20'] = 1.25e-7  # Ion-ion recombination
    k['CH3Plus_HMinus_CH4_H2_cm3_6_21'] = 1.25e-7  # Ion-ion recombination
    k['C2H5Plus_HMinus_C2H6_H_cm3_6_22'] = 1.25e-7  # Ion-ion recombination
    k['ArHPlus_HMinus_Ar_H2_H_cm3_6_23'] = 1.25e-7  # Ion-ion recombination
    k['CH5Plus_CH3Minus_CH4_CH4_H_cm3_6_24'] = 1.25e-7  # Ion-ion recombination
    k['CH4Plus_CH3Minus_CH4_CH3_H_cm3_6_25'] = 1.25e-7  # Ion-ion recombination
    k['CH3Plus_CH3Minus_CH4_CH2_H_cm3_6_26'] = 1.25e-7  # Ion-ion recombination
    k['C2H5Plus_CH3Minus_C2H6_H_cm3_6_27'] = 1.25e-7  # Ion-ion recombination
    k['C2H5Plus_e_C2H4_H_cm3_6_28'] = scale_recombination(3.6e-7, Te, alpha=0.75)
    k['ArHPlus_e_Ar_H_cm3_6_29'] = scale_recombination(2.0e-7, Te, alpha=0.7)  # NEW IN V7: ArH+ recombination
    # H2+ and H3+ recombination
    k['H2Plus_e_H_H_cm3_6_29'] = scale_recombination(2.3e-8, Te, alpha=0.5)  # H2+ + e → H + H
    k['H3Plus_e_H2_H_cm3_6_30'] = scale_recombination(2.3e-7, Te, alpha=0.5)  # H3+ + e → H2 + H (dominant)
    k['H3Plus_e_H_H_H_cm3_6_31'] = scale_recombination(4.8e-8, Te, alpha=0.5)  # H3+ + e → H + H + H
    k['CHPlus_e_CH_cm3_6_32'] = scale_recombination(3.5e-7, Te, alpha=0.7)  # CHPlus + e → CH (NEW!)

    # ===================================================================
    # Group 7: Neutral-Neutral Reactions (Temperature-independent at thermal energies)
    # EXCEPT: H + CH4, H + C2H4, and CH3 + H have activation barriers
    # ===================================================================
    k['CH2_H_CH_H2_cm3_7_1'] = 1.0e-11
    k['CH2_H_C_H2_H_cm3_7_2'] = 1.2e-11
    k['CH_H_C_H2_cm3_7_3'] = 2.0e-10  # Corrected from 1.2e-10 (Baulch 2005)
    k['C_CH_C2_H_cm3_7_4'] = 1.05e-12  # From literature: 6.3e11 cm³/mol/s = 1.05e-12 cm³/s
    k['CH_CH3_C2H4_cm3_7_5'] = 1.5e-10  # Updated from 8e-11 (Baulch 2005)
    # k['C2_H_CH_C_cm3_7_6'] = 9.6e-11  # DISABLED: Endothermic at T<1000K (NOT in Baulch)
    k['C2_H_CH_C_cm3_7_6'] = 0.0  # Was 93% of C2 destruction - but doesn't occur at 570K!
    k['CH_CH2_C2H2_H_cm3_7_7'] = 1.2e-10
    k['C_CH3_C2_H2_H_cm3_7_8'] = 1.2e-10
    k['CH_C_C2_H_cm3_7_9'] = 1.05e-12  # From literature: 6.3e11 cm³/mol/s = 1.05e-12 cm³/s
    k['CH_CH3_C2H3_H_cm3_7_10'] = 8e-11
    k['CH_Ar_Ar_C_H_cm3_7_11'] = 1e-15
    k['C_H_CH_cm3_7_12'] = 8e-11
    k['CH2_CH2_C2H4_cm3_7_13'] = 1e-10
    k['CH3_CH2_C2H5_cm3_7_14'] = 8e-11
    k['CH2_CH2_C2H2_H2_cm3_7_15'] = 1e-10
    k['CH3_CH_C2H2_H2_cm3_7_16'] = 1.2e-10
    k['CH2_C_C2H2_cm3_7_17'] = 1e-10
    k['CH_C2H4_C2H2_CH3_cm3_7_18'] = 1e-10
    k['C2H2_C_C2_CH2_cm3_7_19'] = 2.0e-10  # Corrected from 1.0e-10 (Baulch 2005)
    # CH + CH4 → C2H4 + H with temperature dependence (Thiesemann et al. 1997)
    # k = 6.7e-11 × (T/293)^(-0.4) cm³/s (barrierless addition, 290-700 K)
    k_CH_CH4_ref = 6.7e-11  # cm³/s at T=293K
    T_ref_CH_CH4 = 293.0    # K
    n_CH_CH4 = -0.4         # negative temperature exponent
    k['CH_CH4_C2H4_H_cm3_7_20'] = k_CH_CH4_ref * (Tgas / T_ref_CH_CH4)**n_CH_CH4
    k['CH_H_CH2_cm3_7_21'] = 1e-10
    k['CH_C2H2_C3H2_H_cm3_7_22'] = 1e-10
    k['CH_CH3_C2H2_H2_cm3_7_23'] = 1e-10

    # NEW: Missing CH hydrogenation reactions (Tsang & Hampson 1986, Baulch 2005)
    k['CH_H2_CH2_H_cm3_7_NEW1'] = 1.0e-11  # CH + H2 → CH2 + H (Tsang & Hampson 1986)
    k['CH2_H2_CH3_H_cm3_7_NEW2'] = 1.0e-11  # CH2 + H2 → CH3 + H (Baulch 2005)
    k['CH_C_C2_H2_cm3_7_24'] = 1.05e-12  # From literature: 6.3e11 cm³/mol/s = 1.05e-12 cm³/s (CH + C → C2 + H2)
    # H + CH4 → CH3 + H2 with activation barrier Ea = 0.5 eV
    k_H_CH4_ref = 6e-12  # cm³/s reference rate
    Ea_H_CH4 = 0.5  # eV activation barrier
    kB_eV = 8.617333e-5  # eV/K
    k['H_CH4_CH3_H2_cm3_7_25'] = k_H_CH4_ref * np.exp(-Ea_H_CH4 / (kB_eV * Tgas))
    # k['CH2_CH_C2_H2_H_cm3_7_26'] = 1.2e-10  # DISABLED: NOT in Baulch ("garbage reaction")
    k['CH2_CH_C2_H2_H_cm3_7_26'] = 0.0  # Was 21% of C2 production - but doesn't exist!
    k['CH_C2H2_C3H_H2_cm3_7_27'] = 1e-10
    k['CH_C3H_C4H2_H_cm3_7_28'] = 1e-10
    k['CH_C2H2_C2H_CH2_cm3_7_29'] = 1e-10
    k['CH_H2_CH2_H_cm3_7_30'] = 1e-11
    k['CH_C2H3_C3H3_H_cm3_7_31'] = 1e-10
    k['CH_C2H4_C3H4_H_cm3_7_32'] = 1e-10
    k['CH_C2_C3_H_cm3_7_33'] = 1e-10
    k['CH_C2H5_C3H5_H_cm3_7_34'] = 1e-10
    k['CH_C3H2_C4H_H2_cm3_7_35'] = 1e-10
    # CH3 + H → CH2 + H2 with activation barrier Ea = 0.65 eV
    k_CH3_H_ref = 6e-12  # cm³/s reference rate
    Ea_CH3_H = 0.65  # eV activation barrier
    k['CH3_H_CH2_H2_cm3_7_36'] = k_CH3_H_ref * np.exp(-Ea_CH3_H / (kB_eV * Tgas))
    k['CH_C2H6_C3H6_H_cm3_7_37'] = 1.2e-10
    k['CH3_CH3_CH2_CH4_cm3_7_38'] = 1.2e-11
    k['CH_CH4_CH2_CH3_cm3_7_39'] = 1e-11
    k['CH3_CH3_C2H6_cm3_7_40'] = 1e-11
    k['CH_C2H5_C3H6_cm3_7_41'] = 1e-10
    k['CH2_CH2_C2H2_H2_cm3_7_42'] = 1e-11
    k['CH_C_C2_H_cm3_7_43'] = 1.05e-12  # From literature: 6.3e11 cm³/mol/s = 1.05e-12 cm³/s
    k['CH_CH_C2_H2_cm3_7_44'] = 1.0e-10  # Confirmed (Baulch 2005)
    k['CH_C2H6_C3H6_H_cm3_7_45'] = 1e-11
    k['CH_C2H4_C2H2_CH3_cm3_7_46'] = 1e-10
    k['C2H_H_C2_H2_cm3_7_47'] = 1e-10
    k['CH_CH2_C2H2_H_cm3_7_48'] = 1e-10
    k['CH3_CH3_C2H2_H2_H2_cm3_7_49'] = 1e-11
    # k['C2H2_H_C2_H2_H_cm3_7_50'] = 1e-11  # OLD: 7.5 billion times too high!
    # Correct rate from Baulch (2005): k = 1.67e-14 × T^1.64 × exp(-15250/T)
    k['C2H2_H_C2_H2_H_cm3_7_50'] = 1.67e-14 * Tgas**1.64 * np.exp(-15250/Tgas)  # Arrhenius form
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
    # New reactions from audit (Baulch 2005, Kushner)
    k['C_C_M_C2_M_cm6_7_64'] = 1.0e-32 * n_total  # C + C + M → C2 + M (three-body)
    # H + C2H4 → C2H3 + H2 with activation barrier Ea = 0.6 eV
    k_H_C2H4_ref = 1.0e-11  # cm³/s reference rate
    Ea_H_C2H4 = 0.6  # eV activation barrier (H abstraction from ethylene)
    k['H_C2H4_C2H3_H2_cm3_7_65'] = k_H_C2H4_ref * np.exp(-Ea_H_C2H4 / (kB_eV * Tgas))

    # ===================================================================
    # MISSING CH3 PRODUCTION PATHWAYS (Added to address physical realism)
    # ===================================================================
    # High-priority pathways identified from literature (Baulch 2005, NIST)
    k['e_C2H4_CH3_CH_cm3_7_66'] = scale_electron_impact(2.0e-11, Te, E_threshold=8.0)  # e + C2H4 → CH3 + CH
    k['e_C2H6_CH3_CH3_cm3_7_67'] = scale_electron_impact(3.0e-11, Te, E_threshold=7.5)  # e + C2H6 → CH3 + CH3
    k['e_C2H5_CH3_CH2_cm3_7_68'] = scale_electron_impact(1.5e-11, Te, E_threshold=7.0)  # e + C2H5 → CH3 + CH2
    k['ArStar_C2H4_CH3_CH_cm3_7_69'] = 3.0e-10  # Ar* + C2H4 → CH3 + CH
    k['ArStar_C2H6_CH3_CH3_cm3_7_70'] = 5.0e-10  # Ar* + C2H6 → CH3 + CH3
    k['H_C2H5_CH3_CH2_cm3_7_71'] = 8.0e-11  # H + C2H5 → CH3 + CH2 (important: H is abundant!)
    k['CH2_CH2_CH3_CH_cm3_7_72'] = 6.0e-11  # CH2 + CH2 → CH3 + CH
    k['C2H5Plus_e_CH3_CH2_cm3_7_73'] = scale_recombination(3.0e-7, Te, alpha=0.7)  # C2H5+ + e → CH3 + CH2
    k['CH2_H_M_CH3_M_cm6_7_74'] = 2.0e-30 * n_total  # CH2 + H + M → CH3 + M (three-body)
    k['ArPlus_CH4_CH3Plus_ArH_cm3_7_75'] = 2.0e-9  # Ar+ + CH4 → CH3+ + Ar + H (charge transfer)

    # ===================================================================
    # Group 8: Termolecular Recombination (Temperature-independent)
    # ===================================================================
    # Three-body rate constants (cm⁶/s) are multiplied by n_total (cm⁻³) to give effective rate (cm³/s)
    k['H_H_M_H2_M_cm6_8_1'] = 1.0e-32 * n_total  # H + H + M → H2 + M
    k['CH3_CH3_M_C2H6_M_cm6_8_2'] = 3.6e-29 * n_total  # CH3 + CH3 + M → C2H6 + M
    k['CH3_H_M_CH4_M_cm6_8_3'] = 5.0e-31 * n_total  # CH3 + H + M → CH4 + M

    # Three-body electron-ion recombination (provides stabilization at high densities)
    # These are CRITICAL for preventing runaway ionization
    # Rate constants from literature (Flannery 1969, Bates 1962)
    k['e_ArPlus_M_Ar_M_cm6_8_4'] = 1.0e-25 * n_total  # e + Ar+ + M → Ar + M
    k['e_CH4Plus_M_CH4_M_cm6_8_5'] = 1.0e-25 * n_total  # e + CH4+ + M → CH4 + M
    k['e_CH3Plus_M_CH3_M_cm6_8_6'] = 1.0e-25 * n_total  # e + CH3+ + M → CH3 + M
    k['e_CH5Plus_M_CH5_M_cm6_8_7'] = 1.0e-25 * n_total  # e + CH5+ + M → CH5 + M (then dissociates)
    k['e_ArHPlus_M_ArH_M_cm6_8_8'] = 1.0e-25 * n_total  # e + ArH+ + M → ArH + M
    k['e_C2H5Plus_M_C2H5_M_cm6_8_9'] = 1.0e-26 * n_total  # e + C2H5+ + M → C2H5 + M

    # ===================================================================
    # Group 9: Stick Reactions (Temperature-independent wall sticking)
    # ===================================================================
    k['stick_H_9_1'] = 3.89e2
    k['stick_CH3_9_2'] = 3.51e3
    k['stick_CH_9_3'] = 6.25e3
    k['stick_ArPlus_9_4'] = 7.14e3
    k['stick_ArStar_9_5'] = 3.57e2
    k['stick_CH3Plus_9_6'] = 1.16e4
    k['stick_CH5Plus_9_7'] = 1.09e4
    k['stick_ArHPlus_9_8'] = 7.14e3
    k['stick_C2_9_9'] = 1.25e2  # Changed from 1.25e3 (γ: 0.01 → 0.001, literature-based)
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
    k['stick_H2Plus_9_29'] = 5e3  # H2+ wall loss

    # ===================================================================
    # Group 10: Drift Losses (Temperature-independent)
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
    k['drift_H2Plus_10_14'] = mobilities.get('H2Plus', 5000) * E_field / L_discharge  # H2+ drift (with default)
    k['drift_HMinus_10_15'] = mobilities.get('HMinus', 3000) * E_field / L_discharge  # HMinus drift (with default)

    # ===================================================================
    # Group 11: Loss Reactions (Temperature-independent diffusion losses)
    # Note: These could be made Tgas-dependent via D ~ T^1.75 if needed
    # ===================================================================
    k['loss_CH2_11_1'] = 5e2
    k['loss_H2_11_2'] = 4e2
    k['loss_C2_11_3'] = 2e2
    k['loss_e_11_4'] = 1e3  # Electrons slightly higher
    k['loss_C2H6_11_5'] = 1e2
    k['loss_CH4_11_6'] = 2e2
    k['loss_Ar_11_7'] = 2e2
    k['loss_C_11_8'] = 4e2
    k['loss_CH_11_9'] = 8e2  # Reactive radical
    k['loss_C3H_11_10'] = 2e2
    k['loss_C4H2_11_11'] = 2e2
    k['loss_C2H_11_12'] = 3e2
    k['loss_C3H3_11_13'] = 2e2
    k['loss_C3_11_14'] = 2e2
    k['loss_C4H_11_15'] = 2e2
    k['loss_C3H6_11_16'] = 1e2
    k['loss_C2H5_11_17'] = 3e2
    k['loss_C3H4_11_18'] = 2e2
    k['loss_C2H2_11_19'] = 2e2
    k['loss_C2H4_11_20'] = 2e2
    k['loss_CH3_11_21'] = 3e2
    k['loss_C2H3_11_22'] = 3e2
    k['loss_C3H2_11_23'] = 2e2
    k['loss_C3H5_11_24'] = 2e2
    k['loss_C2H2Star_11_25'] = 3e2

    # ===================================================================
    # DUST/NANOPARTICLE LOSS TERMS
    # ===================================================================
    # Loss rate: k_dust = n_dust × π × r_dust² × v_thermal × α_stick
    # For CH: v_th ≈ 7e4 cm/s (300 K)
    # Configurable via params['dust_density'], params['dust_radius'], params['dust_sticking']
    # Default: moderate dust scenario (n=1e8 cm⁻³, r=50 nm, α=0.5)
    #
    # To enable/disable: set params['enable_dust_loss'] = True/False (default False)
    # To tune: params['dust_multiplier'] = 0.1 to 10.0 (default 1.0)

    enable_dust = params.get('enable_dust_loss', False)
    dust_multiplier = params.get('dust_multiplier', 1.0)

    if enable_dust:
        # Dust parameters (can be customized via params)
        n_dust = params.get('dust_density', 1e8)  # cm⁻³
        r_dust = params.get('dust_radius', 50e-7)  # cm (50 nm default)
        alpha_dust = params.get('dust_sticking', 0.5)  # dimensionless

        # Thermal velocities at 300 K (cm/s)
        v_thermal = {
            'CH': 6.99e4,   # 13 amu
            'CH2': 6.09e4,  # 14 amu
            'CH3': 5.43e4,  # 15 amu
            'C': 8.71e4,    # 12 amu
            'C2': 6.16e4,   # 24 amu
            'H': 1.73e5,    # 1 amu
        }

        # Calculate dust loss coefficients
        for species, v_th in v_thermal.items():
            k_dust_base = n_dust * 3.14159 * r_dust**2 * v_th * alpha_dust
            k_dust = k_dust_base * dust_multiplier

            # Add dust loss terms
            k[f'dust_loss_{species}_12'] = k_dust

    return k
