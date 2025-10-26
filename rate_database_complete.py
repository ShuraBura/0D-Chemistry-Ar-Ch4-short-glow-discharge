"""
rate_database_complete.py - COMPLETE reaction rate database
Auto-generated from MATLAB define_rates.m
All 150+ rates with literature ranges, citations, and flags preserved
"""

class RateConstant:
    """Structure to hold rate constant with metadata."""
    def __init__(self, value, min_val, max_val, source, notes="", flag=""):
        self.value = value
        self.min = min_val
        self.max = max_val
        self.source = source
        self.notes = notes
        self.flag = flag

    def is_within_range(self, val):
        """Check if a value is within literature range."""
        return self.min <= val <= self.max

    def __repr__(self):
        s = f"RateConstant({self.value:.2e}, [{self.min:.2e}, {self.max:.2e}], {self.source})"
        if self.flag:
            s += f" [FLAG: {self.flag}]"
        return s


def get_complete_rate_database():
    """
    Get COMPLETE rate database with ALL rates from MATLAB define_rates.m

    Returns dict of RateConstant objects with:
    - value: Current/nominal value
    - min: Minimum literature value
    - max: Maximum literature value
    - source: Citation
    - notes: Scaling/rationale
    - flag: Validation warnings
    """
    db = {}


    # ========== ELECTRON-IMPACT REACTIONS ==========
    db['e_CH4_CH3_H_cm3_1_1'] = RateConstant(6.00e-11, 4.00e-11, 1.20e-10, "Morgan (1992)", notes="Scaled: 1.2e-10 * exp(-4/1)/exp(-4/2) ≈ 6e-11")
    db['e_CH4_CH2_H2_cm3_1_2'] = RateConstant(3.00e-11, 2.00e-11, 6.00e-11, "Morgan (1992)", notes="Scaled: 6e-11 * exp(-5/1)/exp(-5/2) ≈ 3e-11")
    db['e_CH4_CH_H2_H_vib_cm3_1_3'] = RateConstant(3.00e-11, 2.00e-11, 1.00e-10, "Janev & Reiter (2002)", notes="Scaled: 6e-11 * 0.5 ≈ 3e-11")
    db['e_H2_H_H_cm3_1_4'] = RateConstant(6.00e-12, 4.00e-12, 1.20e-11, "Morgan (1992)", notes="Scaled: 1.2e-11 * exp(-4.5/1)/exp(-4.5/2) ≈ 6e-12")
    db['e_CH3_CH2_H_cm3_1_5'] = RateConstant(3.00e-11, 2.00e-11, 6.00e-11, "Janev & Reiter (2002)", notes="Scaled: 6e-11 * 0.5 ≈ 3e-11")
    db['e_C2H4_C2H2_H2_cm3_1_6'] = RateConstant(1.80e-11, 1.20e-11, 3.60e-11, "Janev & Reiter (2002)", notes="Scaled: 3.6e-11 * 0.5 ≈ 1.8e-11")
    db['e_Ar_ArStar_cm3_1_7'] = RateConstant(6.00e-11, 4.00e-11, 1.80e-10, "Phelps (1999)", notes="Scaled: 1.2e-10 * exp(-11.5/1)/exp(-11.5/2) ≈ 6e-11")
    db['e_C2H6_C2H4_H2_cm3_1_8'] = RateConstant(1.50e-11, 1.20e-11, 3.60e-11, "Janev & Reiter (2002)", notes="Scaled: 3e-11 * 0.5 ≈ 1.5e-11")
    db['e_C2H6_C2H4_H2_e_cm3_1_9'] = RateConstant(1.20e-11, 1.20e-11, 3.60e-11, "Janev & Reiter (2002)", notes="Scaled: 2.4e-11 * 0.5 ≈ 1.2e-11")
    db['e_CH4_CH3Minus_H_cm3_1_10'] = RateConstant(6.00e-18, 4.00e-18, 1.20e-17, "Janev & Reiter (2002)", notes="Scaled: 1.2e-17 * 0.5 ≈ 6e-18", flag="1.10 Rate low, needs LXCat validation")
    db['e_CH4_CH_H_H2_cm3_1_11'] = RateConstant(2.00e-11, 2.00e-11, 6.00e-11, "Janev & Reiter (2002)", notes="Scaled: 4e-11 * 0.5 ≈ 2e-11")
    db['e_CH_CH_C_H_e_cm3_1_12'] = RateConstant(6.00e-11, 4.00e-11, 1.20e-10, "Janev & Reiter (2002)", notes="Scaled: 1.2e-10 * 0.5 ≈ 6e-11")
    db['e_H2_HMinus_H_cm3_1_13'] = RateConstant(6.00e-16, 4.00e-16, 1.20e-15, "Morgan (1992)", notes="Scaled: 1.2e-15 * 0.5 ≈ 6e-16")
    db['e_CH3_CH3Minus_cm3_1_14'] = RateConstant(5.00e-13, 4.00e-13, 1.20e-12, "Janev & Reiter (2002)", notes="Scaled: 1e-12 * 0.5 ≈ 5e-13")
    db['e_C2H4_C2H2_H2_cm3_1_15'] = RateConstant(5.00e-11, 4.00e-11, 1.20e-10, "Janev & Reiter (2002)", notes="Scaled: 1e-10 * 0.5 ≈ 5e-11")
    db['e_C2H2_C2_H2_cm3_1_16'] = RateConstant(5.00e-11, 4.00e-11, 1.20e-10, "Janev & Reiter (2002)", notes="Scaled: 1e-10 * 0.5 ≈ 5e-11")
    db['e_C2H4_C2H2_H_H_cm3_1_17'] = RateConstant(2.50e-11, 2.00e-11, 6.00e-11, "Janev & Reiter (2002)", notes="Scaled: 5e-11 * 0.5 ≈ 2.5e-11")
    db['e_C2H6_C2H2_2H2_cm3_1_18'] = RateConstant(1.50e-11, 1.20e-11, 3.60e-11, "Janev & Reiter (2002)", notes="Scaled: 3e-11 * 0.5 ≈ 1.5e-11")
    db['e_CH4_CH3Plus_H_cm3_2_1'] = RateConstant(1.00e-11, 1.00e-11, 3.84e-11, "Morgan (1992)", notes="Scaled: 2.048e-11 * exp(-12.6/1)/exp(-12.6/2) ≈ 1e-11")
    db['e_CH4_CH4Plus_cm3_2_2'] = RateConstant(1.00e-11, 1.00e-11, 3.84e-11, "Morgan (1992)", notes="Scaled: 2.048e-11 * exp(-12.6/1)/exp(-12.6/2) ≈ 1e-11")
    db['e_Ar_ArPlus_cm3_2_3'] = RateConstant(8.00e-12, 8.00e-12, 2.64e-11, "Phelps (1999)", notes="Scaled: 1.76e-11 * exp(-15.6/1)/exp(-15.6/2) ≈ 8e-12")
    db['e_ArStar_ArPlus_cm3_2_4'] = RateConstant(1.00e-10, 8.00e-11, 2.40e-10, "Phelps (1999)", notes="Scaled: 2e-10 * 0.5 ≈ 1e-10")
    db['e_C2H6_C2H5Plus_H_2e_cm3_2_5'] = RateConstant(8.00e-12, 8.00e-12, 2.40e-11, "Janev & Reiter (2002)", notes="Scaled: 1.6e-11 * 0.5 ≈ 8e-12")
    db['e_C2H4_C2H4Plus_2e_cm3_2_6'] = RateConstant(1.20e-11, 1.20e-11, 3.60e-11, "Janev & Reiter (2002)", notes="Scaled: 2.4e-11 * 0.5 ≈ 1.2e-11")
    db['e_C2H4_C2H3Plus_H_2e_cm3_2_7'] = RateConstant(8.00e-12, 8.00e-12, 2.40e-11, "Janev & Reiter (2002)", notes="Scaled: 1.6e-11 * 0.5 ≈ 8e-12")
    db['e_C2H2_C2HPlus_2e_cm3_2_8'] = RateConstant(8.00e-12, 8.00e-12, 2.40e-11, "Janev & Reiter (2002)", notes="Scaled: 1.6e-11 * 0.5 ≈ 8e-12")

    # ========== Ar* REACTIONS ==========
    db['ArStar_CH4_CH3_H_cm3_3_1'] = RateConstant(5.00e-10, 5.00e-11, 5.00e-10, "Velazco et al. (1978)")
    db['ArStar_H2_H_H_cm3_3_2'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_H2_ArHPlus_H_cm3_3_3'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH2_CH_H_cm3_3_4'] = RateConstant(8.00e-11, 8.00e-11, 1.20e-10, "Phelps (1999)")
    db['ArStar_C_Ar_CStar_cm3_3_5'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Phelps (1999)")
    db['ArStar_H_ArHPlus_e_cm3_3_6'] = RateConstant(5.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH3_Ar_CH2_H_cm3_3_7'] = RateConstant(1.20e-10, 8.00e-11, 1.20e-10, "Phelps (1999)")
    db['ArStar_C2_Ar_C2_cm3_3_8'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Phelps (1999)")
    db['ArStar_C2H4_Ar_C2H4_cm3_3_9'] = RateConstant(1.20e-10, 8.00e-11, 1.20e-10, "Phelps (1999)")
    db['ArStar_CH2_Ar_CH2_cm3_3_10'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Phelps (1999)")
    db['ArStar_CH3_Ar_CH3_cm3_3_11'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Phelps (1999)")
    db['ArStar_M_Ar_3_12'] = RateConstant(1.00e-13, 8.00e-14, 1.20e-13, "Phelps (1999)", flag="3.12 Needs validation")
    db['ArStar_Ar_Ar_Ar_cm3_3_13'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Phelps (1999)")

    # ========== NEUTRAL-NEUTRAL REACTIONS ==========
    db['CH_ArStar_C_H_Ar_cm3_3_14'] = RateConstant(1.00e-11, 1.00e-11, 1.00e-11, "Smith et al. (2018)")
    db['CH_ArStar_CH_Ar_cm3_3_15'] = RateConstant(5.00e-11, 5.00e-11, 5.00e-11, "Tachibana et al. (1989)")

    # ========== Ar* REACTIONS ==========
    db['ArStar_CH3_CH2_H_Ar_cm3_3_16'] = RateConstant(8.00e-11, 8.00e-11, 1.20e-10, "Phelps (1999)")
    db['ArStar_e_Ar_e_cm3_3_17'] = RateConstant(1.00e-09, 8.00e-10, 1.20e-09, "Phelps (1999)")
    db['ArStar_H2_Ar_H2Star_cm3_3_18'] = RateConstant(5.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH4_Ar_CH4Star_cm3_3_19'] = RateConstant(5.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH3_Ar_CH3Star_cm3_3_20'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_H2_Ar_H_H_cm3_3_21'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_C2H2_Ar_C2H2Star_cm3_3_22'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH3_Ar_CH2_H_cm3_3_23'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_C2H4_Ar_C2H4Star_cm3_3_24'] = RateConstant(5.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH4_CH4Plus_cm3_4_1'] = RateConstant(4.80e-11, 4.80e-11, 7.20e-11, "Phelps (1999)")
    db['ArStar_CH4_CH3Plus_H_cm3_4_2'] = RateConstant(4.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_Ar_ArPlus_cm3_4_3'] = RateConstant(4.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH4_ArPlus_CH3_H_e_cm3_4_4'] = RateConstant(2.40e-11, 1.60e-11, 2.40e-11, "Phelps (1999)")
    db['ArStar_CH3_ArPlus_CH2_H_e_cm3_4_5'] = RateConstant(2.40e-11, 1.60e-11, 2.40e-11, "Phelps (1999)")
    db['ArStar_H_ArPlus_HMinus_cm3_4_6'] = RateConstant(4.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH2_ArPlus_CH_H_e_cm3_4_7'] = RateConstant(2.00e-11, 1.60e-11, 2.40e-11, "Phelps (1999)")
    db['ArStar_H2_ArPlus_H2_e_cm3_4_8'] = RateConstant(5.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_C2H2_ArPlus_C2H2_e_cm3_4_9'] = RateConstant(5.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_C2H5_ArPlus_C2H5_e_cm3_4_10'] = RateConstant(5.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_H_ArPlus_H_e_cm3_4_11'] = RateConstant(5.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH4_ArPlus_CH4_e_cm3_4_12'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH3_ArPlus_CH3_e_cm3_4_13'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_C2H4_ArPlus_C2H4_e_cm3_4_14'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_C2H5_ArPlus_C2H5_e_cm3_4_15'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_CH2_ArPlus_CH2_e_cm3_4_16'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")
    db['ArStar_C2H6_ArPlus_C2H6_e_cm3_4_17'] = RateConstant(6.00e-11, 4.00e-11, 6.00e-11, "Phelps (1999)")

    # ========== ION-NEUTRAL REACTIONS ==========
    db['ArPlus_CH4_CH3Plus_H_cm3_5_1'] = RateConstant(9.00e-10, 9.00e-10, 1.30e-09, "Anicich (2003)")
    db['CH3Plus_CH4_CH5Plus_CH2_cm3_5_2'] = RateConstant(7.68e-10, 9.60e-10, 1.44e-09, "Anicich (2003)")
    db['ArPlus_CH3_CH3Plus_cm3_5_3'] = RateConstant(1.20e-09, 8.00e-10, 1.20e-09, "Anicich (2003)")

    # ========== NEUTRAL-NEUTRAL REACTIONS ==========
    db['CH_CH_C2_H2_cm3_5_4'] = RateConstant(2.16e-10, 1.20e-10, 1.80e-10, "Baulch et al. (2005)")

    # ========== ION-NEUTRAL REACTIONS ==========
    db['ArPlus_CH4_Ar_CH4Plus_cm3_5_5'] = RateConstant(5.60e-10, 7.00e-10, 1.40e-09, "Anicich (2003)")
    db['CH4Plus_H2_CH5Plus_H_cm3_5_6'] = RateConstant(1.60e-10, 1.00e-10, 1.60e-10, "Anicich (2003)")
    db['ArPlus_H2_ArHPlus_H_cm3_5_7'] = RateConstant(1.40e-10, 1.40e-10, 2.00e-10, "Anicich (2003)")
    db['ArHPlus_CH4_Ar_CH5Plus_cm3_5_8'] = RateConstant(1.00e-09, 9.00e-10, 1.30e-09, "Anicich (2003)")
    db['ArPlus_CH4_ArHPlus_CH3_cm3_5_9'] = RateConstant(8.00e-10, 8.00e-10, 1.20e-09, "Anicich (2003)")
    db['CH2_CH3Plus_CH3_CH2Plus_cm3_5_10'] = RateConstant(1.00e-09, 8.00e-10, 1.20e-09, "Anicich (2003)")
    db['CH3Plus_CH4_C2H5Plus_H2_cm3_5_11'] = RateConstant(1.00e-09, 8.00e-10, 1.20e-09, "Anicich (2003)")
    db['CH5Plus_C2H4_C2H5Plus_CH4_cm3_5_12'] = RateConstant(8.00e-10, 8.00e-10, 1.20e-09, "Anicich (2003)")
    db['ArPlus_e_Ar_cm3_6_1'] = RateConstant(1.50e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1.2e-7 * (1/2)^(-0.7) ≈ 1.5e-7")
    db['CH3Plus_e_CH3_cm3_6_2'] = RateConstant(4.50e-07, 2.40e-07, 3.60e-07, "UMIST (2012)", notes="Scaled: 3.6e-7 * (1/2)^(-0.7) ≈ 4.5e-7")
    db['CH5Plus_e_CH4_H_cm3_6_3'] = RateConstant(7.50e-07, 4.00e-07, 6.00e-07, "UMIST (2012)", notes="Scaled: 6e-7 * (1/2)^(-0.7) ≈ 7.5e-7")

    # ========== ELECTRON-IMPACT REACTIONS ==========
    # CH4+ dissociative recombination - Updated with Thomas et al. (2013) literature values
    # Total rate: 1.71e-6 cm³/s at 300K, Temperature dependence: (Te/300)^(-0.66)
    # Literature branching: CH3+H(18%), CH2+2H(51%), CH2+H2(6%), CH+H2+H(23%), CH+2H2(2%)
    # Note: Dominant channel CH2+2H (51%) is NOT in model - CH2+H2 used as proxy
    db['e_CH4Plus_CH3_H_cm3_6_4'] = RateConstant(3.08e-07, 1.00e-07, 9.00e-07, "Thomas et al. (2013) + UMIST", notes="Literature: 1.71e-6 * 0.18 = 3.08e-7 (18% branching)")
    db['CH3Minus_ArPlus_CH3_Ar_cm3_6_5'] = RateConstant(1.50e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1.2e-7 * (1/2)^(-0.7) ≈ 1.5e-7", flag="6.5 Needs validation")
    db['CH3Minus_CH4Plus_CH4_CH3_cm3_6_6'] = RateConstant(1.50e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1.2e-7 * (1/2)^(-0.7) ≈ 1.5e-7", flag="6.6 Needs validation")
    db['CH3Minus_CH3Plus_CH4_CH2_cm3_6_7'] = RateConstant(1.50e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1.2e-7 * (1/2)^(-0.7) ≈ 1.5e-7")
    db['CH5Plus_e_CH3_H2_cm3_6_8'] = RateConstant(1.50e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1.2e-7 * (1/2)^(-0.7) ≈ 1.5e-7")
    db['e_CH4Plus_CH2_H2_cm3_6_9'] = RateConstant(9.74e-07, 1.00e-07, 3.00e-06, "Thomas et al. (2013) + UMIST", notes="Literature: 1.71e-6 * (0.06 + 0.51) = 9.74e-7 (CH2+H2 6% + proxy for CH2+2H 51%)")
    db['CH5Plus_e_CH2_H2_H_cm3_6_10'] = RateConstant(1.50e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1.2e-7 * (1/2)^(-0.7) ≈ 1.5e-7")
    db['e_CH4Plus_CH_H2_H_cm3_6_11'] = RateConstant(3.93e-07, 1.00e-07, 1.20e-06, "Thomas et al. (2013) + UMIST", notes="Literature: 1.71e-6 * 0.23 = 3.93e-7 (23% branching)")
    db['CH5Plus_e_CH3_2H_cm3_6_12'] = RateConstant(1.50e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1.2e-7 * (1/2)^(-0.7) ≈ 1.5e-7")
    db['e_CH4Plus_C_2H2_cm3_6_13'] = RateConstant(3.42e-08, 1.00e-08, 1.00e-07, "Thomas et al. (2013) + UMIST", notes="Literature: 1.71e-6 * 0.02 = 3.42e-8 (2% branching)")
    db['C2H5Plus_e_C2H4_H_cm3_6_14'] = RateConstant(3.00e-07, 2.40e-07, 3.60e-07, "UMIST (2012)", notes="Scaled: 2.4e-7 * (1/2)^(-0.7) ≈ 3e-7")
    db['C2H4Plus_e_C2H2_H2_cm3_6_15'] = RateConstant(3.00e-07, 2.40e-07, 3.60e-07, "UMIST (2012)", notes="Scaled: 2.4e-7 * (1/2)^(-0.7) ≈ 3e-7")
    db['C2H3Plus_e_C2H2_H_cm3_6_16'] = RateConstant(3.00e-07, 2.40e-07, 3.60e-07, "UMIST (2012)", notes="Scaled: 2.4e-7 * (1/2)^(-0.7) ≈ 3e-7")
    db['HMinus_ArPlus_H_Ar_cm3_6_17'] = RateConstant(1.80e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1.44e-7 * (1/2)^(-0.7) ≈ 1.8e-7")
    db['C2HPlus_e_C2_H_cm3_6_18'] = RateConstant(3.60e-07, 1.90e-07, 2.90e-07, "UMIST (2012)", notes="Scaled: 2.9e-7 * (1/2)^(-0.7) ≈ 3.6e-7")
    db['HMinus_CH5Plus_CH4_H2_H_cm3_6_19'] = RateConstant(1.25e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1e-7 * (1/2)^(-0.7) ≈ 1.25e-7")
    db['CH4Plus_HMinus_CH4_H_cm3_6_20'] = RateConstant(1.25e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1e-7 * (1/2)^(-0.7) ≈ 1.25e-7")
    db['CH3Plus_HMinus_CH4_H2_cm3_6_21'] = RateConstant(1.25e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1e-7 * (1/2)^(-0.7) ≈ 1.25e-7")
    db['C2H5Plus_HMinus_C2H6_H_cm3_6_22'] = RateConstant(1.25e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1e-7 * (1/2)^(-0.7) ≈ 1.25e-7")
    db['ArHPlus_HMinus_Ar_H2_H_cm3_6_23'] = RateConstant(1.25e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1e-7 * (1/2)^(-0.7) ≈ 1.25e-7")
    db['CH5Plus_CH3Minus_CH4_CH4_H_cm3_6_24'] = RateConstant(1.25e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1e-7 * (1/2)^(-0.7) ≈ 1.25e-7")
    db['CH4Plus_CH3Minus_CH4_CH3_H_cm3_6_25'] = RateConstant(1.25e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1e-7 * (1/2)^(-0.7) ≈ 1.25e-7")
    db['CH3Plus_CH3Minus_CH4_CH2_H_cm3_6_26'] = RateConstant(1.25e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1e-7 * (1/2)^(-0.7) ≈ 1.25e-7")
    db['C2H5Plus_CH3Minus_C2H6_H_cm3_6_27'] = RateConstant(1.25e-07, 8.00e-08, 1.20e-07, "UMIST (2012)", notes="Scaled: 1e-7 * (1/2)^(-0.7) ≈ 1.25e-7")
    db['C2H5Plus_e_C2H4_H_cm3_6_28'] = RateConstant(3.60e-07, 2.40e-07, 3.60e-07, "UMIST (2012)", notes="Scaled: 3.6e-7 * (1/2)^(-0.7) ≈ 4.5e-7")

    # NEW IN V7: ArH+ dissociative recombination (7.8% of ions in V6!)
    db['ArHPlus_e_Ar_H_cm3_6_29'] = RateConstant(2.00e-07, 1.00e-07, 3.00e-07, "Estimated from Ar+ literature", notes="Similar to Ar+ recombination, k~1-3e-7 cm³/s")

    # ========== ATTACHMENT REACTIONS ==========
    # NEW IN V7: Dissociative attachment
    db['e_CH4_CH3_HMinus_cm3_8_1'] = RateConstant(1.00e-13, 1.00e-15, 1.00e-12, "Estimated", notes="Dissociative attachment, very slow at low e⁻ energies")

    db['CH2_H_CH_H2_cm3_7_1'] = RateConstant(1.00e-11, 1.00e-11, 2.25e-11, "Baulch et al. (2005)")
    db['CH2_H_C_H2_H_cm3_7_2'] = RateConstant(1.20e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)", flag="7.2 Needs discussion")

    # ========== NEUTRAL-NEUTRAL REACTIONS ==========
    db['CH_H_C_H2_cm3_7_3'] = RateConstant(1.20e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['C_CH_C2_H_cm3_7_4'] = RateConstant(1.20e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_CH3_C2H4_cm3_7_5'] = RateConstant(8.00e-11, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['C2_H_CH_C_cm3_7_6'] = RateConstant(9.60e-11, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_CH2_C2H2_H_cm3_7_7'] = RateConstant(1.20e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['C_CH3_C2_H2_H_cm3_7_8'] = RateConstant(1.20e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C_C2_H_cm3_7_9'] = RateConstant(1.20e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_CH3_C2H3_H_cm3_7_10'] = RateConstant(8.00e-11, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_Ar_Ar_C_H_cm3_7_11'] = RateConstant(1.00e-15, 1.00e-15, 1.00e-15, "UMIST (2012)", flag="7.11 Needs verification")
    db['C_H_CH_cm3_7_12'] = RateConstant(8.00e-11, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH2_CH2_C2H4_cm3_7_13'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH3_CH2_C2H5_cm3_7_14'] = RateConstant(8.00e-11, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH2_CH2_C2H2_H2_cm3_7_15'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH3_CH_C2H2_H2_cm3_7_16'] = RateConstant(1.20e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH2_C_C2H2_cm3_7_17'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C2H4_C2H2_CH3_cm3_7_18'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['C2H2_C_C2_CH2_cm3_7_19'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_CH4_C2H4_H_cm3_7_20'] = RateConstant(1.20e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH_H_CH2_cm3_7_21'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C2H2_C3H2_H_cm3_7_22'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_CH3_C2H2_H2_cm3_7_23'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C_C2_H2_cm3_7_24'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['H_CH4_CH3_H2_cm3_7_25'] = RateConstant(6.00e-12, 4.00e-12, 8.00e-12, "Baulch et al. (2005)")
    db['CH2_CH_C2_H2_H_cm3_7_26'] = RateConstant(1.20e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C2H2_C3H_H2_cm3_7_27'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C3H_C4H2_H_cm3_7_28'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C2H2_C2H_CH2_cm3_7_29'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_H2_CH2_H_cm3_7_30'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH_C2H3_C3H3_H_cm3_7_31'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C2H4_C3H4_H_cm3_7_32'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C2_C3_H_cm3_7_33'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C2H5_C3H5_H_cm3_7_34'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C3H2_C4H_H2_cm3_7_35'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH3_H_CH2_H2_cm3_7_36'] = RateConstant(6.00e-12, 4.00e-12, 8.00e-12, "Baulch et al. (2005)")
    db['CH_C2H6_C3H6_H_cm3_7_37'] = RateConstant(1.20e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH3_CH3_CH2_CH4_cm3_7_38'] = RateConstant(1.20e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH_CH4_CH2_CH3_cm3_7_39'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH3_CH3_C2H6_cm3_7_40'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH_C2H5_C3H6_cm3_7_41'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH2_CH2_C2H2_H2_cm3_7_42'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH_C_C2_H_cm3_7_43'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_CH_C2_H2_cm3_7_44'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C2H6_C3H6_H_cm3_7_45'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH_C2H4_C2H2_CH3_cm3_7_46'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['C2H_H_C2_H2_cm3_7_47'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_CH2_C2H2_H_cm3_7_48'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH3_CH3_C2H2_H2_H2_cm3_7_49'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['C2H2_H_C2_H2_H_cm3_7_50'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH_H2_CH2_H_cm3_7_51'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['C2_CH_C3_H_cm3_7_52'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH2_C2H3_C2H2_CH3_cm3_7_53'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['C_C2H3_C2_CH3_cm3_7_54'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C2H5_C3H6_cm3_7_55'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['C2H2_CH_C3_H2_H_cm3_7_56'] = RateConstant(1.00e-10, 8.00e-11, 1.20e-10, "Baulch et al. (2005)")
    db['CH_C2H6_C2H2_CH3_H_cm3_7_57'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH2_CH2_C2_H2_H2_cm3_7_58'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH_CH4_C2H4_H_cm3_7_59'] = RateConstant(6.00e-12, 4.00e-12, 6.00e-12, "Baulch et al. (2005)")
    db['CH_C2H4_C3H4_H_cm3_7_60'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH3_C2H5_C2H2_CH3_H2_cm3_7_61'] = RateConstant(1.00e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH2_CH3_C2H2_H_H2_cm3_7_62'] = RateConstant(1.20e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")
    db['CH2_C2H5_C2H2_CH3_H_cm3_7_63'] = RateConstant(1.20e-11, 8.00e-12, 1.20e-11, "Baulch et al. (2005)")

    # ========== TERMOLECULAR REACTIONS ==========
    db['H_H_M_H2_M_cm6_8_1'] = RateConstant(8.00e-33, 8.00e-33, 1.20e-32, "Baulch et al. (2005)")
    db['CH3_CH3_M_C2H6_M_cm6_8_2'] = RateConstant(3.60e-29, 2.40e-29, 3.60e-29, "Baulch et al. (2005)", flag="8.2 Third-body efficiencies need validation")

    # ========== WALL STICKING REACTIONS ==========
    db['stick_H_9_1'] = RateConstant(3.89e+02, 3.89e+02, 3.89e+03, "Perrin (1991)")
    db['stick_CH3_9_2'] = RateConstant(3.51e+03, 1.20e+03, 5.82e+03, "Matsuda et al. (1990)")
    db['stick_CH_9_3'] = RateConstant(6.25e+03, 1.25e+03, 6.25e+03, "Jauberteau et al. (1998)")
    db['stick_ArPlus_9_4'] = RateConstant(7.14e+03, 3.57e+03, 7.14e+03, "Boeuf (1987)")
    db['stick_ArStar_9_5'] = RateConstant(3.57e+02, 7.14e+01, 7.14e+02, "Phelps (1999)")
    db['stick_CH3Plus_9_6'] = RateConstant(1.16e+04, 5.82e+03, 1.16e+04, "Boeuf (1987)")
    db['stick_CH5Plus_9_7'] = RateConstant(1.09e+04, 5.47e+03, 1.09e+04, "Boeuf (1987)")
    db['stick_ArHPlus_9_8'] = RateConstant(7.14e+03, 3.57e+03, 7.14e+03, "Boeuf (1987)")
    db['stick_C2_9_9'] = RateConstant(1.25e+03, 1.25e+03, 6.25e+03, "Jauberteau et al. (1998)")
    db['stick_C_9_10'] = RateConstant(6.25e+03, 1.25e+03, 6.25e+03, "Jauberteau et al. (1998)")
    db['stick_C2H2_9_11'] = RateConstant(5.00e+02, 5.00e+02, 2.00e+03, "Matsuda et al. (1990)")
    db['stick_C2H4_9_12'] = RateConstant(5.00e+02, 5.00e+02, 2.00e+03, "Matsuda et al. (1990)")
    db['stick_CH2_9_13'] = RateConstant(2.00e+03, 1.00e+03, 3.00e+03, "Jauberteau et al. (1998)")
    db['stick_C2H6_9_14'] = RateConstant(3.20e+02, 4.00e+02, 1.60e+03, "Matsuda et al. (1990)")
    db['stick_CH3Minus_9_15'] = RateConstant(1.50e+03, 2.50e+03, 7.50e+03, "Boeuf (1987)")
    db['stick_H2_9_16'] = RateConstant(6.25e+02, 2.50e+02, 1.00e+03, "Perrin (1991)")
    db['stick_C2H5_9_17'] = RateConstant(2.00e+03, 5.00e+02, 2.00e+03, "Matsuda et al. (1990)")
    db['stick_HMinus_9_18'] = RateConstant(5.00e+03, 2.50e+03, 7.50e+03, "Boeuf (1987)")
    db['stick_C3H_9_19'] = RateConstant(1.25e+03, 5.00e+02, 2.00e+03, "Matsuda et al. (1990)")
    db['stick_C4H2_9_20'] = RateConstant(1.25e+03, 5.00e+02, 2.00e+03, "Matsuda et al. (1990)")
    db['stick_C3H3_9_21'] = RateConstant(1.25e+03, 5.00e+02, 2.00e+03, "Matsuda et al. (1990)")
    db['stick_C3H4_9_22'] = RateConstant(1.25e+03, 5.00e+02, 2.00e+03, "Matsuda et al. (1990)")
    db['stick_C3_9_23'] = RateConstant(1.25e+03, 5.00e+02, 2.00e+03, "Matsuda et al. (1990)")
    db['stick_C4H_9_24'] = RateConstant(1.25e+03, 5.00e+02, 2.00e+03, "Matsuda et al. (1990)")
    db['stick_C3H6_9_25'] = RateConstant(1.00e+03, 5.00e+02, 2.00e+03, "Matsuda et al. (1990)")
    db['stick_H3Plus_9_26'] = RateConstant(5.00e+03, 2.50e+03, 7.50e+03, "Boeuf (1987)")
    db['stick_CHPlus_9_27'] = RateConstant(5.00e+03, 2.50e+03, 7.50e+03, "Boeuf (1987)")
    db['stick_C2HPlus_9_28'] = RateConstant(5.00e+03, 2.50e+03, 7.50e+03, "Boeuf (1987)")

    # ========== WALL LOSS REACTIONS ==========
    db['loss_CH2_11_1'] = RateConstant(3.63e+03, 1.21e+03, 6.04e+03, "Jauberteau et al. (1998)")
    db['loss_H2_11_2'] = RateConstant(3.50e+03, 2.00e+03, 5.00e+03, "Matsuda (2004)")
    db['loss_C2_11_3'] = RateConstant(8.00e+02, 1.00e-04, 2.00e+03, "Alman & Ruzic (2003)")
    db['loss_e_11_4'] = RateConstant(7.50e+03, 5.00e+03, 1.00e+04, "Lieberman & Lichtenberg (2005)")
    db['loss_C2H6_11_5'] = RateConstant(5.00e+02, 5.00e+02, 1.50e+03, "Alman & Ruzic (2003)")
    db['loss_CH4_11_6'] = RateConstant(1.50e+03, 1.00e+03, 2.00e+03, "Alman & Ruzic (2003)")
    db['loss_Ar_11_7'] = RateConstant(2.00e+03, 1.00e+03, 3.00e+03, "Estimated", flag="11.7 Needs validation")
    db['loss_C_11_8'] = RateConstant(2.00e+03, 5.00e+02, 2.00e+03, "Alman & Ruzic (2003)")
    db['loss_CH_11_9'] = RateConstant(1.00e+04, 1.00e+03, 1.00e+04, "Estimated", flag="11.9 Needs validation")
    db['loss_C3H_11_10'] = RateConstant(5.50e+02, 1.00e+02, 1.00e+03, "Alman & Ruzic (2003)")
    db['loss_C4H2_11_11'] = RateConstant(5.50e+02, 1.00e+02, 1.00e+03, "Alman & Ruzic (2003)")
    db['loss_C2H_11_12'] = RateConstant(5.50e+02, 1.00e+02, 1.00e+03, "Alman & Ruzic (2003)")
    db['loss_C3H3_11_13'] = RateConstant(5.50e+02, 1.00e+02, 1.00e+03, "Alman & Ruzic (2003)")
    db['loss_C3_11_14'] = RateConstant(5.50e+02, 1.00e+02, 1.00e+03, "Alman & Ruzic (2003)")
    db['loss_C4H_11_15'] = RateConstant(5.50e+02, 1.00e+02, 1.00e+03, "Alman & Ruzic (2003)")
    db['loss_C3H6_11_16'] = RateConstant(1.00e+02, 1.00e+02, 1.00e+03, "Alman & Ruzic (2003)")
    db['loss_C2H5_11_17'] = RateConstant(1.50e+03, 5.00e+02, 1.50e+03, "Alman & Ruzic (2003)")
    db['loss_C3H4_11_18'] = RateConstant(8.25e+02, 1.00e+02, 1.50e+03, "Alman & Ruzic (2003)")
    db['loss_C2H2_11_19'] = RateConstant(1.00e+03, 1.00e+03, 2.00e+03, "Alman & Ruzic (2003)")
    db['loss_C2H4_11_20'] = RateConstant(1.50e+03, 1.50e+03, 1.50e+03, "Alman & Ruzic")
    db['loss_CH3_11_21'] = RateConstant(1.20e+03, 1.20e+03, 1.20e+03, "Jauberteau et al.")
    db['loss_C2H3_11_22'] = RateConstant(1.50e+03, 1.50e+03, 1.50e+03, "Alman & Ruzic")
    db['loss_C3H2_11_23'] = RateConstant(1.00e+03, 1.00e+03, 1.00e+03, "Alman & Ruzic")
    db['loss_C3H5_11_24'] = RateConstant(1.00e+03, 1.00e+03, 1.00e+03, "Alman & Ruzic")
    db['loss_C2H2Star_11_25'] = RateConstant(1.00e+03, 1.00e+02, 1.00e+04, "Estimated")

    return db


def get_rates_by_species(db, species_name):
    """Get all rates involving a specific species."""
    return {k: v for k, v in db.items() if species_name in k}


def get_rates_by_source(db, source_name):
    """Get all rates from a specific literature source."""
    return {k: v for k, v in db.items() if source_name.lower() in v.source.lower()}


def get_flagged_rates(db):
    """Get all rates that have validation flags."""
    return {k: v for k, v in db.items() if v.flag}


def get_tunable_rates_for_target(target_species):
    """
    Get rates most relevant for tuning specific target species.

    Parameters:
    -----------
    target_species : str
        'H', 'CH', or 'C2'

    Returns:
    --------
    dict : {rate_name: importance_note}
    """

    if target_species == 'H':
        return {
            'e_CH4_CH3_H_cm3_1_1': 'e + CH4 → CH3 + H (primary H source)',
            'e_H2_H_H_cm3_1_4': 'e + H2 → 2H (H2 dissociation)',
            'e_CH4_CH_H2_H_vib_cm3_1_3': 'e + CH4 → CH + H2 + H',
            'e_CH4_CH_H_H2_cm3_1_11': 'e + CH4 → CH + H + H2',
            'ArStar_CH4_CH3_H_cm3_3_1': 'Ar* + CH4 → CH3 + H + Ar',
            'H_CH4_CH3_H2_cm3_7_25': 'H + CH4 → CH3 + H2 (H consumption)',
            'CH2_H_CH_H2_cm3_7_1': 'CH2 + H → CH + H2 (H consumption)',
            'CH_H_C_H2_cm3_7_3': 'CH + H → C + H2 (H consumption)',
        }

    elif target_species == 'CH':
        return {
            'e_CH4_CH_H2_H_vib_cm3_1_3': 'e + CH4 → CH + H2 + H (primary CH source)',
            'e_CH4_CH_H_H2_cm3_1_11': 'e + CH4 → CH + H + H2 (secondary CH source)',
            'ArStar_CH2_CH_H_cm3_3_4': 'Ar* + CH2 → CH + H + Ar (CH production)',
            'CH2_H_CH_H2_cm3_7_1': 'CH2 + H → CH + H2 (CH production)',
            'CH_CH_C2_H2_cm3_5_4': 'CH + CH → C2 + H2 (CH loss, C2 gain) ***KEY***',
            'CH_CH3_C2H2_H2_cm3_7_23': 'CH + CH3 → C2H2 + H2 (CH loss)',
            'CH_CH3_C2H3_H_cm3_7_10': 'CH + CH3 → C2H3 + H (CH loss)',
            'CH_H_C_H2_cm3_7_3': 'CH + H → C + H2 (CH loss)',
            'loss_CH_11_9': 'CH wall loss ***CRITICAL - 10x range***',
            'stick_CH_9_3': 'CH wall sticking (5x range)',
        }

    elif target_species == 'C2':
        return {
            'CH_CH_C2_H2_cm3_5_4': 'CH + CH → C2 + H2 ***PRIMARY C2 SOURCE***',
            'C_CH_C2_H_cm3_7_4': 'C + CH → C2 + H (C2 formation)',
            'CH_C_C2_H_cm3_7_9': 'CH + C → C2 + H (C2 formation)',
            'e_C2H2_C2_H2_cm3_1_16': 'e + C2H2 → C2 + H2 (C2 from C2H2)',
            'C2H2_C_C2_CH2_cm3_7_19': 'C2H2 + C → C2 + CH2 (C2 formation)',
            'C2_H_CH_C_cm3_7_6': 'C2 + H → CH + C (C2 loss)',
            'loss_C2_11_3': 'C2 wall loss ***HUGE range 1e-4 to 2e3!***',
            'stick_C2_9_9': 'C2 wall sticking (5x range)',
        }

    else:
        return {}


if __name__ == '__main__':
    # Test the database
    db = get_complete_rate_database()
    print(f"Complete database loaded: {len(db)} rates")

    # Show flagged rates
    flagged = get_flagged_rates(db)
    print(f"\nFlagged rates needing validation: {len(flagged)}")
    for name, rate in list(flagged.items())[:5]:
        print(f"  {name}: {rate.flag}")

    # Show CH-related rates
    ch_rates = get_rates_by_species(db, 'CH')
    print(f"\nRates involving CH: {len(ch_rates)}")
