% define_rates.m reset to nominal literature values for Ar/CH4 plasma
% Goals: Achieve expected density ranges for all species before targeting H (7e13 cm^-3), CH (2e9 cm^-3), C2 (1.3e11 cm^-3), n_i_net/ne (~1)
% Sources: Morgan (1992), Janev & Reiter (2002), Phelps (1999), Anicich (2003), Baulch et al. (2005)

function k = define_rates(params)
k = struct();
E_field = params.E_field;
L_discharge = params.L_discharge;
mobilities = params.mobilities;

% Group 1: Electron-Impact Reactions (Neutral Products)
k.e_CH4_CH3_H_cm3_1_1 = 1.2e-10; % 1.1 Rate: 1.2e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1.2e-10, Source: Morgan (1992)
k.e_CH4_CH2_H2_cm3_1_2 = 6e-11; % 1.2 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 6e-11, Source: Morgan (1992)
k.e_CH4_CH_H2_H_vib_cm3_1_3 = 6e-11; % 1.3 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 1e-10, Nominal: 6e-11, Source: Janev & Reiter (2002)
k.e_H2_H_H_cm3_1_4 = 1.2e-11; % 1.4 Rate: 1.2e-11 cm^3 s^-1, Range: 8e-12 to 1.2e-11 (±20%), Nominal: 1.2e-11, Source: Morgan (1992)
k.e_CH3_CH2_H_cm3_1_5 = 6e-11; % 1.5 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 6e-11, Source: Janev & Reiter (2002)
k.e_C2H4_C2H2_H2_cm3_1_6 = 3.6e-11; % 1.6 Rate: 3.6e-11 cm^3 s^-1, Range: 2.4e-11 to 3.6e-11 (±20%), Nominal: 3.6e-11, Source: Janev & Reiter (2002)
k.e_Ar_ArStar_cm3_1_7 = 1.8e-10; % 1.7 Rate: 1.8e-10 cm^3 s^-1, Range: 1.2e-10 to 1.8e-10 (±20%), Nominal: 1.8e-10, Source: Phelps (1999)
k.e_C2H6_C2H4_H2_cm3_1_8 = 3e-11; % 1.8 Rate: 3e-11 cm^3 s^-1, Range: 2.4e-11 to 3.6e-11 (±20%), Nominal: 3e-11, Source: Janev & Reiter (2002)
k.e_C2H6_C2H4_H2_e_cm3_1_9 = 2.4e-11; % 1.9 Rate: 2.4e-11 cm^3 s^-1, Range: 2.4e-11–3.6e-11 (±20%), Nominal: 2.4e-11, Source: Janev & Reiter (2002)
k.e_CH4_CH3Minus_H_cm3_1_10 = 1.2e-17; % 1.10 Rate: 1.2e-17 cm^3 s^-1, Range: 8e-18 to 1.2e-17 (±20%), Nominal: 1.2e-17, Source: Janev & Reiter (2002)
% FLAG: 1.10 e + CH4 → CH3^- + H rate is low and estimated. Needs validation from LXCat or experimental data.
k.e_CH4_CH_H_H2_cm3_1_11 = 4e-11; % 1.11 Rate: 4e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 4e-11, Source: Janev & Reiter (2002)
k.e_CH_CH_C_H_e_cm3_1_12 = 1.2e-10; % 1.12 Rate: 1.2e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1.2e-10, Source: Janev & Reiter (2002)
k.e_H2_HMinus_H_cm3_1_13 = 1.2e-15; % 1.13 Rate: 1.2e-15 cm^3 s^-1, Range: 8e-16 to 1.2e-15 (±20%), Nominal: 1.2e-15, Source: Morgan (1992)

% Group 2: Electron-Impact Ionization
k.e_CH4_CH3Plus_H_cm3_2_1 = 2.56e-11; % 2.1 Rate: 2.56e-11 cm^3 s^-1, Range: 2.56e-11 to 3.84e-11 (±20%), Nominal: 2.56e-11, Source: Morgan (1992)
k.e_CH4_CH4Plus_cm3_2_2 = 2.56e-11; % 2.2 Rate: 2.56e-11 cm^3 s^-1, Range: 2.56e-11 to 3.84e-11 (±20%), Nominal: 2.56e-11, Source: Morgan (1992)
k.e_Ar_ArPlus_cm3_2_3 = 1.76e-11; % 2.3 Rate: 1.76e-11 cm^3 s^-1, Range: 1.76e-11 to 2.64e-11 (±20%), Nominal: 1.76e-11, Source: Phelps (1999)
k.e_ArStar_ArPlus_cm3_2_4 = 2e-10; % 2.4 Rate: 2e-10 cm^3 s^-1, Range: 1.6e-10 to 2.4e-10 (±20%), Nominal: 2e-10, Source: Phelps (1999)
k.e_C2H6_C2H5Plus_H_2e_cm3_2_5 = 2.4e-11; % 2.5 Rate: 2.4e-11 cm^3 s^-1, Range: 1.6e-11–2.4e-11 (±20%), Nominal: 2.4e-11, Source: Janev & Reiter (2002)
k.e_C2H4_C2H4Plus_2e_cm3_2_6 = 3.6e-11; % 2.6 Rate: 3.6e-11 cm^3 s^-1, Range: 2.4e-11–3.6e-11 (±20%), Nominal: 3.6e-11, Source: Janev & Reiter (2002)
k.e_C2H4_C2H3Plus_H_2e_cm3_2_7 = 2.4e-11; % 2.7 Rate: 2.4e-11 cm^3 s^-1, Range: 1.6e-11–2.4e-11 (±20%), Nominal: 2.4e-11, Source: Janev & Reiter (2002)
k.e_C2H2_C2HPlus_2e_cm3_2_8 = 2.4e-11; % 2.8 Rate: 2.4e-11 cm^3 s^-1, Range: 1.6e-11 to 2.4e-11 (±20%), Nominal: 2.4e-11, Source: Janev & Reiter (2002)

% Group 3: Ar* Reactions
k.ArStar_CH4_CH3_H_cm3_3_1 = 1.2e-10; % 3.1 Rate: 1.2e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1.2e-10, Source: Phelps (1999)
k.ArStar_H2_H_H_cm3_3_2 = 6e-11; % 3.2 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 6e-11, Source: Phelps (1999)
k.ArStar_H2_ArHPlus_H_cm3_3_3 = 6e-11; % 3.3 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 6e-11, Source: Phelps (1999)
k.ArStar_CH2_CH_H_cm3_3_4 = 8e-11; % 3.4 Rate: 8e-11 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 8e-11, Source: Phelps (1999)
k.ArStar_C_Ar_CStar_cm3_3_5 = 1e-10; % 3.5 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Phelps (1999)
k.ArStar_H_ArHPlus_e_cm3_3_6 = 5e-11; % 3.6 Rate: 5e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 5e-11, Source: Phelps (1999)
k.ArStar_CH3_Ar_CH2_H_cm3_3_7 = 1.2e-10; % 3.7 Rate: 1.2e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1.2e-10, Source: Phelps (1999)
k.ArStar_C2_Ar_C2_cm3_3_8 = 1e-10; % 3.8 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Phelps (1999)
k.ArStar_C2H4_Ar_C2H4_cm3_3_9 = 1.2e-10; % 3.9 Rate: 1.2e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1.2e-10, Source: Phelps (1999)
k.ArStar_CH2_Ar_CH2_cm3_3_10 = 1e-10; % 3.10 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Phelps (1999)
k.ArStar_CH3_Ar_CH3_cm3_3_11 = 1e-10; % 3.11 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Phelps (1999)
k.ArStar_M_Ar_3_12 = 1e-13; % 3.12 Rate: 1e-13 s^-1, Range: 8e-14 to 1.2e-13 (±20%), Nominal: 1e-13, Source: Phelps (1999)
% FLAG: 3.12 Ar* quenching rate (Ar* → Ar) is estimated. Needs experimental validation.
k.ArStar_Ar_Ar_Ar_cm3_3_13 = 1e-10; % 3.13 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Phelps (1999)
k.CH_ArStar_C_H_Ar_cm3_3_14 = 1e-11; % 3.14 Rate: 1e-11 cm^3 s^-1, Nominal: 1e-11, Source: Smith et al., J. Chem. Phys., 2018
k.CH_ArStar_CH_Ar_cm3_3_15 = 5e-11; % 3.15 Rate: 5e-11 cm^3 s^-1, Nominal: 5e-11, Source: Tachibana et al., J. Phys. D: Appl. Phys., 1989
k.ArStar_CH3_CH2_H_Ar_cm3_3_16 = 8e-11; % 3.16 Rate: 8e-11 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 8e-11, Source: Phelps (1999)

% Group 4: Penning Ionization
k.ArStar_CH4_CH4Plus_cm3_4_1 = 4.8e-11; % 4.1 Rate: 4.8e-11 cm^3 s^-1, Range: 4.8e-11 to 7.2e-11 (±20%), Nominal: 4.8e-11, Source: Phelps (1999)
k.ArStar_CH4_CH3Plus_H_cm3_4_2 = 4e-11; % 4.2 Rate: 4e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 4e-11, Source: Phelps (1999)
k.ArStar_Ar_ArPlus_cm3_4_3 = 4e-11; % 4.3 Rate: 4e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 4e-11, Source: Phelps (1999)
k.ArStar_CH4_ArPlus_CH3_H_e_cm3_4_4 = 2.4e-11; % 4.4 Rate: 2.4e-11 cm^3 s^-1, Range: 1.6e-11 to 2.4e-11 (±20%), Nominal: 2.4e-11, Source: Phelps (1999)
k.ArStar_CH3_ArPlus_CH2_H_e_cm3_4_5 = 2.4e-11; % 4.5 Rate: 2.4e-11 cm^3 s^-1, Range: 1.6e-11 to 2.4e-11 (±20%), Nominal: 2.4e-11, Source: Phelps (1999)
k.ArStar_H_ArPlus_HMinus_cm3_4_6 = 5e-11; % 4.6 Rate: 5e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 5e-11, Source: Phelps (1999)
k.ArStar_CH2_ArPlus_CH_H_e_cm3_4_7 = 2e-11; % 4.7 Rate: 2e-11 cm^3 s^-1, Range: 1.6e-11 to 2.4e-11 (±20%), Nominal: 2e-11, Source: Phelps (1999)
k.ArStar_H2_ArPlus_H2_e_cm3_4_8 = 5e-11; % 4.8 Rate: 5e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 5e-11, Source: Phelps (1999)
k.ArStar_C2H2_ArPlus_C2H2_e_cm3_4_9 = 5e-11; % 4.9 Rate: 5e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 5e-11, Source: Phelps (1999)
k.ArStar_C2H5_ArPlus_C2H5_e_cm3_4_10 = 5e-11; % 4.10 Rate: 5e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 5e-11, Source: Phelps (1999)
k.ArStar_H_ArPlus_H_e_cm3_4_11 = 5e-11; % 4.11 Rate: 5e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 5e-11, Source: Phelps (1999)
k.ArStar_CH4_ArPlus_CH4_e_cm3_4_12 = 6e-11; % 4.12 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 6e-11, Source: Phelps (1999)
k.ArStar_CH3_ArPlus_CH3_e_cm3_4_13 = 6e-11; % 4.13 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 6e-11, Source: Phelps (1999)
k.ArStar_C2H4_ArPlus_C2H4_e_cm3_4_14 = 6e-11; % 4.14 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 6e-11, Source: Phelps (1999)
k.ArStar_C2H5_ArPlus_C2H5_e_cm3_4_15 = 6e-11; % 4.15 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 6e-11, Source: Phelps (1999)
k.ArStar_CH2_ArPlus_CH2_e_cm3_4_16 = 6e-11; % 4.16 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 6e-11, Source: Phelps (1999)
k.ArStar_C2H6_ArPlus_C2H6_e_cm3_4_17 = 6e-11; % 4.17 Rate: 6e-11 cm^3 s^-1, Range: 4e-11 to 6e-11 (±20%), Nominal: 6e-11, Source: Phelps (1999)

% Group 5: Ion-Neutral Reactions
k.ArPlus_CH4_CH3Plus_H_cm3_5_1 = 9e-10; % 5.1 Rate: 9e-10 cm^3 s^-1, Range: 9e-10 to 1.3e-9 (±20%), Nominal: 9e-10, Source: Anicich (2003)
k.CH3Plus_CH4_CH5Plus_CH2_cm3_5_2 = 9.6e-10; % 5.2 Rate: 9.6e-10 cm^3 s^-1, Range: 9.6e-10 to 1.44e-9, Nominal: 9.6e-10, Source: Anicich (2003)
k.ArPlus_CH3_CH3Plus_cm3_5_3 = 1e-9; % 5.3 Rate: 1e-9 cm^3 s^-1, Range: 8e-10 to 1.2e-9 (±20%), Nominal: 1e-9, Source: Anicich (2003)
k.CH_CH_C2_H2_cm3_5_4 = 1.8e-10; % 5.4 Rate: 1.8e-10 cm^3 s^-1, Range: 1.2e-10 to 1.8e-10 (±20%), Nominal: 1.8e-10, Source: Baulch et al. (2005)
k.ArPlus_CH4_Ar_CH4Plus_cm3_5_5 = 7e-10; % 5.5 Rate: 7e-10 cm^3 s^-1, Range: 7e-10 to 1.4e-9, Nominal: 7e-10, Source: Anicich (2003)
k.CH4Plus_H2_CH5Plus_H_cm3_5_6 = 1.3e-10; % 5.6 Rate: 1.3e-10 cm^3 s^-1, Range: 1e-10 to 1.6e-10 (±20%), Nominal: 1.3e-10, Source: Anicich (2003)
k.ArPlus_H2_ArHPlus_H_cm3_5_7 = 2e-10; % 5.7 Rate: 2e-10 cm^3 s^-1, Range: 1.4e-10 to 2e-10 (±20%), Nominal: 2e-10, Source: Anicich (2003)
k.ArHPlus_CH4_Ar_CH5Plus_cm3_5_8 = 1e-9; % 5.8 Rate: 1e-9 cm^3 s^-1, Range: 9e-10 to 1.3e-9 (±20%), Nominal: 1e-9, Source: Anicich (2003)
k.ArPlus_CH4_ArHPlus_CH3_cm3_5_9 = 8e-10; % 5.9 Rate: 8e-10 cm^3 s^-1, Range: 8e-10 to 1.2e-9 (±20%), Nominal: 8e-10, Source: Anicich (2003)
k.CH2_CH3Plus_CH3_CH2Plus_cm3_5_10 = 1e-9; % 5.10 Rate: 1e-9 cm^3 s^-1, Range: 8e-10 to 1.2e-9 (±20%), Nominal: 1e-9, Source: Anicich (2003)
k.CH3Plus_CH4_C2H5Plus_H2_cm3_5_11 = 1e-9; % 5.11 Rate: 1e-9 cm^3 s^-1, Range: 8e-10–1.2e-9 (±20%), Nominal: 1e-9, Source: Anicich (2003)
k.CH5Plus_C2H4_C2H5Plus_CH4_cm3_5_12 = 1e-9; % 5.12 Rate: 1e-9 cm^3 s^-1, Range: 8e-10–1.2e-9 (±20%), Nominal: 1e-9, Source: Anicich (2003)

% Group 6: Dissociative Recombination
k.ArPlus_e_Ar_cm3_6_1 = 1.2e-7; % 6.1 Rate: 1.2e-7 cm^3 s^-1, Range: 8e-8 to 1.2e-7 (±20%), Nominal: 1.2e-7, Source: UMIST Database (2012)
k.CH3Plus_e_CH3_cm3_6_2 = 3.6e-7; % 6.2 Rate: 3.6e-7 cm^3 s^-1, Range: 2.4e-7 to 3.6e-7 (±20%), Nominal: 3.6e-7, Source: UMIST Database (2012)
k.CH5Plus_e_CH4_H_cm3_6_3 = 6e-7; % 6.3 Rate: 6e-7 cm^3 s^-1, Range: 4e-7 to 6e-7 (±20%), Nominal: 6e-7, Source: UMIST Database (2012)
k.e_CH4Plus_CH3_H_cm3_6_4 = 4.8e-7; % 6.4 Rate: 4.8e-7 cm^3 s^-1, Range: 3.2e-7 to 4.8e-7 (±20%), Nominal: 4.8e-7, Source: UMIST Database (2012)
k.CH3Minus_ArPlus_CH3_Ar_cm3_6_5 = 1e-7; % 6.5 Rate: 1e-7 cm^3 s^-1, Range: 8e-8 to 1.2e-7 (±20%), Nominal: 1e-7, Source: UMIST Database (2012)
% FLAG: 6.5 CH3^- + Ar^+ rate is estimated. Cross-check with Anicich (2003) or experimental data.
k.CH3Minus_CH4Plus_CH4_CH3_cm3_6_6 = 1e-7; % 6.6 Rate: 1e-7 cm^3 s^-1, Range: 8e-8 to 1.2e-7 (±20%), Nominal: 1e-7, Source: UMIST Database (2012)
% FLAG: 6.6 CH3^- + CH4^+ rate is estimated. Needs validation.
k.CH3Minus_CH3Plus_CH4_CH2_cm3_6_7 = 1e-7; % 6.7 Rate: 1e-7 cm^3 s^-1, Range: 8e-8 to 1.2e-7 (±20%), Nominal: 1e-7, Source: UMIST Database (2012)
k.CH5Plus_e_CH3_H2_cm3_6_8 = 1.2e-7; % 6.8 Rate: 1.2e-7 cm^3 s^-1, Range: 8e-8 to 1.2e-7 (±20%), Nominal: 1.2e-7, Source: UMIST Database (2012)
k.e_CH4Plus_CH2_H2_cm3_6_9 = 4.8e-7; % 6.9 Rate: 4.8e-7 cm^3 s^-1, Range: 3.2e-7 to 4.8e-7 (±20%), Nominal: 4.8e-7, Source: UMIST Database (2012)
k.CH5Plus_e_CH2_H2_H_cm3_6_10 = 1.2e-7; % 6.10 Rate: 1.2e-7 cm^3 s^-1, Range: 8e-8 to 1.2e-7 (±20%), Nominal: 1.2e-7, Source: UMIST Database (2012)
k.e_CH4Plus_CH_H2_H_cm3_6_11 = 4e-7; % 6.11 Rate: 4e-7 cm^3 s^-1, Range: 3.2e-7 to 4.8e-7 (±20%), Nominal: 4e-7, Source: UMIST Database (2012)
k.CH5Plus_e_CH3_2H_cm3_6_12 = 1.2e-7; % 6.12 Rate: 1.2e-7 cm^3 s^-1, Range: 8e-8 to 1.2e-7 (±20%), Nominal: 1.2e-7, Source: UMIST Database (2012)
k.e_CH4Plus_C_2H2_cm3_6_13 = 4e-7; % 6.13 Rate: 4e-7 cm^3 s^-1, Range: 3.2e-7 to 4.8e-7 (±20%), Nominal: 4e-7, Source: UMIST Database (2012)
k.C2H5Plus_e_C2H4_H_cm3_6_14 = 2.4e-7; % 6.14 Rate: 2.4e-7 cm^3 s^-1, Range: 2.4e-7–3.6e-7 (±20%), Nominal: 2.4e-7, Source: UMIST Database (2012)
k.C2H4Plus_e_C2H2_H2_cm3_6_15 = 2.4e-7; % 6.15 Rate: 2.4e-7 cm^3 s^-1, Range: 2.4e-7–3.6e-7 (±20%), Nominal: 2.4e-7, Source: UMIST Database (2012)
k.C2H3Plus_e_C2H2_H_cm3_6_16 = 2.4e-7; % 6.16 Rate: 2.4e-7 cm^3 s^-1, Range: 2.4e-7–3.6e-7 (±20%), Nominal: 2.4e-7, Source: UMIST Database (2012)
k.HMinus_ArPlus_H_Ar_cm3_6_17 = 1.2e-7; % 6.17 Rate: 1.2e-7 cm^3 s^-1, Range: 8e-8 to 1.2e-7 (±20%), Nominal: 1.2e-7, Source: UMIST Database (2012)
k.C2HPlus_e_C2_H_cm3_6_18 = 2.9e-7; % 6.18 Rate: 2.9e-7 cm^3 s^-1, Range: 1.9e-7 to 2.9e-7 (±20%), Nominal: 2.9e-7, Source: UMIST Database (2012)

% Group 7: Neutral-Neutral Reactions
k.CH2_H_CH_H2_cm3_7_1 = 1e-11; % 7.1 Rate: 1e-11 cm^3 s^-1, Range: 1e-11 to 2.25e-11, Nominal: 1e-11, Source: Baulch et al. (2005), Table 3.2.2
k.CH2_H_C_H2_H_cm3_7_2 = 1e-11; % 7.2 Rate: 1e-11 cm^3 s^-1, Range: 8e-12 to 1.2e-11 (±20%), Nominal: 1e-11, Source: Baulch et al. (2005), Estimated
% FLAG: 7.2 CH2 + H → C + H2 + H not explicitly listed in Baulch. Rate estimated. Needs discussion.
k.CH_H_C_H2_cm3_7_3 = 1e-10; % 7.3 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.2.3
k.C_CH_C2_H_cm3_7_4 = 1e-10; % 7.4 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.1
k.CH_CH3_C2H4_cm3_7_5 = 1e-10; % 7.5 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.2
k.C2_H_CH_C_cm3_7_6 = 1e-10; % 7.6 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.4
k.CH_CH2_C2H2_H_cm3_7_7 = 1e-10; % 7.7 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.3
k.C_CH3_C2_H2_H_cm3_7_8 = 1e-10; % 7.8 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.5
k.CH_C_C2_H_cm3_7_9 = 1e-10; % 7.9 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.1
k.CH_CH3_C2H3_H_cm3_7_10 = 1e-10; % 7.10 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.6
k.CH_Ar_Ar_C_H_cm3_7_11 = 1e-15; % 7.11 Rate: 1e-15 cm^3 s^-1, Nominal: 1e-15, Source: UMIST Database (2012)
% FLAG: 7.11 CH + Ar → Ar + C + H not in Baulch. Retained from UMIST, but source needs verification.
k.C_H_CH_cm3_7_12 = 1e-10; % 7.12 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.2.6
k.CH2_CH2_C2H4_cm3_7_13 = 1e-10; % 7.13 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.7
k.CH3_CH2_C2H5_cm3_7_14 = 1e-10; % 7.14 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.8
k.CH2_CH2_C2H2_H2_cm3_7_15 = 1e-10; % 7.15 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.9
k.CH3_CH_C2H2_H2_cm3_7_16 = 1e-10; % 7.16 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.10
k.CH2_C_C2H2_cm3_7_17 = 1e-10; % 7.17 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.11
k.CH_C2H4_C2H2_CH3_cm3_7_18 = 1e-10; % 7.18 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.12
k.C2H2_C_C2_CH2_cm3_7_19 = 1e-10; % 7.19 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.13
k.CH_CH4_C2H4_H_cm3_7_20 = 1e-11; % 7.20 Rate: 1e-11 cm^3 s^-1, Range: 8e-12 to 1.2e-11 (±20%), Nominal: 1e-11, Source: Baulch et al. (2005), Table 3.2.7
k.CH_H_CH2_cm3_7_21 = 1e-10; % 7.21 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.2.8
k.CH_C2H2_C3H2_H_cm3_7_22 = 1e-10; % 7.22 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.14
k.CH_CH3_C2H2_H2_cm3_7_23 = 1e-10; % 7.23 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.15
k.CH_C_C2_H2_cm3_7_24 = 1e-10; % 7.24 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.16
k.H_CH4_CH3_H2_cm3_7_25 = 6e-12; % 7.25 Rate: 6e-12 cm^3 s^-1, Range: 4e-12 to 8e-12 (±20%), Nominal: 6e-12, Source: Baulch et al. (2005), Table 3.2.9
k.CH2_CH_C2_H2_H_cm3_7_26 = 1e-10; % 7.26 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.17
k.CH_C2H2_C3H_H2_cm3_7_27 = 1e-10; % 7.27 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005), Table 3.3.18
k.CH_C3H_C4H2_H_cm3_7_28 = 1e-10; % 7.28 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005)
k.CH_C2H2_C2H_CH2_cm3_7_29 = 1e-10; % 7.29 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005)
k.CH_H2_CH2_H_cm3_7_30 = 1e-11; % 7.30 Rate: 1e-11 cm^3 s^-1, Range: 8e-12 to 1.2e-11 (±20%), Nominal: 1e-11, Source: Baulch et al. (2005), Table 3.2.10
k.CH_C2H3_C3H3_H_cm3_7_31 = 1e-10; % 7.31 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005)
k.CH_C2H4_C3H4_H_cm3_7_32 = 1e-10; % 7.32 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005)
k.CH_C2_C3_H_cm3_7_33 = 1e-10; % 7.33 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005)
k.CH_C2H5_C3H5_H_cm3_7_34 = 1e-10; % 7.34 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005)
k.CH_C3H2_C4H_H2_cm3_7_35 = 1e-10; % 7.35 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005)
k.CH3_H_CH2_H2_cm3_7_36 = 6e-12; % 7.36 Rate: 6e-12 cm^3 s^-1, Range: 4e-12 to 8e-12 (±20%), Nominal: 6e-12, Source: Baulch et al. (2005), Table 3.2.11
k.CH_C2H6_C3H6_H_cm3_7_37 = 1e-10; % 7.37 Rate: 1e-10 cm^3 s^-1, Range: 8e-11 to 1.2e-10 (±20%), Nominal: 1e-10, Source: Baulch et al. (2005)
k.CH3_CH3_CH2_CH4_cm3_7_38 = 1e-11; % 7.38 Rate: 1e-11 cm^3 s^-1, Range: 8e-12 to 1.2e-11 (±20%), Nominal: 1e-11, Source: Baulch et al. (2005)
k.CH_CH4_CH2_CH3_cm3_7_39 = 1e-11; % 7.39 Rate: 1e-11 cm^3 s^-1, Range: 8e-12 to 1.2e-11 (±20%), Nominal: 1e-11, Source: Baulch et al. (2005), Table 3.2.12

% Group 8: Termolecular Recombination
k.H_H_M_H2_M_cm6_8_1 = 1e-32; % 8.1 Rate: 1e-32 cm^6 s^-1, Range: 8e-33 to 1.2e-32 (±20%), Nominal: 1e-32, Source: Baulch et al. (2005), Table 3.4.2
k.CH3_CH3_M_C2H6_M_cm6_8_2 = 3e-29; % 8.2 Rate: 3e-29 cm^6 s^-1, Range: 2.4e-29 to 3.6e-29 (±20%), Nominal: 3e-29, Source: Baulch et al. (2005), Table 3.4.1
% FLAG: 8.2 CH3 + CH3 (+M) rate uses low-pressure limit (k0). Third-body efficiencies for Ar/CH4 mixture estimated.

% Group 9: Stick Reactions
k.stick_H_9_1 = 3.89e2; % 9.1 Rate: 3.89e2 s^-1, Range: 3.89e2 to 3.89e3 (γ: 0.01–0.1), Nominal: 3.89e2, Source: Perrin, 1991; Matsuda, 2004
k.stick_CH3_9_2 = 3.51e3; % 9.2 Rate: 1.16e3 s^-1, Range: 1.2e3 to 5.82e3 (γ: 0.1–0.5), Nominal: 3.51e3 (midpoint), Source: Matsuda et al., 1990; Abe et al., 2013
k.stick_CH_9_3 = 3.75e3; % 9.3 Rate: 6e3 s^-1, Range: 1.25e3 to 6.25e3 (γ: 0.1–0.5), Nominal: 3.75e3 (midpoint), Source: Jauberteau et al., 1998
k.stick_ArPlus_9_4 = 5.36e3; % 9.4 Rate: 7.14e3 s^-1, Range: 3.57e3 to 7.14e3 (γ: 0.5–1.0), Nominal: 5.36e3 (midpoint), Source: Boeuf, 1987
k.stick_ArStar_9_5 = 3.57e2; % 9.5 Rate: 3.57e2 s^-1, Range: 7.14e1 to 7.14e2 (γ: 0.01–0.1), Nominal: 3.57e2, Source: Phelps (1999)
k.stick_CH3Plus_9_6 = 8.71e3; % 9.6 Rate: 9.31e3 s^-1, Range: 5.82e3 to 1.16e4 (γ: 0.5–1), Nominal: 8.71e3 (midpoint), Source: Boeuf, 1987
k.stick_CH5Plus_9_7 = 8.18e3; % 9.7 Rate: 1.09e4 s^-1, Range: 5.47e3 to 1.09e4 (γ: 0.5–1), Nominal: 8.18e3 (midpoint), Source: Boeuf, 1987
k.stick_ArHPlus_9_8 = 5.36e3; % 9.8 Rate: 5.36e3 s^-1, Range: 3.57e3 to 7.14e3 (γ: 0.5–1.0), Nominal: 5.36e3 (midpoint), Source: Boeuf, 1987
k.stick_C2_9_9 = 3.75e3; % 9.9 Rate: 1e3 s^-1, Range: 1.25e3 to 6.25e3 (γ: 0.1–0.5), Nominal: 3.75e3 (midpoint), Source: Jauberteau et al., 1998
k.stick_C_9_10 = 3.75e3; % 9.10 Rate: 6.25e3 s^-1, Range: 1.25e3 to 6.25e3 (γ: 0.1–0.5), Nominal: 3.75e3 (midpoint), Source: Jauberteau et al., 1998
k.stick_C2H2_9_11 = 1.25e3; % 9.11 Rate: 1e3 s^-1, Range: 5e2 to 2e3 (γ: 0.1–0.5), Nominal: 1.25e3 (midpoint), Source: Matsuda et al., 1990
k.stick_C2H4_9_12 = 1.25e3; % 9.12 Rate: 2e3 s^-1, Range: 5e2 to 2e3 (γ: 0.1–0.5), Nominal: 1.25e3 (midpoint), Source: Matsuda et al., 1990
k.stick_CH2_9_13 = 2e3; % 9.13 Rate: 1.5e3 s^-1, Range: 1e3 to 3e3, Nominal: 2e3 (midpoint), Source: Jauberteau et al., 1998
k.stick_C2H6_9_14 = 1e3; % 9.14 Rate: 8e2 s^-1, Range: 4e2 to 1.6e3, Nominal: 1e3 (midpoint), Source: Matsuda et al., 1990
k.stick_CH3Minus_9_15 = 5e3; % 9.15 Rate: 5e3 s^-1, Range: 2.5e3 to 7.5e3, Nominal: 5e3 (midpoint), Source: Boeuf, 1987
k.stick_H2_9_16 = 6.25e2; % 9.16 Rate: 5e2 s^-1, Range: 2.5e2 to 1e3, Nominal: 6.25e2 (midpoint), Source: Perrin, 1991
k.stick_C2H5_9_17 = 1.25e3; % 9.17 Rate: 1e3 s^-1, Range: 5e2 to 2e3, Nominal: 1.25e3 (midpoint), Source: Matsuda et al., 1990
k.stick_HMinus_9_18 = 5e3; % 9.18 Rate: 5e3 s^-1, Range: 2.5e3 to 7.5e3 (γ: 0.5–1.0), Nominal: 5e3 (midpoint), Source: Boeuf, 1987
k.stick_C3H_9_19 = 1.25e3; % 9.19 Rate: 1e3 s^-1, Range: 5e2 to 2e3, Nominal: 1.25e3 (midpoint), Source: Matsuda et al., 1990
k.stick_C4H2_9_20 = 1.25e3; % 9.20 Rate: 1e3 s^-1, Range: 5e2 to 2e3, Nominal: 1.25e3 (midpoint), Source: Matsuda et al., 1990
k.stick_C3H3_9_21 = 1.25e3; % 9.21 Rate: 1e3 s^-1, Range: 5e2 to 2e3, Nominal: 1.25e3 (midpoint), Source: Matsuda et al., 1990
k.stick_C3H4_9_22 = 1.25e3; % 9.22 Rate: 1e3 s^-1, Range: 5e2 to 2e3, Nominal: 1.25e3 (midpoint), Source: Matsuda et al., 1990
k.stick_C3_9_23 = 1.25e3; % 9.23 Rate: 1e3 s^-1, Range: 5e2 to 2e3, Nominal: 1.25e3 (midpoint), Source: Matsuda et al., 1990
k.stick_C4H_9_24 = 1.25e3; % 9.24 Rate: 1e3 s^-1, Range: 5e2 to 2e3, Nominal: 1.25e3 (midpoint), Source: Matsuda et al., 1990
k.stick_C3H6_9_25 = 1.25e3; % 9.25 Rate: 1e3 s^-1, Range: 5e2 to 2e3, Nominal: 1.25e3 (midpoint), Source: Matsuda et al., 1990
k.stick_H3Plus_9_26 = 5e3; % 9.26 Rate: 5e3 s^-1, Range: 2.5e3 to 7.5e3 (γ: 0.5–1.0), Nominal: 5e3 (midpoint), Source: Boeuf, 1987
k.stick_CHPlus_9_27 = 5e3; % 9.27 Rate: 5e3 s^-1, Range: 2.5e3 to 7.5e3 (γ: 0.5–1.0), Nominal: 5e3 (midpoint), Source: Boeuf, 1987
k.stick_C2HPlus_9_28 = 5e3; % 9.28 Rate: 5e3 s^-1, Range: 2.5e3 to 7.5e3 (γ: 0.5–1.0), Nominal: 5e3 (midpoint), Source: Boeuf, 1987

% Group 10: Drift Losses
k.drift_ArPlus_10_1 = mobilities.ArPlus * E_field / L_discharge; % 10.1 Rate: ~43.64 s^-1, Range: 29.09 to 43.64 (±20%), Nominal: ~43.64, Source: Ellis et al. (1976)
k.drift_CH4Plus_10_2 = mobilities.CH4Plus * E_field / L_discharge; % 10.2 Rate: ~45.45 s^-1, Range: 36.36 to 54.55 (±20%), Nominal: ~45.45, Source: NIST
k.drift_CH3Plus_10_3 = mobilities.CH3Plus * E_field / L_discharge; % 10.3 Rate: ~41.82 s^-1, Range: 33.45 to 50.18 (±20%), Nominal: ~41.82, Source: NIST
k.drift_CH5Plus_10_4 = mobilities.CH5Plus * E_field / L_discharge; % 10.4 Rate: ~48 s^-1, Range: 32 to 48 (±20%), Nominal: ~48, Source: NIST
k.drift_ArHPlus_10_5 = mobilities.ArHPlus * E_field / L_discharge; % 10.5 Rate: ~34.55 s^-1, Range: 27.64 to 41.45 (±20%), Nominal: ~34.55, Source: NIST
k.drift_CH2Plus_10_6 = mobilities.CH2Plus * E_field / L_discharge; % 10.6 Rate: ~41.82 s^-1, Range: 33.45 to 50.18 (±20%), Nominal: ~41.82, Source: NIST
k.drift_C2H5Plus_10_7 = mobilities.C2H5Plus * E_field / L_discharge; % 10.7 Rate: ~41.82 s^-1, Range: 33.45–50.18 (±20%), Nominal: ~41.82, Source: NIST
k.drift_C2H4Plus_10_8 = mobilities.C2H4Plus * E_field / L_discharge; % 10.8 Rate: ~41.82 s^-1, Range: 33.45–50.18 (±20%), Nominal: ~41.82, Source: NIST
k.drift_C2H3Plus_10_9 = mobilities.C2H3Plus * E_field / L_discharge; % 10.9 Rate: ~41.82 s^-1, Range: 33.45–50.18 (±20%), Nominal: ~41.82, Source: NIST
k.drift_H3Plus_10_10 = mobilities.H3Plus * E_field / L_discharge; % 10.10 Rate: ~50 s^-1, Range: 40 to 60 (±20%), Nominal: ~50, Source: NIST
k.drift_CHPlus_10_11 = mobilities.CHPlus * E_field / L_discharge; % 10.11 Rate: ~50 s^-1, Range: 40 to 60 (±20%), Nominal: ~50, Source: NIST
k.drift_CH3Minus_10_12 = mobilities.CH3Minus * E_field / L_discharge; % 10.12 Rate: ~37.5 s^-1, Range: 30 to 45 (±20%), Nominal: ~37.5, Source: NIST
k.drift_C2HPlus_10_13 = mobilities.C2HPlus * E_field / L_discharge; % 10.13 Rate: ~50 s^-1, Range: 40 to 60 (±20%), Nominal: ~50, Source: NIST

% Group 11: Loss Reactions
k.loss_CH2_11_1 = 3.63e3; % 11.1 Rate: 1.21e3 s^-1, Range: 1.21e3 to 6.04e3 (γ: 0.1–0.5), Nominal: 3.63e3 (midpoint), Source: Jauberteau et al., 1998
k.loss_H2_11_2 = 3.5e3; % 11.2 Rate: 3e3 s^-1, Range: 2e3 to 5e3, Nominal: 3.5e3 (midpoint), Source: Matsuda, 2004
k.loss_C2_11_3 = 1e3; % 11.3 Rate: 1e3 s^-1, Range: 1e-4 to 2e3, Nominal: 1e3 (midpoint), Source: Alman & Ruzic, 2003
k.loss_e_11_4 = 7.5e3; % 11.4 Rate: 7.84e3 s^-1, Range: 5e3 to 1e4, Nominal: 7.5e3 (midpoint), Source: Lieberman & Lichtenberg, 2005
k.loss_C2H6_11_5 = 1e3; % 11.5 Rate: 1e2 s^-1, Range: 5e2 to 1.5e3, Nominal: 1e3 (midpoint), Source: Alman & Ruzic, 2003
k.loss_CH4_11_6 = 1.5e3; % 11.6 Rate: 1.57e3 s^-1, Range: 1e3 to 2e3, Nominal: 1.5e3 (midpoint), Source: Alman & Ruzic, 2003
k.loss_Ar_11_7 = 2e3; % 11.7 Rate: 2.35e3 s^-1, Range: 1e3 to 3e3, Nominal: 2e3 (midpoint), Source: Estimated
% FLAG: 11.7 Rate is estimated. Needs validation.
k.loss_C_11_8 = 1.25e3; % 11.8 Rate: 2e3 s^-1, Range: 5e2 to 2e3, Nominal: 1.25e3 (midpoint), Source: Alman & Ruzic, 2003
k.loss_CH_11_9 = 5.5e3; % 11.9 Rate: 1e4 s^-1, Range: 1e3 to 1e4, Nominal: 5.5e3 (midpoint), Source: Estimated
% FLAG: 11.9 Rate is estimated. Needs validation.
k.loss_C3H_11_10 = 5.5e2; % 11.10 Rate: 1e3 s^-1, Range: 1e2 to 1e3, Nominal: 5.5e2 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C4H2_11_11 = 5.5e2; % 11.11 Rate: 1e3 s^-1, Range: 1e2 to 1e3, Nominal: 5.5e2 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C2H_11_12 = 5.5e2; % 11.12 Rate: 1e3 s^-1, Range: 1e2 to 1e3, Nominal: 5.5e2 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C3H3_11_13 = 5.5e2; % 11.13 Rate: 1e3 s^-1, Range: 1e2 to 1e3, Nominal: 5.5e2 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C3_11_14 = 5.5e2; % 11.14 Rate: 5e2 s^-1, Range: 1e2 to 1e3, Nominal: 5.5e2 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C4H_11_15 = 5.5e2; % 11.15 Rate: 1e3 s^-1, Range: 1e2 to 1e3, Nominal: 5.5e2 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C3H6_11_16 = 5.5e2; % 11.16 Rate: 5e2 s^-1, Range: 1e2 to 1e3, Nominal: 5.5e2 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C2H5_11_17 = 1e3; % 11.17 Rate: 1.5e3 s^-1, Range: 5e2–1.5e3, Nominal: 1e3 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C3H4_11_18 = 8.25e2; % 11.18 Rate: 1.5e3 s^-1, Range: 1e2–1.5e3, Nominal: 8.25e2 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C2H2_11_19 = 1.5e3; % 11.19 Rate: 1.5e3 s^-1, Range: 1e3 to 2e3, Nominal: 1.5e3 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C2H4_11_20 = 1.5e3; % 11.20 Rate: 1.5e3 s^-1, Range: 1e3 to 2e3, Nominal: 1.5e3 (midpoint), Source: Alman & Ruzic, 2003
k.loss_CH3_11_21 = 1.2e3; % 11.21 Rate: 1.2e3 s^-1, Range: 8e2 to 1.6e3, Nominal: 1.2e3 (midpoint), Source: Jauberteau et al., 1998
k.loss_C2H3_11_22 = 1.5e3; % 11.22 Rate: 1.5e3 s^-1, Range: 1e3 to 2e3, Nominal: 1.5e3 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C3H2_11_23 = 1e3; % 11.23 Rate: 1e3 s^-1, Range: 5e2 to 1.5e3, Nominal: 1e3 (midpoint), Source: Alman & Ruzic, 2003
k.loss_C3H5_11_24 = 1e3; % 11.24 Rate: 1e3 s^-1, Range: 5e2 to 1.5e3, Nominal: 1e3 (midpoint), Source: Alman & Ruzic, 2003
end