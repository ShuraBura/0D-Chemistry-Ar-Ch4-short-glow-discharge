function [R, tags] = build_reactions(params)
species = params.species;
ns = length(species);
k = params.k;
R = struct('reactants', {}, 'products', {}, 'rate', {}, 'tag', {});
tags = {};

function R = push(R, reactants, products, rate, tag)
    R(end+1) = struct('reactants', reactants, 'products', products, 'rate', rate, 'tag', tag);
    tags{end+1} = tag;
end

function vec = sto(species, varargin)
    vec = zeros(length(species), 1);
    for i = 1:2:length(varargin)
        idx = find(strcmp(species, varargin{i}));
        if isempty(idx)
            error('Species %s not found', varargin{i});
        end
        vec(idx) = varargin{i+1};
    end
end

% Group 1: Electron-Impact Reactions (Neutral Products)
R = push(R, sto(species, 'e', 1, 'CH4', 1), sto(species, 'e', 1, 'CH3', 1, 'H', 1), k.e_CH4_CH3_H_cm3_1_1, 'e_CH4_CH3_H_cm3_1_1'); % 1.1
R = push(R, sto(species, 'e', 1, 'CH4', 1), sto(species, 'e', 1, 'CH2', 1, 'H2', 1), k.e_CH4_CH2_H2_cm3_1_2, 'e_CH4_CH2_H2_cm3_1_2'); % 1.2
R = push(R, sto(species, 'e', 1, 'CH4', 1), sto(species, 'e', 1, 'CH', 1, 'H2', 1, 'H', 1), k.e_CH4_CH_H2_H_vib_cm3_1_3, 'e_CH4_CH_H2_H_vib_cm3_1_3'); % 1.3
R = push(R, sto(species, 'e', 1, 'H2', 1), sto(species, 'e', 1, 'H', 2), k.e_H2_H_H_cm3_1_4, 'e_H2_H_H_cm3_1_4'); % 1.4
R = push(R, sto(species, 'e', 1, 'CH3', 1), sto(species, 'e', 1, 'CH2', 1, 'H', 1), k.e_CH3_CH2_H_cm3_1_5, 'e_CH3_CH2_H_cm3_1_5'); % 1.5
R = push(R, sto(species, 'e', 1, 'C2H4', 1), sto(species, 'e', 1, 'C2H2', 1, 'H2', 1), k.e_C2H4_C2H2_H2_cm3_1_6, 'e_C2H4_C2H2_H2_cm3_1_6'); % 1.6
R = push(R, sto(species, 'e', 1, 'Ar', 1), sto(species, 'e', 1, 'ArStar', 1), k.e_Ar_ArStar_cm3_1_7, 'e_Ar_ArStar_cm3_1_7'); % 1.7
R = push(R, sto(species, 'e', 1, 'C2H6', 1), sto(species, 'e', 1, 'C2H4', 1, 'H2', 1), k.e_C2H6_C2H4_H2_cm3_1_8, 'e_C2H6_C2H4_H2_cm3_1_8'); % 1.8
R = push(R, sto(species, 'e', 1, 'C2H6', 1), sto(species, 'e', 1, 'C2H4', 1, 'H2', 1), k.e_C2H6_C2H4_H2_e_cm3_1_9, 'e_C2H6_C2H4_H2_e_cm3_1_9'); % 1.9
R = push(R, sto(species, 'e', 1, 'CH4', 1), sto(species, 'CH3Minus', 1, 'H', 1), k.e_CH4_CH3Minus_H_cm3_1_10, 'e_CH4_CH3Minus_H_cm3_1_10'); % 1.10
% FLAG: 1.10 Rate is low and estimated. Needs validation.
R = push(R, sto(species, 'e', 1, 'CH4', 1), sto(species, 'e', 1, 'CH', 1, 'H', 1, 'H2', 1), k.e_CH4_CH_H_H2_cm3_1_11, 'e_CH4_CH_H_H2_cm3_1_11'); % 1.11
R = push(R, sto(species, 'e', 1, 'CH', 1), sto(species, 'e', 1, 'C', 1, 'H', 1), k.e_CH_CH_C_H_e_cm3_1_12, 'e_CH_CH_C_H_e_cm3_1_12'); % 1.12
R = push(R, sto(species, 'e', 1, 'H2', 1), sto(species, 'HMinus', 1, 'H', 1), k.e_H2_HMinus_H_cm3_1_13, 'e_H2_HMinus_H_cm3_1_13'); % 1.13
R = push(R, sto(species, 'e', 1, 'CH3', 1), sto(species, 'CH3Minus', 1), k.e_CH3_CH3Minus_cm3_1_14, 'e_CH3_CH3Minus_cm3_1_14'); % 1.14
R = push(R, sto(species, 'e', 1, 'C2H4', 1), sto(species, 'C2H2', 1, 'H2', 1), k.e_C2H4_C2H2_H2_cm3_1_15, 'e_C2H4_C2H2_H2_cm3_1_15'); % 1.15
R = push(R, sto(species, 'e', 1, 'C2H2', 1), sto(species, 'C2', 1, 'H2', 1), k.e_C2H2_C2_H2_cm3_1_16, 'e_C2H2_C2_H2_cm3_1_16'); % 1.16
R = push(R, sto(species, 'e', 1, 'C2H4', 1), sto(species, 'C2H2', 1, 'H', 2), k.e_C2H4_C2H2_H_H_cm3_1_17, 'e_C2H4_C2H2_H_H_cm3_1_17'); % 1.17
R = push(R, sto(species, 'e', 1, 'C2H6', 1), sto(species, 'C2H2', 1, 'H2', 2), k.e_C2H6_C2H2_2H2_cm3_1_18, 'e_C2H6_C2H2_2H2_cm3_1_18'); % 1.18

% Group 2: Electron-Impact Ionization
R = push(R, sto(species, 'e', 1, 'CH4', 1), sto(species, 'CH3Plus', 1, 'H', 1, 'e', 2), k.e_CH4_CH3Plus_H_cm3_2_1, 'e_CH4_CH3Plus_H_cm3_2_1'); % 2.1
R = push(R, sto(species, 'e', 1, 'CH4', 1), sto(species, 'e', 2, 'CH4Plus', 1), k.e_CH4_CH4Plus_cm3_2_2, 'e_CH4_CH4Plus_cm3_2_2'); % 2.2
R = push(R, sto(species, 'e', 1, 'Ar', 1), sto(species, 'e', 2, 'ArPlus', 1), k.e_Ar_ArPlus_cm3_2_3, 'e_Ar_ArPlus_cm3_2_3'); % 2.3
R = push(R, sto(species, 'e', 1, 'ArStar', 1), sto(species, 'e', 2, 'ArPlus', 1), k.e_ArStar_ArPlus_cm3_2_4, 'e_ArStar_ArPlus_cm3_2_4'); % 2.4
R = push(R, sto(species, 'e', 1, 'C2H6', 1), sto(species, 'e', 2, 'C2H5Plus', 1, 'H', 1), k.e_C2H6_C2H5Plus_H_2e_cm3_2_5, 'e_C2H6_C2H5Plus_H_2e_cm3_2_5'); % 2.5
R = push(R, sto(species, 'e', 1, 'C2H4', 1), sto(species, 'e', 2, 'C2H4Plus', 1), k.e_C2H4_C2H4Plus_2e_cm3_2_6, 'e_C2H4_C2H4Plus_2e_cm3_2_6'); % 2.6
R = push(R, sto(species, 'e', 1, 'C2H4', 1), sto(species, 'e', 2, 'C2H3Plus', 1, 'H', 1), k.e_C2H4_C2H3Plus_H_2e_cm3_2_7, 'e_C2H4_C2H3Plus_H_2e_cm3_2_7'); % 2.7
R = push(R, sto(species, 'e', 1, 'C2H2', 1), sto(species, 'e', 2, 'C2HPlus', 1, 'H', 1), k.e_C2H2_C2HPlus_2e_cm3_2_8, 'e_C2H2_C2HPlus_2e_cm3_2_8'); % 2.8

% Group 3: ArStar Reactions
R = push(R, sto(species, 'ArStar', 1, 'CH4', 1), sto(species, 'Ar', 1, 'CH3', 1, 'H', 1), k.ArStar_CH4_CH3_H_cm3_3_1, 'ArStar_CH4_CH3_H_cm3_3_1'); % 3.1
R = push(R, sto(species, 'ArStar', 1, 'H2', 1), sto(species, 'Ar', 1, 'H', 2), k.ArStar_H2_H_H_cm3_3_2, 'ArStar_H2_H_H_cm3_3_2'); % 3.2
R = push(R, sto(species, 'ArStar', 1, 'H2', 1), sto(species, 'ArHPlus', 1, 'H', 1, 'e', 1), k.ArStar_H2_ArHPlus_H_cm3_3_3, 'ArStar_H2_ArHPlus_H_cm3_3_3'); % 3.3
R = push(R, sto(species, 'ArStar', 1, 'CH2', 1), sto(species, 'Ar', 1, 'CH', 1, 'H', 1), k.ArStar_CH2_CH_H_cm3_3_4, 'ArStar_CH2_CH_H_cm3_3_4'); % 3.4
R = push(R, sto(species, 'ArStar', 1, 'C', 1), sto(species, 'Ar', 1, 'C', 1), k.ArStar_C_Ar_CStar_cm3_3_5, 'ArStar_C_Ar_CStar_cm3_3_5'); % 3.5
R = push(R, sto(species, 'ArStar', 1, 'H', 1), sto(species, 'ArHPlus', 1, 'e', 1), k.ArStar_H_ArHPlus_e_cm3_3_6, 'ArStar_H_ArHPlus_e_cm3_3_6'); % 3.6
R = push(R, sto(species, 'ArStar', 1, 'CH3', 1), sto(species, 'Ar', 1, 'CH2', 1, 'H', 1), k.ArStar_CH3_Ar_CH2_H_cm3_3_7, 'ArStar_CH3_Ar_CH2_H_cm3_3_7'); % 3.7
R = push(R, sto(species, 'ArStar', 1, 'C2', 1), sto(species, 'Ar', 1, 'C2', 1), k.ArStar_C2_Ar_C2_cm3_3_8, 'ArStar_C2_Ar_C2_cm3_3_8'); % 3.8
R = push(R, sto(species, 'ArStar', 1, 'C2H4', 1), sto(species, 'Ar', 1, 'C2H4', 1), k.ArStar_C2H4_Ar_C2H4_cm3_3_9, 'ArStar_C2H4_Ar_C2H4_cm3_3_9'); % 3.9
R = push(R, sto(species, 'ArStar', 1, 'CH2', 1), sto(species, 'Ar', 1, 'CH2', 1), k.ArStar_CH2_Ar_CH2_cm3_3_10, 'ArStar_CH2_Ar_CH2_cm3_3_10'); % 3.10
R = push(R, sto(species, 'ArStar', 1, 'CH3', 1), sto(species, 'Ar', 1, 'CH3', 1), k.ArStar_CH3_Ar_CH3_cm3_3_11, 'ArStar_CH3_Ar_CH3_cm3_3_11'); % 3.11
R = push(R, sto(species, 'ArStar', 1), sto(species, 'Ar', 1), k.ArStar_M_Ar_3_12, 'ArStar_M_Ar_3_12'); % 3.12
% FLAG: 3.12 Rate is estimated. Needs validation.
R = push(R, sto(species, 'ArStar', 1, 'Ar', 1), sto(species, 'Ar', 2), k.ArStar_Ar_Ar_Ar_cm3_3_13, 'ArStar_Ar_Ar_Ar_cm3_3_13'); % 3.13
R = push(R, sto(species, 'CH', 1, 'ArStar', 1), sto(species, 'C', 1, 'H', 1, 'Ar', 1), k.CH_ArStar_C_H_Ar_cm3_3_14, 'CH_ArStar_C_H_Ar_cm3_3_14'); % 3.14
R = push(R, sto(species, 'CH', 1, 'ArStar', 1), sto(species, 'CH', 1, 'Ar', 1), k.CH_ArStar_CH_Ar_cm3_3_15, 'CH_ArStar_CH_Ar_cm3_3_15'); % 3.15
R = push(R, sto(species, 'ArStar', 1, 'CH3', 1), sto(species, 'Ar', 1, 'CH2', 1, 'H', 1), k.ArStar_CH3_CH2_H_Ar_cm3_3_16, 'ArStar_CH3_CH2_H_Ar_cm3_3_16'); % 3.16
R = push(R, sto(species, 'ArStar', 1, 'e', 1), sto(species, 'Ar', 1, 'e', 1), k.ArStar_e_Ar_e_cm3_3_17, 'ArStar_e_Ar_e_cm3_3_17'); % 3.17
R = push(R, sto(species, 'ArStar', 1, 'H2', 1), sto(species, 'Ar', 1, 'H2', 1), k.ArStar_H2_Ar_H2Star_cm3_3_18, 'ArStar_H2_Ar_H2Star_cm3_3_18'); % 3.18
R = push(R, sto(species, 'ArStar', 1, 'CH4', 1), sto(species, 'Ar', 1, 'CH4', 1), k.ArStar_CH4_Ar_CH4Star_cm3_3_19, 'ArStar_CH4_Ar_CH4Star_cm3_3_19'); % 3.19
R = push(R, sto(species, 'ArStar', 1, 'CH3', 1), sto(species, 'Ar', 1, 'CH3', 1), k.ArStar_CH3_Ar_CH3Star_cm3_3_20, 'ArStar_CH3_Ar_CH3Star_cm3_3_20'); % 3.20
R = push(R, sto(species, 'ArStar', 1, 'H2', 1), sto(species, 'Ar', 1, 'H', 2), k.ArStar_H2_Ar_H_H_cm3_3_21, 'ArStar_H2_Ar_H_H_cm3_3_21'); % 3.21
R = push(R, sto(species, 'ArStar', 1, 'C2H2', 1), sto(species, 'Ar', 1, 'C2H2Star', 1), k.ArStar_C2H2_Ar_C2H2Star_cm3_3_22, 'ArStar_C2H2_Ar_C2H2Star_cm3_3_22'); % 3.22
R = push(R, sto(species, 'ArStar', 1, 'CH3', 1), sto(species, 'Ar', 1, 'CH2', 1, 'H', 1), k.ArStar_CH3_Ar_CH2_H_cm3_3_23, 'ArStar_CH3_Ar_CH2_H_cm3_3_23'); % 3.23
R = push(R, sto(species, 'ArStar', 1, 'C2H4', 1), sto(species, 'Ar', 1, 'C2H4', 1), k.ArStar_C2H4_Ar_C2H4Star_cm3_3_24, 'ArStar_C2H4_Ar_C2H4Star_cm3_3_24'); % 3.24

% Group 4: Penning Ionization
R = push(R, sto(species, 'ArStar', 1, 'CH4', 1), sto(species, 'Ar', 1, 'CH4Plus', 1, 'e', 1), k.ArStar_CH4_CH4Plus_cm3_4_1, 'ArStar_CH4_CH4Plus_cm3_4_1'); % 4.1
R = push(R, sto(species, 'ArStar', 1, 'CH4', 1), sto(species, 'Ar', 1, 'CH3Plus', 1, 'H', 1, 'e', 1), k.ArStar_CH4_CH3Plus_H_cm3_4_2, 'ArStar_CH4_CH3Plus_H_cm3_4_2'); % 4.2
R = push(R, sto(species, 'ArStar', 1, 'Ar', 1), sto(species, 'ArPlus', 1, 'Ar', 1, 'e', 1), k.ArStar_Ar_ArPlus_cm3_4_3, 'ArStar_Ar_ArPlus_cm3_4_3'); % 4.3
R = push(R, sto(species, 'ArStar', 1, 'CH4', 1), sto(species, 'ArPlus', 1, 'CH3', 1, 'H', 1, 'e', 1), k.ArStar_CH4_ArPlus_CH3_H_e_cm3_4_4, 'ArStar_CH4_ArPlus_CH3_H_e_cm3_4_4'); % 4.4
R = push(R, sto(species, 'ArStar', 1, 'CH3', 1), sto(species, 'ArPlus', 1, 'CH2', 1, 'H', 1, 'e', 1), k.ArStar_CH3_ArPlus_CH2_H_e_cm3_4_5, 'ArStar_CH3_ArPlus_CH2_H_e_cm3_4_5'); % 4.5
R = push(R, sto(species, 'ArStar', 1, 'H', 1), sto(species, 'ArPlus', 1, 'HMinus', 1, 'e', 1), k.ArStar_H_ArPlus_HMinus_cm3_4_6, 'ArStar_H_ArPlus_HMinus_cm3_4_6'); % 4.6
R = push(R, sto(species, 'ArStar', 1, 'CH2', 1), sto(species, 'ArPlus', 1, 'CH', 1, 'H', 1, 'e', 1), k.ArStar_CH2_ArPlus_CH_H_e_cm3_4_7, 'ArStar_CH2_ArPlus_CH_H_e_cm3_4_7'); % 4.7
R = push(R, sto(species, 'ArStar', 1, 'H2', 1), sto(species, 'ArPlus', 1, 'H2', 1, 'e', 1), k.ArStar_H2_ArPlus_H2_e_cm3_4_8, 'ArStar_H2_ArPlus_H2_e_cm3_4_8'); % 4.8
R = push(R, sto(species, 'ArStar', 1, 'C2H2', 1), sto(species, 'ArPlus', 1, 'C2H2', 1, 'e', 1), k.ArStar_C2H2_ArPlus_C2H2_e_cm3_4_9, 'ArStar_C2H2_ArPlus_C2H2_e_cm3_4_9'); % 4.9
R = push(R, sto(species, 'ArStar', 1, 'C2H5', 1), sto(species, 'ArPlus', 1, 'C2H5', 1, 'e', 1), k.ArStar_C2H5_ArPlus_C2H5_e_cm3_4_10, 'ArStar_C2H5_ArPlus_C2H5_e_cm3_4_10'); % 4.10
R = push(R, sto(species, 'ArStar', 1, 'H', 1), sto(species, 'ArPlus', 1, 'H', 1, 'e', 1), k.ArStar_H_ArPlus_H_e_cm3_4_11, 'ArStar_H_ArPlus_H_e_cm3_4_11'); % 4.11
R = push(R, sto(species, 'ArStar', 1, 'CH4', 1), sto(species, 'ArPlus', 1, 'CH4', 1, 'e', 1), k.ArStar_CH4_ArPlus_CH4_e_cm3_4_12, 'ArStar_CH4_ArPlus_CH4_e_cm3_4_12'); % 4.12
R = push(R, sto(species, 'ArStar', 1, 'CH3', 1), sto(species, 'ArPlus', 1, 'CH3', 1, 'e', 1), k.ArStar_CH3_ArPlus_CH3_e_cm3_4_13, 'ArStar_CH3_ArPlus_CH3_e_cm3_4_13'); % 4.13
R = push(R, sto(species, 'ArStar', 1, 'C2H4', 1), sto(species, 'ArPlus', 1, 'C2H4', 1, 'e', 1), k.ArStar_C2H4_ArPlus_C2H4_e_cm3_4_14, 'ArStar_C2H4_ArPlus_C2H4_e_cm3_4_14'); % 4.14
R = push(R, sto(species, 'ArStar', 1, 'C2H5', 1), sto(species, 'ArPlus', 1, 'C2H5', 1, 'e', 1), k.ArStar_C2H5_ArPlus_C2H5_e_cm3_4_15, 'ArStar_C2H5_ArPlus_C2H5_e_cm3_4_15'); % 4.15
R = push(R, sto(species, 'ArStar', 1, 'CH2', 1), sto(species, 'ArPlus', 1, 'CH2', 1, 'e', 1), k.ArStar_CH2_ArPlus_CH2_e_cm3_4_16, 'ArStar_CH2_ArPlus_CH2_e_cm3_4_16'); % 4.16
R = push(R, sto(species, 'ArStar', 1, 'C2H6', 1), sto(species, 'ArPlus', 1, 'C2H6', 1, 'e', 1), k.ArStar_C2H6_ArPlus_C2H6_e_cm3_4_17, 'ArStar_C2H6_ArPlus_C2H6_e_cm3_4_17'); % 4.17

% Group 5: Ion-Neutral Reactions
R = push(R, sto(species, 'ArPlus', 1, 'CH4', 1), sto(species, 'CH3Plus', 1, 'H', 1, 'Ar', 1), k.ArPlus_CH4_CH3Plus_H_cm3_5_1, 'ArPlus_CH4_CH3Plus_H_cm3_5_1'); % 5.1
R = push(R, sto(species, 'CH3Plus', 1, 'CH4', 1), sto(species, 'CH5Plus', 1, 'CH2', 1), k.CH3Plus_CH4_CH5Plus_CH2_cm3_5_2, 'CH3Plus_CH4_CH5Plus_CH2_cm3_5_2'); % 5.2
R = push(R, sto(species, 'ArPlus', 1, 'CH3', 1), sto(species, 'CH3Plus', 1, 'Ar', 1), k.ArPlus_CH3_CH3Plus_cm3_5_3, 'ArPlus_CH3_CH3Plus_cm3_5_3'); % 5.3
R = push(R, sto(species, 'CH', 2), sto(species, 'C2', 1, 'H2', 1), k.CH_CH_C2_H2_cm3_5_4, 'CH_CH_C2_H2_cm3_5_4'); % 5.4
R = push(R, sto(species, 'ArPlus', 1, 'CH4', 1), sto(species, 'Ar', 1, 'CH4Plus', 1), k.ArPlus_CH4_Ar_CH4Plus_cm3_5_5, 'ArPlus_CH4_Ar_CH4Plus_cm3_5_5'); % 5.5
R = push(R, sto(species, 'CH4Plus', 1, 'H2', 1), sto(species, 'CH5Plus', 1, 'H', 1), k.CH4Plus_H2_CH5Plus_H_cm3_5_6, 'CH4Plus_H2_CH5Plus_H_cm3_5_6'); % 5.6
R = push(R, sto(species, 'ArPlus', 1, 'H2', 1), sto(species, 'ArHPlus', 1, 'H', 1), k.ArPlus_H2_ArHPlus_H_cm3_5_7, 'ArPlus_H2_ArHPlus_H_cm3_5_7'); % 5.7
R = push(R, sto(species, 'ArHPlus', 1, 'CH4', 1), sto(species, 'Ar', 1, 'CH5Plus', 1), k.ArHPlus_CH4_Ar_CH5Plus_cm3_5_8, 'ArHPlus_CH4_Ar_CH5Plus_cm3_5_8'); % 5.8
R = push(R, sto(species, 'ArPlus', 1, 'CH4', 1), sto(species, 'ArHPlus', 1, 'CH3', 1), k.ArPlus_CH4_ArHPlus_CH3_cm3_5_9, 'ArPlus_CH4_ArHPlus_CH3_cm3_5_9'); % 5.9
R = push(R, sto(species, 'CH2', 1, 'CH3Plus', 1), sto(species, 'CH3', 1, 'CH2Plus', 1), k.CH2_CH3Plus_CH3_CH2Plus_cm3_5_10, 'CH2_CH3Plus_CH3_CH2Plus_cm3_5_10'); % 5.10
R = push(R, sto(species, 'CH3Plus', 1, 'CH4', 1), sto(species, 'C2H5Plus', 1, 'H2', 1), k.CH3Plus_CH4_C2H5Plus_H2_cm3_5_11, 'CH3Plus_CH4_C2H5Plus_H2_cm3_5_11'); % 5.11
R = push(R, sto(species, 'CH5Plus', 1, 'C2H4', 1), sto(species, 'C2H5Plus', 1, 'CH4', 1), k.CH5Plus_C2H4_C2H5Plus_CH4_cm3_5_12, 'CH5Plus_C2H4_C2H5Plus_CH4_cm3_5_12'); % 5.12

% Group 6: Dissociative Recombination
R = push(R, sto(species, 'ArPlus', 1, 'e', 1), sto(species, 'Ar', 1), k.ArPlus_e_Ar_cm3_6_1, 'ArPlus_e_Ar_cm3_6_1'); % 6.1
R = push(R, sto(species, 'CH3Plus', 1, 'e', 1), sto(species, 'CH3', 1), k.CH3Plus_e_CH3_cm3_6_2, 'CH3Plus_e_CH3_cm3_6_2'); % 6.2
R = push(R, sto(species, 'CH5Plus', 1, 'e', 1), sto(species, 'CH4', 1, 'H', 1), k.CH5Plus_e_CH4_H_cm3_6_3, 'CH5Plus_e_CH4_H_cm3_6_3'); % 6.3
R = push(R, sto(species, 'e', 1, 'CH4Plus', 1), sto(species, 'CH3', 1, 'H', 1), k.e_CH4Plus_CH3_H_cm3_6_4, 'e_CH4Plus_CH3_H_cm3_6_4'); % 6.4
R = push(R, sto(species, 'CH3Minus', 1, 'ArPlus', 1), sto(species, 'CH3', 1, 'Ar', 1), k.CH3Minus_ArPlus_CH3_Ar_cm3_6_5, 'CH3Minus_ArPlus_CH3_Ar_cm3_6_5'); % 6.5
% FLAG: 6.5 Rate is estimated. Needs validation.
R = push(R, sto(species, 'CH3Minus', 1, 'CH4Plus', 1), sto(species, 'CH4', 1, 'CH3', 1), k.CH3Minus_CH4Plus_CH4_CH3_cm3_6_6, 'CH3Minus_CH4Plus_CH4_CH3_cm3_6_6'); % 6.6
% FLAG: 6.6 Rate is estimated. Needs validation.
R = push(R, sto(species, 'CH3Minus', 1, 'CH3Plus', 1), sto(species, 'CH4', 1, 'CH2', 1), k.CH3Minus_CH3Plus_CH4_CH2_cm3_6_7, 'CH3Minus_CH3Plus_CH4_CH2_cm3_6_7'); % 6.7
% FLAG: 6.7 Rate is estimated. Needs validation.
R = push(R, sto(species, 'CH5Plus', 1, 'e', 1), sto(species, 'CH3', 1, 'H2', 1), k.CH5Plus_e_CH3_H2_cm3_6_8, 'CH5Plus_e_CH3_H2_cm3_6_8'); % 6.8
R = push(R, sto(species, 'CH4Plus', 1, 'e', 1), sto(species, 'CH2', 1, 'H2', 1), k.e_CH4Plus_CH2_H2_cm3_6_9, 'e_CH4Plus_CH2_H2_cm3_6_9'); % 6.9
R = push(R, sto(species, 'CH5Plus', 1, 'e', 1), sto(species, 'CH2', 1, 'H2', 1, 'H', 1), k.CH5Plus_e_CH2_H2_H_cm3_6_10, 'CH5Plus_e_CH2_H2_H_cm3_6_10'); % 6.10
R = push(R, sto(species, 'CH4Plus', 1, 'e', 1), sto(species, 'CH', 1, 'H2', 1, 'H', 1), k.e_CH4Plus_CH_H2_H_cm3_6_11, 'e_CH4Plus_CH_H2_H_cm3_6_11'); % 6.11
R = push(R, sto(species, 'CH5Plus', 1, 'e', 1), sto(species, 'CH3', 1, 'H', 2), k.CH5Plus_e_CH3_2H_cm3_6_12, 'CH5Plus_e_CH3_2H_cm3_6_12'); % 6.12
R = push(R, sto(species, 'CH4Plus', 1, 'e', 1), sto(species, 'C', 1, 'H2', 2), k.e_CH4Plus_C_2H2_cm3_6_13, 'e_CH4Plus_C_2H2_cm3_6_13'); % 6.13
R = push(R, sto(species, 'C2H5Plus', 1, 'e', 1), sto(species, 'C2H4', 1, 'H', 1), k.C2H5Plus_e_C2H4_H_cm3_6_14, 'C2H5Plus_e_C2H4_H_cm3_6_14'); % 6.14
R = push(R, sto(species, 'C2H4Plus', 1, 'e', 1), sto(species, 'C2H2', 1, 'H2', 1), k.C2H4Plus_e_C2H2_H2_cm3_6_15, 'C2H4Plus_e_C2H2_H2_cm3_6_15'); % 6.15
R = push(R, sto(species, 'C2H3Plus', 1, 'e', 1), sto(species, 'C2H2', 1, 'H', 1), k.C2H3Plus_e_C2H2_H_cm3_6_16, 'C2H3Plus_e_C2H2_H_cm3_6_16'); % 6.16
R = push(R, sto(species, 'HMinus', 1, 'ArPlus', 1), sto(species, 'H', 1, 'Ar', 1), k.HMinus_ArPlus_H_Ar_cm3_6_17, 'HMinus_ArPlus_H_Ar_cm3_6_17'); % 6.17
R = push(R, sto(species, 'C2HPlus', 1, 'e', 1), sto(species, 'C2', 1, 'H', 1), k.C2HPlus_e_C2_H_cm3_6_18, 'C2HPlus_e_C2_H_cm3_6_18'); % 6.18
R = push(R, sto(species, 'HMinus', 1, 'CH5Plus', 1), sto(species, 'CH4', 1, 'H2', 1, 'H', 1), k.HMinus_CH5Plus_CH4_H2_H_cm3_6_19, 'HMinus_CH5Plus_CH4_H2_H_cm3_6_19'); % 6.19
R = push(R, sto(species, 'CH4Plus', 1, 'HMinus', 1), sto(species, 'CH4', 1, 'H', 1), k.CH4Plus_HMinus_CH4_H_cm3_6_20, 'CH4Plus_HMinus_CH4_H_cm3_6_20'); % 6.20
R = push(R, sto(species, 'CH3Plus', 1, 'HMinus', 1), sto(species, 'CH4', 1, 'H2', 1), k.CH3Plus_HMinus_CH4_H2_cm3_6_21, 'CH3Plus_HMinus_CH4_H2_cm3_6_21'); % 6.21
R = push(R, sto(species, 'C2H5Plus', 1, 'HMinus', 1), sto(species, 'C2H6', 1, 'H', 1), k.C2H5Plus_HMinus_C2H6_H_cm3_6_22, 'C2H5Plus_HMinus_C2H6_H_cm3_6_22'); % 6.22
R = push(R, sto(species, 'ArHPlus', 1, 'HMinus', 1), sto(species, 'Ar', 1, 'H2', 1, 'H', 1), k.ArHPlus_HMinus_Ar_H2_H_cm3_6_23, 'ArHPlus_HMinus_Ar_H2_H_cm3_6_23'); % 6.23
R = push(R, sto(species, 'CH5Plus', 1, 'CH3Minus', 1), sto(species, 'CH4', 2, 'H', 1), k.CH5Plus_CH3Minus_CH4_CH4_H_cm3_6_24, 'CH5Plus_CH3Minus_CH4_CH4_H_cm3_6_24'); % 6.24
R = push(R, sto(species, 'CH4Plus', 1, 'CH3Minus', 1), sto(species, 'CH4', 1, 'CH3', 1, 'H', 1), k.CH4Plus_CH3Minus_CH4_CH3_H_cm3_6_25, 'CH4Plus_CH3Minus_CH4_CH3_H_cm3_6_25'); % 6.25
R = push(R, sto(species, 'CH3Plus', 1, 'CH3Minus', 1), sto(species, 'CH4', 1, 'CH2', 1, 'H', 1), k.CH3Plus_CH3Minus_CH4_CH2_H_cm3_6_26, 'CH3Plus_CH3Minus_CH4_CH2_H_cm3_6_26'); % 6.26
R = push(R, sto(species, 'C2H5Plus', 1, 'CH3Minus', 1), sto(species, 'C2H6', 1, 'H', 1), k.C2H5Plus_CH3Minus_C2H6_H_cm3_6_27, 'C2H5Plus_CH3Minus_C2H6_H_cm3_6_27'); % 6.27
R = push(R, sto(species, 'C2H5Plus', 1, 'e', 1), sto(species, 'C2H4', 1, 'H', 1), k.C2H5Plus_e_C2H4_H_cm3_6_28, 'C2H5Plus_e_C2H4_H_cm3_6_28'); % 6.28

% Group 7: Neutral-Neutral Reactions
R = push(R, sto(species, 'CH2', 1, 'H', 1), sto(species, 'CH', 1, 'H2', 1), k.CH2_H_CH_H2_cm3_7_1, 'CH2_H_CH_H2_cm3_7_1'); % 7.1
R = push(R, sto(species, 'CH2', 1, 'H', 1), sto(species, 'C', 1, 'H2', 1, 'H', 1), k.CH2_H_C_H2_H_cm3_7_2, 'CH2_H_C_H2_H_cm3_7_2'); % 7.2
% FLAG: 7.2 Rate is estimated. Needs discussion.
R = push(R, sto(species, 'CH', 1, 'H', 1), sto(species, 'C', 1, 'H2', 1), k.CH_H_C_H2_cm3_7_3, 'CH_H_C_H2_cm3_7_3'); % 7.3
R = push(R, sto(species, 'C', 1, 'CH', 1), sto(species, 'C2', 1, 'H', 1), k.C_CH_C2_H_cm3_7_4, 'C_CH_C2_H_cm3_7_4'); % 7.4
R = push(R, sto(species, 'CH', 1, 'CH3', 1), sto(species, 'C2H4', 1), k.CH_CH3_C2H4_cm3_7_5, 'CH_CH3_C2H4_cm3_7_5'); % 7.5
R = push(R, sto(species, 'C2', 1, 'H', 1), sto(species, 'CH', 1, 'C', 1), k.C2_H_CH_C_cm3_7_6, 'C2_H_CH_C_cm3_7_6'); % 7.6
R = push(R, sto(species, 'CH', 1, 'CH2', 1), sto(species, 'C2H2', 1, 'H', 1), k.CH_CH2_C2H2_H_cm3_7_7, 'CH_CH2_C2H2_H_cm3_7_7'); % 7.7
R = push(R, sto(species, 'C', 1, 'CH3', 1), sto(species, 'C2H2', 1, 'H', 1), k.C_CH3_C2_H2_H_cm3_7_8, 'C_CH3_C2_H2_H_cm3_7_8'); % 7.8
R = push(R, sto(species, 'CH', 1, 'C', 1), sto(species, 'C2', 1, 'H', 1), k.CH_C_C2_H_cm3_7_9, 'CH_C_C2_H_cm3_7_9'); % 7.9
R = push(R, sto(species, 'CH', 1, 'CH3', 1), sto(species, 'C2H3', 1, 'H', 1), k.CH_CH3_C2H3_H_cm3_7_10, 'CH_CH3_C2H3_H_cm3_7_10'); % 7.10
R = push(R, sto(species, 'CH', 1, 'Ar', 1), sto(species, 'Ar', 1, 'C', 1, 'H', 1), k.CH_Ar_Ar_C_H_cm3_7_11, 'CH_Ar_Ar_C_H_cm3_7_11'); % 7.11
% FLAG: 7.11 Rate from UMIST, not Baulch. Needs verification.
R = push(R, sto(species, 'C', 1, 'H', 1), sto(species, 'CH', 1), k.C_H_CH_cm3_7_12, 'C_H_CH_cm3_7_12'); % 7.12
R = push(R, sto(species, 'CH2', 2), sto(species, 'C2H4', 1), k.CH2_CH2_C2H4_cm3_7_13, 'CH2_CH2_C2H4_cm3_7_13'); % 7.13
R = push(R, sto(species, 'CH3', 1, 'CH2', 1), sto(species, 'C2H5', 1), k.CH3_CH2_C2H5_cm3_7_14, 'CH3_CH2_C2H5_cm3_7_14'); % 7.14
R = push(R, sto(species, 'CH2', 2), sto(species, 'C2H2', 1, 'H2', 1), k.CH2_CH2_C2H2_H2_cm3_7_15, 'CH2_CH2_C2H2_H2_cm3_7_15'); % 7.15
R = push(R, sto(species, 'CH3', 1, 'CH', 1), sto(species, 'C2H2', 1, 'H2', 1), k.CH3_CH_C2H2_H2_cm3_7_16, 'CH3_CH_C2H2_H2_cm3_7_16'); % 7.16
R = push(R, sto(species, 'CH2', 1, 'C', 1), sto(species, 'C2H2', 1), k.CH2_C_C2H2_cm3_7_17, 'CH2_C_C2H2_cm3_7_17'); % 7.17
R = push(R, sto(species, 'CH', 1, 'C2H4', 1), sto(species, 'C2H2', 1, 'CH3', 1), k.CH_C2H4_C2H2_CH3_cm3_7_18, 'CH_C2H4_C2H2_CH3_cm3_7_18'); % 7.18
R = push(R, sto(species, 'C2H2', 1, 'C', 1), sto(species, 'C2', 1, 'CH2', 1), k.C2H2_C_C2_CH2_cm3_7_19, 'C2H2_C_C2_CH2_cm3_7_19'); % 7.19
R = push(R, sto(species, 'CH', 1, 'CH4', 1), sto(species, 'C2H4', 1, 'H', 1), k.CH_CH4_C2H4_H_cm3_7_20, 'CH_CH4_C2H4_H_cm3_7_20'); % 7.20
R = push(R, sto(species, 'CH', 1, 'H', 1), sto(species, 'CH2', 1), k.CH_H_CH2_cm3_7_21, 'CH_H_CH2_cm3_7_21'); % 7.21
R = push(R, sto(species, 'CH', 1, 'C2H2', 1), sto(species, 'C3H2', 1, 'H', 1), k.CH_C2H2_C3H2_H_cm3_7_22, 'CH_C2H2_C3H2_H_cm3_7_22'); % 7.22
R = push(R, sto(species, 'CH', 1, 'CH3', 1), sto(species, 'C2H2', 1, 'H2', 1), k.CH_CH3_C2H2_H2_cm3_7_23, 'CH_CH3_C2H2_H2_cm3_7_23'); % 7.23
R = push(R, sto(species, 'CH', 1, 'C', 1), sto(species, 'C2', 1, 'H2', 1), k.CH_C_C2_H2_cm3_7_24, 'CH_C_C2_H2_cm3_7_24'); % 7.24
R = push(R, sto(species, 'H', 1, 'CH4', 1), sto(species, 'CH3', 1, 'H2', 1), k.H_CH4_CH3_H2_cm3_7_25, 'H_CH4_CH3_H2_cm3_7_25'); % 7.25
R = push(R, sto(species, 'CH2', 1, 'CH', 1), sto(species, 'C2H2', 1, 'H', 1), k.CH2_CH_C2_H2_H_cm3_7_26, 'CH2_CH_C2_H2_H_cm3_7_26'); % 7.26
R = push(R, sto(species, 'CH', 1, 'C2H2', 1), sto(species, 'C3H', 1, 'H2', 1), k.CH_C2H2_C3H_H2_cm3_7_27, 'CH_C2H2_C3H_H2_cm3_7_27'); % 7.27
R = push(R, sto(species, 'CH', 1, 'C3H', 1), sto(species, 'C4H2', 1, 'H', 1), k.CH_C3H_C4H2_H_cm3_7_28, 'CH_C3H_C4H2_H_cm3_7_28'); % 7.28
R = push(R, sto(species, 'CH', 1, 'C2H2', 1), sto(species, 'C2H', 1, 'CH2', 1), k.CH_C2H2_C2H_CH2_cm3_7_29, 'CH_C2H2_C2H_CH2_cm3_7_29'); % 7.29
R = push(R, sto(species, 'CH', 1, 'H2', 1), sto(species, 'CH2', 1, 'H', 1), k.CH_H2_CH2_H_cm3_7_30, 'CH_H2_CH2_H_cm3_7_30'); % 7.30
R = push(R, sto(species, 'CH', 1, 'C2H3', 1), sto(species, 'C3H3', 1, 'H', 1), k.CH_C2H3_C3H3_H_cm3_7_31, 'CH_C2H3_C3H3_H_cm3_7_31'); % 7.31
R = push(R, sto(species, 'CH', 1, 'C2H4', 1), sto(species, 'C3H4', 1, 'H', 1), k.CH_C2H4_C3H4_H_cm3_7_32, 'CH_C2H4_C3H4_H_cm3_7_32'); % 7.32
R = push(R, sto(species, 'CH', 1, 'C2', 1), sto(species, 'C3', 1, 'H', 1), k.CH_C2_C3_H_cm3_7_33, 'CH_C2_C3_H_cm3_7_33'); % 7.33
R = push(R, sto(species, 'CH', 1, 'C2H5', 1), sto(species, 'C3H5', 1, 'H', 1), k.CH_C2H5_C3H5_H_cm3_7_34, 'CH_C2H5_C3H5_H_cm3_7_34'); % 7.34
R = push(R, sto(species, 'CH', 1, 'C3H2', 1), sto(species, 'C4H', 1, 'H2', 1), k.CH_C3H2_C4H_H2_cm3_7_35, 'CH_C3H2_C4H_H2_cm3_7_35'); % 7.35
R = push(R, sto(species, 'CH3', 1, 'H', 1), sto(species, 'CH2', 1, 'H2', 1), k.CH3_H_CH2_H2_cm3_7_36, 'CH3_H_CH2_H2_cm3_7_36'); % 7.36
R = push(R, sto(species, 'CH', 1, 'C2H6', 1), sto(species, 'C3H6', 1, 'H', 1), k.CH_C2H6_C3H6_H_cm3_7_37, 'CH_C2H6_C3H6_H_cm3_7_37'); % 7.37
R = push(R, sto(species, 'CH3', 2), sto(species, 'CH2', 1, 'CH4', 1), k.CH3_CH3_CH2_CH4_cm3_7_38, 'CH3_CH3_CH2_CH4_cm3_7_38'); % 7.38
R = push(R, sto(species, 'CH', 1, 'CH4', 1), sto(species, 'CH2', 1, 'CH3', 1), k.CH_CH4_CH2_CH3_cm3_7_39, 'CH_CH4_CH2_CH3_cm3_7_39'); % 7.39
R = push(R, sto(species, 'CH3', 2), sto(species, 'C2H6', 1), k.CH3_CH3_C2H6_cm3_7_40, 'CH3_CH3_C2H6_cm3_7_40'); % 7.40
R = push(R, sto(species, 'CH', 1, 'C2H5', 1), sto(species, 'C3H6', 1), k.CH_C2H5_C3H6_cm3_7_41, 'CH_C2H5_C3H6_cm3_7_41'); % 7.41
R = push(R, sto(species, 'CH2', 2), sto(species, 'C2H2', 1, 'H2', 1), k.CH2_CH2_C2H2_H2_cm3_7_42, 'CH2_CH2_C2H2_H2_cm3_7_42'); % 7.42
R = push(R, sto(species, 'CH', 1, 'C', 1), sto(species, 'C2', 1, 'H', 1), k.CH_C_C2_H_cm3_7_43, 'CH_C_C2_H_cm3_7_43'); % 7.43
R = push(R, sto(species, 'CH', 2), sto(species, 'C2', 1, 'H2', 1), k.CH_CH_C2_H2_cm3_7_44, 'CH_CH_C2_H2_cm3_7_44'); % 7.44
R = push(R, sto(species, 'CH', 1, 'C2H6', 1), sto(species, 'C3H6', 1, 'H', 1), k.CH_C2H6_C3H6_H_cm3_7_45, 'CH_C2H6_C3H6_H_cm3_7_45'); % 7.45
R = push(R, sto(species, 'CH', 1, 'C2H4', 1), sto(species, 'C2H2', 1, 'CH3', 1), k.CH_C2H4_C2H2_CH3_cm3_7_46, 'CH_C2H4_C2H2_CH3_cm3_7_46'); % 7.46
R = push(R, sto(species, 'C2H', 1, 'H', 1), sto(species, 'C2', 1, 'H2', 1), k.C2H_H_C2_H2_cm3_7_47, 'C2H_H_C2_H2_cm3_7_47'); % 7.47
R = push(R, sto(species, 'CH', 1, 'CH2', 1), sto(species, 'C2H2', 1, 'H', 1), k.CH_CH2_C2H2_H_cm3_7_48, 'CH_CH2_C2H2_H_cm3_7_48'); % 7.48
R = push(R, sto(species, 'CH3', 2), sto(species, 'C2H2', 1, 'H2', 2), k.CH3_CH3_C2H2_H2_H2_cm3_7_49, 'CH3_CH3_C2H2_H2_H2_cm3_7_49'); % 7.49
R = push(R, sto(species, 'C2H2', 1, 'H', 1), sto(species, 'C2', 1, 'H2', 1, 'H', 1), k.C2H2_H_C2_H2_H_cm3_7_50, 'C2H2_H_C2_H2_H_cm3_7_50'); % 7.50
R = push(R, sto(species, 'CH', 1, 'H2', 1), sto(species, 'CH2', 1, 'H', 1), k.CH_H2_CH2_H_cm3_7_51, 'CH_H2_CH2_H_cm3_7_51'); % 7.51
R = push(R, sto(species, 'C2', 1, 'CH', 1), sto(species, 'C3', 1, 'H', 1), k.C2_CH_C3_H_cm3_7_52, 'C2_CH_C3_H_cm3_7_52'); % 7.52
R = push(R, sto(species, 'CH2', 1, 'C2H3', 1), sto(species, 'C2H2', 1, 'CH3', 1), k.CH2_C2H3_C2H2_CH3_cm3_7_53, 'CH2_C2H3_C2H2_CH3_cm3_7_53'); % 7.53
R = push(R, sto(species, 'C', 1, 'C2H3', 1), sto(species, 'C2', 1, 'CH3', 1), k.C_C2H3_C2_CH3_cm3_7_54, 'C_C2H3_C2_CH3_cm3_7_54'); % 7.54
R = push(R, sto(species, 'CH', 1, 'C2H5', 1), sto(species, 'C3H6', 1), k.CH_C2H5_C3H6_cm3_7_55, 'CH_C2H5_C3H6_cm3_7_55'); % 7.55
R = push(R, sto(species, 'C2H2', 1, 'CH', 1), sto(species, 'C3', 1, 'H2', 1, 'H', 1), k.C2H2_CH_C3_H2_H_cm3_7_56, 'C2H2_CH_C3_H2_H_cm3_7_56'); % 7.56
R = push(R, sto(species, 'CH', 1, 'C2H6', 1), sto(species, 'C2H2', 1, 'CH3', 1, 'H', 1), k.CH_C2H6_C2H2_CH3_H_cm3_7_57, 'CH_C2H6_C2H2_CH3_H_cm3_7_57'); % 7.57
R = push(R, sto(species, 'CH2', 2), sto(species, 'C2', 1, 'H2', 2), k.CH2_CH2_C2_H2_H2_cm3_7_58, 'CH2_CH2_C2_H2_H2_cm3_7_58'); % 7.58
R = push(R, sto(species, 'CH', 1, 'CH4', 1), sto(species, 'C2H4', 1, 'H', 1), k.CH_CH4_C2H4_H_cm3_7_59, 'CH_CH4_C2H4_H_cm3_7_59'); % 7.59
R = push(R, sto(species, 'CH', 1, 'C2H4', 1), sto(species, 'C3H4', 1, 'H', 1), k.CH_C2H4_C3H4_H_cm3_7_60, 'CH_C2H4_C3H4_H_cm3_7_60'); % 7.60
R = push(R, sto(species, 'CH3', 1, 'C2H5', 1), sto(species, 'C2H2', 1, 'CH3', 1, 'H2', 1), k.CH3_C2H5_C2H2_CH3_H2_cm3_7_61, 'CH3_C2H5_C2H2_CH3_H2_cm3_7_61'); % 7.61
R = push(R, sto(species, 'CH2', 1, 'CH3', 1), sto(species, 'C2H2', 1, 'H', 1, 'H2', 1), k.CH2_CH3_C2H2_H_H2_cm3_7_62, 'CH2_CH3_C2H2_H_H2_cm3_7_62'); % 7.62
R = push(R, sto(species, 'CH2', 1, 'C2H5', 1), sto(species, 'C2H2', 1, 'CH3', 1, 'H', 1), k.CH2_C2H5_C2H2_CH3_H_cm3_7_63, 'CH2_C2H5_C2H2_CH3_H_cm3_7_63'); % 7.63

% Group 8: Termolecular Recombination
R = push(R, sto(species, 'H', 2), sto(species, 'H2', 1), k.H_H_M_H2_M_cm6_8_1, 'H_H_M_H2_M_cm6_8_1'); % 8.1
R = push(R, sto(species, 'CH3', 2), sto(species, 'C2H6', 1), k.CH3_CH3_M_C2H6_M_cm6_8_2, 'CH3_CH3_M_C2H6_M_cm6_8_2'); % 8.2
% FLAG: 8.2 Third-body efficiencies estimated. Needs validation.

% Group 9: Stick Reactions
R = push(R, sto(species, 'H', 1), sto(species), k.stick_H_9_1, 'stick_H_9_1'); % 9.1
R = push(R, sto(species, 'CH3', 1), sto(species), k.stick_CH3_9_2, 'stick_CH3_9_2'); % 9.2
R = push(R, sto(species, 'CH', 1), sto(species), k.stick_CH_9_3, 'stick_CH_9_3'); % 9.3
R = push(R, sto(species, 'ArPlus', 1), sto(species), k.stick_ArPlus_9_4, 'stick_ArPlus_9_4'); % 9.4
R = push(R, sto(species, 'ArStar', 1), sto(species), k.stick_ArStar_9_5, 'stick_ArStar_9_5'); % 9.5
R = push(R, sto(species, 'CH3Plus', 1), sto(species), k.stick_CH3Plus_9_6, 'stick_CH3Plus_9_6'); % 9.6
R = push(R, sto(species, 'CH5Plus', 1), sto(species), k.stick_CH5Plus_9_7, 'stick_CH5Plus_9_7'); % 9.7
R = push(R, sto(species, 'ArHPlus', 1), sto(species), k.stick_ArHPlus_9_8, 'stick_ArHPlus_9_8'); % 9.8
R = push(R, sto(species, 'C2', 1), sto(species), k.stick_C2_9_9, 'stick_C2_9_9'); % 9.9
R = push(R, sto(species, 'C', 1), sto(species), k.stick_C_9_10, 'stick_C_9_10'); % 9.10
R = push(R, sto(species, 'C2H2', 1), sto(species), k.stick_C2H2_9_11, 'stick_C2H2_9_11'); % 9.11
R = push(R, sto(species, 'C2H4', 1), sto(species), k.stick_C2H4_9_12, 'stick_C2H4_9_12'); % 9.12
R = push(R, sto(species, 'CH2', 1), sto(species), k.stick_CH2_9_13, 'stick_CH2_9_13'); % 9.13
R = push(R, sto(species, 'C2H6', 1), sto(species), k.stick_C2H6_9_14, 'stick_C2H6_9_14'); % 9.14
R = push(R, sto(species, 'CH3Minus', 1), sto(species), k.stick_CH3Minus_9_15, 'stick_CH3Minus_9_15'); % 9.15
R = push(R, sto(species, 'H2', 1), sto(species), k.stick_H2_9_16, 'stick_H2_9_16'); % 9.16
R = push(R, sto(species, 'C2H5', 1), sto(species), k.stick_C2H5_9_17, 'stick_C2H5_9_17'); % 9.17
R = push(R, sto(species, 'HMinus', 1), sto(species), k.stick_HMinus_9_18, 'stick_HMinus_9_18'); % 9.18
R = push(R, sto(species, 'C3H', 1), sto(species), k.stick_C3H_9_19, 'stick_C3H_9_19'); % 9.19
R = push(R, sto(species, 'C4H2', 1), sto(species), k.stick_C4H2_9_20, 'stick_C4H2_9_20'); % 9.20
R = push(R, sto(species, 'C3H3', 1), sto(species), k.stick_C3H3_9_21, 'stick_C3H3_9_21'); % 9.21
R = push(R, sto(species, 'C3H4', 1), sto(species), k.stick_C3H4_9_22, 'stick_C3H4_9_22'); % 9.22
R = push(R, sto(species, 'C3', 1), sto(species), k.stick_C3_9_23, 'stick_C3_9_23'); % 9.23
R = push(R, sto(species, 'C4H', 1), sto(species), k.stick_C4H_9_24, 'stick_C4H_9_24'); % 9.24
R = push(R, sto(species, 'C3H6', 1), sto(species), k.stick_C3H6_9_25, 'stick_C3H6_9_25'); % 9.25
R = push(R, sto(species, 'H3Plus', 1), sto(species), k.stick_H3Plus_9_26, 'stick_H3Plus_9_26'); % 9.26
R = push(R, sto(species, 'CHPlus', 1), sto(species), k.stick_CHPlus_9_27, 'stick_CHPlus_9_27'); % 9.27
R = push(R, sto(species, 'C2HPlus', 1), sto(species), k.stick_C2HPlus_9_28, 'stick_C2HPlus_9_28'); % 9.28

% Group 10: Drift Losses
R = push(R, sto(species, 'ArPlus', 1), sto(species), k.drift_ArPlus_10_1, 'drift_ArPlus_10_1'); % 10.1
R = push(R, sto(species, 'CH4Plus', 1), sto(species), k.drift_CH4Plus_10_2, 'drift_CH4Plus_10_2'); % 10.2
R = push(R, sto(species, 'CH3Plus', 1), sto(species), k.drift_CH3Plus_10_3, 'drift_CH3Plus_10_3'); % 10.3
R = push(R, sto(species, 'CH5Plus', 1), sto(species), k.drift_CH5Plus_10_4, 'drift_CH5Plus_10_4'); % 10.4
R = push(R, sto(species, 'ArHPlus', 1), sto(species), k.drift_ArHPlus_10_5, 'drift_ArHPlus_10_5'); % 10.5
R = push(R, sto(species, 'CH2Plus', 1), sto(species), k.drift_CH2Plus_10_6, 'drift_CH2Plus_10_6'); % 10.6
R = push(R, sto(species, 'C2H5Plus', 1), sto(species), k.drift_C2H5Plus_10_7, 'drift_C2H5Plus_10_7'); % 10.7
R = push(R, sto(species, 'C2H4Plus', 1), sto(species), k.drift_C2H4Plus_10_8, 'drift_C2H4Plus_10_8'); % 10.8
R = push(R, sto(species, 'C2H3Plus', 1), sto(species), k.drift_C2H3Plus_10_9, 'drift_C2H3Plus_10_9'); % 10.9
R = push(R, sto(species, 'H3Plus', 1), sto(species), k.drift_H3Plus_10_10, 'drift_H3Plus_10_10'); % 10.10
R = push(R, sto(species, 'CHPlus', 1), sto(species), k.drift_CHPlus_10_11, 'drift_CHPlus_10_11'); % 10.11
R = push(R, sto(species, 'CH3Minus', 1), sto(species), k.drift_CH3Minus_10_12, 'drift_CH3Minus_10_12'); % 10.12
R = push(R, sto(species, 'C2HPlus', 1), sto(species), k.drift_C2HPlus_10_13, 'drift_C2HPlus_10_13'); % 10.13

% Group 11: Loss Reactions
R = push(R, sto(species, 'CH2', 1), sto(species), k.loss_CH2_11_1, 'loss_CH2_11_1'); % 11.1
R = push(R, sto(species, 'H2', 1), sto(species), k.loss_H2_11_2, 'loss_H2_11_2'); % 11.2
R = push(R, sto(species, 'C2', 1), sto(species), k.loss_C2_11_3, 'loss_C2_11_3'); % 11.3
R = push(R, sto(species, 'e', 1), sto(species), k.loss_e_11_4, 'loss_e_11_4'); % 11.4
R = push(R, sto(species, 'C2H6', 1), sto(species), k.loss_C2H6_11_5, 'loss_C2H6_11_5'); % 11.5
R = push(R, sto(species, 'CH4', 1), sto(species), k.loss_CH4_11_6, 'loss_CH4_11_6'); % 11.6
R = push(R, sto(species, 'Ar', 1), sto(species), k.loss_Ar_11_7, 'loss_Ar_11_7'); % 11.7
% FLAG: 11.7 Rate is estimated. Needs validation.
R = push(R, sto(species, 'C', 1), sto(species), k.loss_C_11_8, 'loss_C_11_8'); % 11.8
R = push(R, sto(species, 'CH', 1), sto(species), k.loss_CH_11_9, 'loss_CH_11_9'); % 11.9
% FLAG: 11.9 Rate is estimated. Needs validation.
R = push(R, sto(species, 'C3H', 1), sto(species), k.loss_C3H_11_10, 'loss_C3H_11_10'); % 11.10
R = push(R, sto(species, 'C4H2', 1), sto(species), k.loss_C4H2_11_11, 'loss_C4H2_11_11'); % 11.11
R = push(R, sto(species, 'C2H', 1), sto(species), k.loss_C2H_11_12, 'loss_C2H_11_12'); % 11.12
R = push(R, sto(species, 'C3H3', 1), sto(species), k.loss_C3H3_11_13, 'loss_C3H3_11_13'); % 11.13
R = push(R, sto(species, 'C3', 1), sto(species), k.loss_C3_11_14, 'loss_C3_11_14'); % 11.14
R = push(R, sto(species, 'C4H', 1), sto(species), k.loss_C4H_11_15, 'loss_C4H_11_15'); % 11.15
R = push(R, sto(species, 'C3H6', 1), sto(species), k.loss_C3H6_11_16, 'loss_C3H6_11_16'); % 11.16
R = push(R, sto(species, 'C2H5', 1), sto(species), k.loss_C2H5_11_17, 'loss_C2H5_11_17'); % 11.17
R = push(R, sto(species, 'C3H4', 1), sto(species), k.loss_C3H4_11_18, 'loss_C3H4_11_18'); % 11.18
R = push(R, sto(species, 'C2H2', 1), sto(species), k.loss_C2H2_11_19, 'loss_C2H2_11_19'); % 11.19
R = push(R, sto(species, 'C2H4', 1), sto(species), k.loss_C2H4_11_20, 'loss_C2H4_11_20'); % 11.20
R = push(R, sto(species, 'CH3', 1), sto(species), k.loss_CH3_11_21, 'loss_CH3_11_21'); % 11.21
R = push(R, sto(species, 'C2H3', 1), sto(species), k.loss_C2H3_11_22, 'loss_C2H3_11_22'); % 11.22
R = push(R, sto(species, 'C3H2', 1), sto(species), k.loss_C3H2_11_23, 'loss_C3H2_11_23'); % 11.23
R = push(R, sto(species, 'C3H5', 1), sto(species), k.loss_C3H5_11_24, 'loss_C3H5_11_24'); % 11.24
R = push(R, sto(species, 'C2H2Star', 1), sto(species, 'C2H2', 1), k.loss_C2H2Star_11_25, 'loss_C2H2Star_11_25'); % 11.25


for i = 1:length(R)
    if ~isfield(k, R(i).tag)
        error('Rate constant for reaction %s not found in params.k', R(i).tag);
    end
end
end