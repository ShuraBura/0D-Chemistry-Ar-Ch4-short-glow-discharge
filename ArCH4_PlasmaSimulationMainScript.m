% Ar/CH4 Plasma Simulation Main Script
% 0-D model for short glow discharge. Conditions are: P=400 mTorr. T=400 K.  Gas - Ar 85% ,  CH4 15% . 
% The anode, 20 mm diameter , brass , at the top. V=1000 V is applied to anode
% Cathode is at the bottom, 40 mm diameter, copper, grounded.
% Cathode-anode gap is 17 mm. Measured current is 0.736 mA, cross section A for current a circle with a diameter of 36 mm. 
% Ne seeded according to the measured current.  Ar and CH4 densities are also held constant ( constant gas flow). 
% The model uses the chemical network to calculate densities of all species and ions. 
% The rates for all reactions are given in define_rates file citing sources and acceptable ranges.
% The rates can be tuned in this literature range only.  
% We are modeling ~4.5 mm from cathode, where measured nH is ~5.18e13, nCH~1e9 cm, nC2~1.3e11.
% There is a contribution from H drifting from its peak (nH=2.4e14 cm^-3) at cathode glow (0.5 mm from cathode) introduced as H_drift_gain term in odefun_diagnostics. 
clc
clear all;
close all;
% Clear ch2_diagnostics.txt
if exist('ch2_diagnostics.txt', 'file')
    delete('ch2_diagnostics.txt');
end
tic
params = struct();
params.P = 0.4; % Pressure in Torr
params.Tg = 400; % Gas temperature in K
params.n_tot = 9.66e15; % Total density in cm^-3
params.ne = 1e10; % Electron density in cm^-3
params.Te = 1; % Electron temperature in eV
params.tspan = [0, 0.1, 1, 10, 100]; % Simulation time span
params.verbose = 1; % Enable debug output
params.E_field = 50; % V/cm;
params.L_discharge = 0.45; % cm
params.species = {'e', 'Ar', 'CH4', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', ...
                  'CH3Minus', 'H', 'C2', 'CH', 'H2', 'ArStar', 'C2H4', 'C2H6', 'CH2', ...
                  'C2H2', 'C2H5', 'CH3', 'C', 'H3Plus', 'C2H3', 'C3H2', 'CHPlus', 'C3H', ...
                  'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6', 'CH2Plus',...
                  'C2H5Plus','C2H4Plus','C2H3Plus','HMinus','C2HPlus','C2H2Star'};

params.ion_species = {'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH3Minus', 'H3Plus', 'CHPlus', ...
                       'CH2Plus', 'C2H5Plus','C2H4Plus','C2H3Plus','HMinus','C2HPlus'};

params.mobilities = struct(...
    'ArPlus', 3057.28, ...
    'CH4Plus', 6432, ... % Modified: 5360 * 1.2 to achieve k_drift ~54.55 s^-1
    'CH3Plus', 4949.6, ...
    'CH5Plus', 4761.6, ... % Unchanged: Already achieves k_drift ~48 s^-1
    'ArHPlus', 2969.6, ...
    'CH2Plus', 4949.6, ...
    'C2H5Plus', 4949.6, ...
    'C2H4Plus', 4949.6, ...
    'C2H3Plus', 4949.6, ...
    'C2HPlus', 5000, ...
    'H3Plus', 5000, ... % Estimated, similar to CH3Plus, NIST
    'CHPlus', 5000, ... % Estimated, similar to CH3Plus, NIST
    'CH3Minus', 3000, ... % Estimated, similar to ArHPlus, NIST
    'HMinus', 3000 ... % Estimated, similar to ArHPlus, NIST
);

params.k = define_rates(params);
[params.R, params.tags] = build_reactions(params);
y0 = calculate_initial_densities(params);

AbsTol = 1e-8 * ones(length(params.species), 1); % Base tolerance
% Looser tolerance for high-density species
high_density_species = {'Ar', 'CH4'};
for sp = high_density_species
    idx = find(strcmp(params.species, sp{1}));
    if ~isempty(idx)
        AbsTol(idx) = 1e-6; % High-density species
    end
end
% Tighter tolerance for low-density species
low_density_species = {'ArStar', 'C2H5Plus', 'CH'};
for sp = low_density_species
    idx = find(strcmp(params.species, sp{1}));
    if ~isempty(idx)
        AbsTol(idx) = 1e-10; % Low-density species
    end
end
options = odeset(...
    'RelTol', 1e-7, ... % Tighter for accuracy
    'AbsTol', AbsTol, ... % Species-specific tolerances
    'JPattern', build_JPattern(params.species, params.R, params.ion_species), ...
    'OutputFcn', @(t,y,flag) progress_output_func(t, y, flag, params), ...
    'NonNegative', 1:length(params.species), ...
    'Stats', 'on', ...
    'OutputSel', 1:length(params.species), ...
    'Refine', 1, ...
    'InitialStep', 1e-6 ... % Stable initial step
);
% Rationale: RelTol = 1e-7, AbsTol = 1e-10 for ArStar, C2H5Plus, CH improve rate accuracy.

% Rationale: RelTol = 1e-7, AbsTol = 1e-10 for ArStar, C2H5Plus, CH improve rate accuracy.
[t, y] = ode15s(@(t,y) odefun_diagnostic(t, y, params), params.tspan, y0, options);

toc

% Extract final densities for all species
disp(' ');
disp('--- Steady-State Densities at Final Time (t = 100 s) ---');

% Neutral Species
neutral_species_ordered = {'H', 'CH', 'C2', 'H2', 'C2H2', 'CH3', 'CH2', 'Ar', 'CH4', 'ArStar', 'C2H4', 'C2H6', 'C2H5', 'C', 'C2H3', 'C3H2', 'C3H', 'C4H2', 'C2H', 'C3H3', 'C3H4', 'C3', 'C3H5', 'C4H', 'C3H6'};
disp('Neutral Species Densities (cm^-3):');
for species = neutral_species_ordered
    idx = find(strcmp(params.species, species{1}));
    if ~isempty(idx)
        fprintf('%s: %.4e\n', species{1}, y(end, idx));
    end
end

% Ion Species
ion_species = {'e', 'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', 'CH3Minus', 'H3Plus', 'CHPlus', 'CH2Plus', 'C2H5Plus', 'C2H4Plus', 'C2H3Plus', 'HMinus', 'C2HPlus'};
disp('Ion Species Densities (cm^-3):');
for species = ion_species
    idx = find(strcmp(params.species, species{1}));
    if ~isempty(idx)
        fprintf('%s: %.4e\n', species{1}, y(end, idx));
    end
end

% Key Species and Charge Balance
disp(' ');
disp('--- Key Species and Charge Balance ---');
final_H_density = y(end, strcmp(params.species, 'H'));
final_CH_density = y(end, strcmp(params.species, 'CH'));
final_C2_density = y(end, strcmp(params.species, 'C2'));
final_e_density = y(end, strcmp(params.species, 'e'));
[~, final_diagnostics] = odefun_diagnostic(t(end), y(end,:)', params);
final_n_positive_ions_sum = final_diagnostics.charge_balance.n_positive_ions_sum;
final_n_negative_ions_sum = final_diagnostics.charge_balance.n_negative_ions_sum;
n_i_net = final_diagnostics.charge_balance.n_i_net;
fprintf('Final H density: %.4e cm^-3\n', final_H_density);
fprintf('Final CH density: %.4e cm^-3\n', final_CH_density);
fprintf('Final C2 density: %.4e cm^-3\n', final_C2_density);
fprintf('Final Electron density (fixed): %.4e cm^-3\n', final_e_density);
fprintf('Final Total Positive Ion density: %.4e cm^-3\n', final_n_positive_ions_sum);
fprintf('Final Total Negative Ion density (CH3Minus + HMinus): %.4e cm^-3\n', final_n_negative_ions_sum);
fprintf('Net Ion Density (n_i_net): %.4e cm^-3\n', n_i_net);
fprintf('Charge Neutrality Check (n_i_net / ne): %.4f\n', n_i_net / final_e_density);

disp(' ');
disp('--- Diagnostic Rates for H (at final time) ---');
fprintf('H Total Production Rate: %.4e cm^-3 s^-1\n', final_diagnostics.H.production_total);
fprintf('H Total Consumption Rate: %.4e cm^-3 s^-1\n', final_diagnostics.H.consumption_total);

disp(' ');
disp('--- Diagnostic Rates for CH (at final time) ---');
fprintf('CH Total Production Rate: %.4e cm^-3 s^-1\n', final_diagnostics.CH.production_total);
fprintf('CH Total Consumption Rate: %.4e cm^-3 s^-1\n', final_diagnostics.CH.consumption_total);

disp(' ');
disp('--- Diagnostic Rates for C2 (at final time) ---');
fprintf('C2 Total Production Rate: %.4e cm^-3 s^-1\n', final_diagnostics.C2.production_total);
fprintf('C2 Total Consumption Rate: %.4e cm^-3 s^-1\n', final_diagnostics.C2.consumption_total);


function status = progress_output_func(t, y, flag, ~)
    status = 0; % Continue integration
    if strcmp(flag, 'done')
        fprintf('\nIntegration finished.\n');
    end
    % Diagnostics handled by odefun_diagnostic
end

function S = build_JPattern(species, R, ion_species)
    ns = length(species); % Number of species (70)
    S = sparse(ns, ns); % Initialize sparse matrix
    
    % Loop through reactions
    for r = 1:length(R)
        % Get indices of reactants and products with non-zero stoichiometry
        reactants_idx = find(R(r).reactants > 0);
        products_idx = find(R(r).products > 0);
        
        % Vectorized: Product species i depend on all reactants j
        if ~isempty(products_idx) && ~isempty(reactants_idx)
            [i, j] = meshgrid(products_idx, reactants_idx);
            S(sub2ind([ns, ns], i(:), j(:))) = 1;
        end
        
        % Vectorized: Reactant species i depend only on reactants j in rate
        if ~isempty(reactants_idx)
            [i, j] = meshgrid(reactants_idx, reactants_idx);
            S(sub2ind([ns, ns], i(:), j(:))) = 1;
        end
    end    
    % Add ion cap dependencies: all ions (including CH3Minus) affect each other
    ion_indices = zeros(length(ion_species) + 1, 1);
    for s = 1:length(ion_species)
        ion_indices(s) = find(strcmp(species, ion_species{s}));
    end
    ion_indices(end) = find(strcmp(species, 'CH3Minus'));
    [i, j] = meshgrid(ion_indices, ion_indices);
    S(sub2ind([ns, ns], i(:), j(:))) = 1;    
    % Ensure diagonal (self-dependencies)
    S = S + speye(ns);    
    % Set rows for fixed species (e, Ar, CH4) to zero, as dydt = 0
    fixed_indices = [find(strcmp(species, 'e')), ...
                     find(strcmp(species, 'Ar')), ...
                     find(strcmp(species, 'CH4'))];
    S(fixed_indices, :) = 0;
end

function y0 = calculate_initial_densities(params)
    species = params.species;
    ns = length(species);
    y0 = zeros(ns, 1);    
    % Set specific initial density
    y0(strcmp(species, 'e')) = params.ne; % Electron density at CDS center (cm^-3)
    y0(strcmp(species, 'Ar')) = 0.85 * 9.66e15; % 85% of n_tot
    y0(strcmp(species, 'CH4')) = 0.15 * 9.66e15; % 15% of n_tot
    y0(strcmp(species, 'ArPlus')) = 1e7; % 
    y0(strcmp(species, 'CH4Plus')) = 1e5; % 
    y0(strcmp(species, 'CH3Plus')) = 1e5; % 
    y0(strcmp(species, 'CH5Plus')) = 1e3; % 
    y0(strcmp(species, 'ArHPlus')) = 5e5; % Reduced contribution
    y0(strcmp(species, 'CH3Minus')) = 5e4; % Minimized to reduce negative ions
    y0(strcmp(species, 'H2')) = 1e12; % Lower end of target [1e13, 1e15]
    y0(strcmp(species, 'ArStar')) = 5e6; % Within target
    y0(strcmp(species, 'H')) = 1e11; % Below target
    y0(strcmp(species, 'C2')) = 5e7; % Below target
    y0(strcmp(species, 'CH')) = 5e4; % Near target
    y0(strcmp(species, 'C2H4')) = 5e7;
    y0(strcmp(species, 'C2H6')) = 1e6;
    y0(strcmp(species, 'CH2')) = 1e11;
    y0(strcmp(species, 'C2H2')) = 1e12;
    y0(strcmp(species, 'C2H5')) = 1e6;
    y0(strcmp(species, 'CH3')) = 5e7;
    y0(strcmp(species, 'C')) = 5e7;
    y0(strcmp(species, 'H3Plus')) = 1e4; % Reduced contribution
    y0(strcmp(species, 'HMinus')) = 1e4; % Low initial density, similar to CH3Minus
    y0(strcmp(species, 'C2HPlus')) = 1e5; % Similar to CH3Plus
    y0(strcmp(species, 'CHPlus')) = 1e4; % Similar to CH3Plus
    y0(strcmp(species, 'C2H6')) = 1e6; % Keep or adjust to 1e7 if needed
    y0(strcmp(species, 'C3H6')) = 1e3; % Keep or adjust to 1e4
    y0(strcmp(species, 'C2H2Star')) = 1e6; % cm^-3
    y0 = max(y0, 1e3); % Allow CH3Minus at 5e3
end