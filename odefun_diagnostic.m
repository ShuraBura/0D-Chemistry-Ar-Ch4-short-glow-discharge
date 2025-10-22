% ODE function for Ar/CH4 plasma simulation with diagnostics and enhanced debugging
% Goals: 
% - Limit debug logs for 3.22, 3.23, 6.28, 7.62, 7.63 to t = [0, 0.1, 1, 10, 100] s
% - Log density for ArStar, C2H2, CH3, C2H5Plus, CH2, C2H5 to diagnose rate scaling
% - Focus diagnostics on out-of-range species (H, CH, C2, C2H2, ArStar)
% - Debug rate calculations and dydt contributions to resolve ~1000x and ~10x rate discrepancies
% Notes: Uses persistent debug_printed_times, time_tolerance = 1e-5, includes 9.11, 11.19 debugging

function [dydt, diagnostics, varargout] = odefun_diagnostic(t, y, params)
    % ODE function for Ar/CH4 plasma simulation with corrected diagnostics.
    species = params.species;
    ns = length(species);
    R = params.R;
    nr = length(R);
    k = params.k;
    tags = params.tags;
    H_drift_gain = 3.2e17; % Increased to boost H density to ~5.18e13 cm^-3
    % Initialize outputs
    dydt = zeros(ns, 1);
    diagnostics = struct();

    % Open diagnostics file early
    try, fid = fopen('diagnostics.txt', 'a'); catch, fid = -1; end

    % Hardcode species indices
    e_idx = find(strcmp(species, 'e'));
    Ar_idx = find(strcmp(species, 'Ar'));
    CH4_idx = find(strcmp(species, 'CH4'));
    
    % Prevent negative density
    y(y<0) = 1e-6;

    % Compute reaction rates
    reaction_rates_values = zeros(nr, 1);
    target_times = [0, 0.1, 1, 10, 100];
    time_tolerance = 1e-5; % Stricter tolerance
    
    % Persistent flag for debug logging
    persistent debug_printed_times;
    if isempty(debug_printed_times), debug_printed_times = []; end
    is_target_time = any(abs(t - target_times) < time_tolerance) && ~any(abs(debug_printed_times - t) < time_tolerance);
    
    % Debug density for key species
    if is_target_time
        fprintf('--- Density Debug Triggered at t = %.3e s ---\n', t);
        if fid ~= -1
            fprintf(fid, '--- Density Debug Triggered at t = %.3e s ---\n', t);
        end
        species_to_debug = {'ArStar', 'C2H2', 'CH3', 'C2H5Plus', 'CH2', 'C2H5'};
        for s = species_to_debug
            idx = find(strcmp(species, s{1}));
            if ~isempty(idx)
                fprintf('Debug Density %s t=%.3e: y=%e cm^-3\n', s{1}, t, y(idx));
                if fid ~= -1
                    fprintf(fid, 'Debug Density %s t=%.3e: y=%e cm^-3\n', s{1}, t, y(idx));
                end
            end
        end
    end

    for i = 1:nr
        current_rate = k.(tags{i});
        reactant_indices = find(R(i).reactants > 0);
        for j = 1:length(reactant_indices)
            species_idx = reactant_indices(j);
            stoich_coeff = R(i).reactants(species_idx);
            current_rate = current_rate * (y(species_idx)^stoich_coeff);
        end
        if is_target_time && any(strcmp(tags{i}, {'ArStar_C2H2_Ar_C2H2Star_cm3_3_22', 'ArStar_CH3_Ar_CH2_H_cm3_3_23', 'C2H5Plus_e_C2H4_H_cm3_6_28', 'CH2_CH3_C2H2_H_H2_cm3_7_62', 'CH2_C2H5_C2H2_CH3_H_cm3_7_63'}))
            fprintf('Debug %s t=%.3e: k=%e, reactants=%s, rate=%e\n', tags{i}, t, k.(tags{i}), mat2str(y(reactant_indices)), current_rate);
            if fid ~= -1
                fprintf(fid, 'Debug %s t=%.3e: k=%e, reactants=%s, rate=%e\n', tags{i}, t, k.(tags{i}), mat2str(y(reactant_indices)), current_rate);
            end
        end
        reaction_rates_values(i) = current_rate;
    end

    % Enhanced rate debugging
    if is_target_time
        fprintf('--- Rate Debug at t = %.3e s ---\n', t);
        if fid ~= -1
            fprintf(fid, '--- Rate Debug at t = %.3e s ---\n', t);
        end
        focus_reactions = {'ArStar_C2H2_Ar_C2H2Star_cm3_3_22', 'ArStar_CH3_Ar_CH2_H_cm3_3_23', ...
                           'C2H5Plus_e_C2H4_H_cm3_6_28', 'CH2_CH3_C2H2_H_H2_cm3_7_62', ...
                           'CH2_C2H5_C2H2_CH3_H_cm3_7_63', 'stick_C2H2_9_11', 'loss_C2H2_11_19'};
        for i = 1:nr
            reactant_indices = find(R(i).reactants > 0);
            expected_rate = k.(tags{i});
            reactant_str = '';
            for j = 1:length(reactant_indices)
                species_idx = reactant_indices(j);
                stoich_coeff = R(i).reactants(species_idx);
                expected_rate = expected_rate * (y(species_idx)^stoich_coeff);
                reactant_str = [reactant_str, sprintf('%s^%d (%.3e) ', species{species_idx}, stoich_coeff, y(species_idx))];
            end
            ratio = reaction_rates_values(i) / max(1e-30, expected_rate);
            if abs(ratio - 1) > 1e-6 || any(strcmp(tags{i}, focus_reactions))
                fprintf('Rate Debug %s: computed=%e, expected=%e, ratio=%e, reactants=%s\n', ...
                    tags{i}, reaction_rates_values(i), expected_rate, ratio, reactant_str);
                if fid ~= -1
                    fprintf(fid, 'Rate Debug %s: computed=%e, expected=%e, ratio=%e, reactants=%s\n', ...
                        tags{i}, reaction_rates_values(i), expected_rate, ratio, reactant_str);
                end
            end
        end
    end

    if is_target_time
        debug_printed_times = [debug_printed_times, t];
    end

    % Compute density changes (dydt)
    for i = 1:ns
        for j = 1:nr
            net_coeffs = R(j).products(i) - R(j).reactants(i);
            dydt(i) = dydt(i) + net_coeffs * reaction_rates_values(j);
        end
    end

    % Debug dydt contributions
    if is_target_time
        species_to_check = {'ArStar', 'C2H5Plus', 'CH2', 'C2H2'};
        fprintf('--- dydt Debug at t = %.3e s ---\n', t);
        if fid ~= -1
            fprintf(fid, '--- dydt Debug at t = %.3e s ---\n', t);
        end
        for s = species_to_check
            idx = find(strcmp(species, s{1}));
            if ~isempty(idx)
                total_rate = 0;
                for j = 1:nr
                    net_coeff = R(j).products(idx) - R(j).reactants(idx);
                    total_rate = total_rate + net_coeff * reaction_rates_values(j);
                end
                fprintf('dydt Debug %s: dydt=%e, sum_rates=%e\n', s{1}, dydt(idx), total_rate);
                if fid ~= -1
                    fprintf(fid, 'dydt Debug %s: dydt=%e, sum_rates=%e\n', s{1}, dydt(idx), total_rate);
                end
            end
        end
    end

    % Add diffused H source term
    H_idx_for_gain = find(strcmp(species, 'H'));
    if ~isempty(H_idx_for_gain)
        dydt(H_idx_for_gain) = dydt(H_idx_for_gain) + H_drift_gain;
    end

    % Ensure fixed species' derivatives are zero
    if ~isempty(e_idx), dydt(e_idx) = 0; end
    if ~isempty(Ar_idx), dydt(Ar_idx) = 0; end
    if ~isempty(CH4_idx), dydt(CH4_idx) = 0; end

    % --- DETAILED DIAGNOSTICS ---
    species_to_diagnose = {'H', 'CH', 'C2', 'C2H2', 'ArStar'};
    diagnostic_threshold_pct = 5.0;

    function pct = get_pct(val, total)
        if total > 0, pct = 100 * val / total; else, pct = 0; end
    end

    persistent printed_times;
    if isempty(printed_times), printed_times = []; end

    % Specific reactions to log
    reactions_to_log = {...
        'ArStar_C2H2_Ar_C2H2Star_cm3_3_22', ...
        'ArStar_CH3_Ar_CH2_H_cm3_3_23', ...
        'CH3_C2H5_C2H2_CH3_H2_cm3_7_61', ...
        'CH2_CH3_C2H2_H_H2_cm3_7_62', ...
        'C2H5Plus_e_C2H4_H_cm3_6_28', ...
        'CH2_C2H5_C2H2_CH3_H_cm3_7_63'};
    
    for s_idx = 1:length(species_to_diagnose)
        species_name = species_to_diagnose{s_idx};
        species_idx = find(strcmp(species, species_name));
        if isempty(species_idx), continue; end
        
        total_prod = 0;
        total_cons = 0;
        if strcmp(species_name, 'H'), total_prod = total_prod + H_drift_gain; end

        prod_reactions = struct('tag', {}, 'rate', {}, 'percent', {});
        cons_reactions = struct('tag', {}, 'rate', {}, 'percent', {});
        specific_reactions = struct('tag', {}, 'rate', {}, 'percent', {}, 'type', {});

        for j = 1:nr
            net_coeff = R(j).products(species_idx) - R(j).reactants(species_idx);
            rate_val = reaction_rates_values(j);
            if net_coeff > 0
                total_prod = total_prod + net_coeff * rate_val;
                prod_reactions(end+1) = struct('tag', tags{j}, 'rate', net_coeff * rate_val, 'percent', 0);
                if any(strcmp(tags{j}, reactions_to_log))
                    specific_reactions(end+1) = struct('tag', tags{j}, 'rate', net_coeff * rate_val, 'percent', 0, 'type', 'Production');
                end
            elseif net_coeff < 0
                total_cons = total_cons + abs(net_coeff) * rate_val;
                cons_reactions(end+1) = struct('tag', tags{j}, 'rate', abs(net_coeff) * rate_val, 'percent', 0);
                if any(strcmp(tags{j}, reactions_to_log))
                    specific_reactions(end+1) = struct('tag', tags{j}, 'rate', abs(net_coeff) * rate_val, 'percent', 0, 'type', 'Consumption');
                end
            end
        end

        diagnostics.(species_name).production_total = total_prod;
        diagnostics.(species_name).consumption_total = total_cons;
        diagnostics.(species_name).density = y(species_idx);

        for k = 1:length(prod_reactions), prod_reactions(k).percent = get_pct(prod_reactions(k).rate, total_prod); end
        for k = 1:length(cons_reactions), cons_reactions(k).percent = get_pct(cons_reactions(k).rate, total_cons); end
        for k = 1:length(specific_reactions), specific_reactions(k).percent = get_pct(specific_reactions(k).rate, strcmp(specific_reactions(k).type, 'Production') * total_prod + strcmp(specific_reactions(k).type, 'Consumption') * total_cons); end
        
        for target_t = target_times
            if abs(t - target_t) < time_tolerance && ~any(abs(printed_times - target_t) < time_tolerance)
                fprintf('\n--- %s Diagnostics at t = %.3e s ---\n', species_name, t);
                fprintf('%s Density: %.3e cm^-3\n', species_name, y(species_idx));
                fprintf('Total Production Rate: %.3e cm^-3 s^-1\n', total_prod);
                fprintf('Total Consumption Rate: %.3e cm^-3 s^-1\n', total_cons);

                fprintf('\n%s Major Production Reactions (>= %.1f%%):\n', species_name, diagnostic_threshold_pct);
                fprintf('%-35s Rate (cm^-3 s^-1)\tPercent\n', 'Reaction Tag');
                found_major_prod = false;
                for k = 1:length(prod_reactions)
                    if prod_reactions(k).percent >= diagnostic_threshold_pct
                        fprintf('%-35s %.3e\t\t%.2f%%\n', prod_reactions(k).tag, prod_reactions(k).rate, prod_reactions(k).percent);
                        found_major_prod = true;
                    end
                end
                if ~found_major_prod, fprintf('None above threshold.\n'); end

                fprintf('\n%s Major Consumption Reactions (>= %.1f%%):\n', species_name, diagnostic_threshold_pct);
                fprintf('%-35s Rate (cm^-3 s^-1)\tPercent\n', 'Reaction Tag');
                found_major_cons = false;
                for k = 1:length(cons_reactions)
                    if cons_reactions(k).percent >= diagnostic_threshold_pct
                        fprintf('%-35s %.3e\t\t%.2f%%\n', cons_reactions(k).tag, cons_reactions(k).rate, cons_reactions(k).percent);
                        found_major_cons = true;
                    end
                end
                if ~found_major_cons, fprintf('None above threshold.\n'); end

                % Log specific reactions
                if any(strcmp(species_name, {'ArStar', 'C2H2', 'C2H5Plus'}))
                    fprintf('\n%s Specific Reaction Rates:\n', species_name);
                    fprintf('%-35s %-12s Rate (cm^-3 s^-1)\tPercent\n', 'Reaction Tag', 'Type');
                    found_specific = false;
                    for k = 1:length(specific_reactions)
                        fprintf('%-35s %-12s %.3e\t\t%.2f%%\n', specific_reactions(k).tag, specific_reactions(k).type, specific_reactions(k).rate, specific_reactions(k).percent);
                        found_specific = true;
                        if fid ~= -1
                            fprintf(fid, '%s t=%.3e %s %s Rate: %.3e cm^-3 s^-1 (%.2f%%)\n', species_name, t, specific_reactions(k).tag, specific_reactions(k).type, specific_reactions(k).rate, specific_reactions(k).percent);
                        end
                    end
                    if ~found_specific
                        fprintf('None of the specified reactions contribute.\n');
                        if fid ~= -1, fprintf(fid, '%s t=%.3e No specified reactions contribute.\n', species_name, t); end
                    end
                end

                % Log species indices for C2H2Star check
                if strcmp(species_name, 'C2H2')
                    c2h2star_idx = find(strcmp(species, 'C2H2Star'));
                    fprintf('\nC2H2 Species Index: %d, C2H2Star Index: %s\n', species_idx, mat2str(c2h2star_idx));
                    if fid ~= -1
                        fprintf(fid, 'C2H2 t=%.3e C2H2 Index: %d, C2H2Star Index: %s\n', t, species_idx, mat2str(c2h2star_idx));
                    end
                end

                fprintf('--- End %s Diagnostics ---\n', species_name);
                if fid ~= -1
                    fprintf(fid, '--- End %s Diagnostics at t=%.3e ---\n', species_name, t);
                end
            end
        end
    end

    for target_t = target_times
        if abs(t - target_t) < time_tolerance && ~any(abs(printed_times - target_t) < time_tolerance)
            printed_times = [printed_times, target_t];
            break; 
        end
    end

    % --- Charge balance ---
    positive_ion_list = {'ArPlus', 'CH4Plus', 'CH3Plus', 'CH5Plus', 'ArHPlus', ...
                         'H3Plus', 'CHPlus', 'CH2Plus', 'C2H5Plus','C2H4Plus', ...
                         'C2H3Plus','C2HPlus'};
    
    negative_ion_list = {'CH3Minus','HMinus'};

    n_positive_ions_sum = 0;
    for i = 1:length(positive_ion_list)
        ion_idx = find(strcmp(species, positive_ion_list{i}));
        if ~isempty(ion_idx)
            n_positive_ions_sum = n_positive_ions_sum + y(ion_idx);
        end
    end

    n_negative_ions_sum = 0;
    for i = 1:length(negative_ion_list)
        ion_idx = find(strcmp(species, negative_ion_list{i}));
        if ~isempty(ion_idx)
            n_negative_ions_sum = n_negative_ions_sum + y(ion_idx);
        end
    end
    
    diagnostics.charge_balance.n_e = y(e_idx);
    diagnostics.charge_balance.n_i_net = n_positive_ions_sum - n_negative_ions_sum;
    diagnostics.charge_balance.n_positive_ions_sum = n_positive_ions_sum;
    diagnostics.charge_balance.n_negative_ions_sum = n_negative_ions_sum;

    if fid ~= -1, fclose(fid); end
end