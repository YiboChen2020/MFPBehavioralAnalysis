function [subsections, other] = choose_opt_subdivision(neural_data, labels, other_label, n_subsections_vals_to_try, loss_parameters)


if isempty(loss_parameters)
    loss_parameters.max_duration = 0.25;
    loss_parameters.total_duration = 0.25;
    loss_parameters.purity = 0.25;
    loss_parameters.diversity = 0.25;
end


% Try each number of subsections (each value for n_subsections).
n_n_subsections_vals_to_try = length(n_subsections_vals_to_try);
subsections_trial = cell(n_n_subsections_vals_to_try, 1);
other_trial = cell(n_n_subsections_vals_to_try, 1);
utility = nan(n_n_subsections_vals_to_try, 1);
for i_n_subsections = 1:n_n_subsections_vals_to_try
    % Try current value for n_subsections.
    [subsections_trial{i_n_subsections}, other_trial{i_n_subsections}] ...
        = subdivide_session(neural_data, ...
                            labels, ...
                            other_label, ...
                            n_subsections_vals_to_try(i_n_subsections));

    % Print results to console.
    if true
        fprintf('\nn_subsections = %d \n', n_subsections_vals_to_try(i_n_subsections))
        for i_subsection = 1:n_subsections_vals_to_try(i_n_subsections)
            fprintf('duration = %d. purity = %f. diversity = %d bhvs. \n', ...
                    length(subsections_trial{i_n_subsections}{i_subsection}.labels), ...
                    subsections_trial{i_n_subsections}{i_subsection}.purity, ...
                    subsections_trial{i_n_subsections}{i_subsection}.n_bhvs);
    
        end
    end

    % Compute loss for current value of n_subsections.
    utility(i_n_subsections) = compute_utility(subsections_trial{i_n_subsections}, size(neural_data, 2), loss_parameters);
end

[~, i_max] = max(utility);
n_subsections_final = n_subsections_vals_to_try(i_max);

[subsections, other] = subdivide_session(neural_data, ...
                            labels, ...
                            other_label, ...
                            n_subsections_final); 

end



function utility = compute_utility(subsections, total_length, loss_params)

n_subsections = length(subsections);

durations = nan(n_subsections, 1);
purity = nan(n_subsections, 1);
diversity = nan(n_subsections, 1);

for i_subsection = 1:n_subsections
    durations(i_subsection) = length(subsections{i_subsection}.labels);
    purity(i_subsection) = subsections{i_subsection}.purity;
    diversity(i_subsection) = subsections{i_subsection}.n_bhvs;
end

durations = durations ./ total_length;
purity = median(purity);
diversity = median(diversity);

max_duration = max(durations);
total_duration = sum(durations);

utility = loss_params.max_duration * max_duration ...
          + loss_params.total_duration * total_duration ...
          + loss_params.purity * purity ...
          + loss_params.diversity * diversity;

end
