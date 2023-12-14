function [subsections, other] = choose_opt_subdivision(neural_data, labels, other_label, n_subsections_vals_to_try, loss_parameters)
% XXX: This function is incomplete. The idea is to create a utility
% function based on a weighted combination of richness (proportion of
% timepoints in a segment that are not "Other", i.e., proportion of
% timepoints that are "meaningfully" labeled) and diversity (proportion of
% all behaviors in session that occurred in the given segment) max duration
% (length of longest segment) and total duration (sum of lengths of all
% segments) for automatic segmentation. Very short segments can have very
% high richness and diversity scores, so a good segmentation should have
% some notion of length of segments (longer is better). The main thing left
% is to devise a suitable way to convert lengths varying across very wide
% scales (less than 10 timepoints up to thousands or tens of thousands) to
% a reliable/meaningful value on [0, 1], then validate on data etc.
% Currently the lengths enter directly into the utility function which will
% cause them to completely dominate the other terms. For example, one could
% set p = sum of lengths of all interactions in a session and define the
% max duration utility parameter as the length of the longest segment
% divided by p (similarly, total duration parameter is sum of segment
% lengths, divided by p). If each parameter's value is bounded on [0, 1],
% then the final utility should be bounded on [0, 1] if the weights sum to
% 1.
%
% Author: Jonathan Chien

if isempty(loss_parameters)
    loss_parameters.max_duration = 0.25;
    loss_parameters.total_duration = 0.25;
    loss_parameters.richness = 0.25;
    loss_parameters.diversity = 0.25;
end

error("Please read function documentation before proceeding.")

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
            fprintf('duration = %d. richness = %f. diversity = %d bhvs. \n', ...
                    length(subsections_trial{i_n_subsections}{i_subsection}.labels), ...
                    subsections_trial{i_n_subsections}{i_subsection}.richness, ...
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
richness = nan(n_subsections, 1);
diversity = nan(n_subsections, 1);

for i_subsection = 1:n_subsections
    durations(i_subsection) = length(subsections{i_subsection}.labels);
    richness(i_subsection) = subsections{i_subsection}.richness;
    diversity(i_subsection) = subsections{i_subsection}.n_bhvs;
end

durations = durations ./ total_length;
richness = median(richness);
diversity = median(diversity);

max_duration = max(durations);
total_duration = sum(durations);

utility = loss_params.max_duration * max_duration ...
          + loss_params.total_duration * total_duration ...
          + loss_params.richness * richness ...
          + loss_params.diversity * diversity;

end
