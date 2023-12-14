function mean_activity ...
    = calc_mean_activity(timeseries, bhv_labels, bhvs_to_exclude, bhv_name_codes)
% Calculates mean activity of each behavior (across all pooled instances of
% that behavior). 
% 
% PARAMETERS
% ----------
% timeseries      : n_regions x n_timepoints array of neural activity.
% bhv_labels      : n_timepoints x 1 or 1 x n_timepoints array whose i_th
%                   element is the behavioral label of the i_th timepoint.
% bhvs_to_exclude : Cell array whose i_th element contains the name of one
%                   of the behaviors that we would like to exclude. Note
%                   that it is possible to include behaviors that are not
%                   actually present in a given session's data (they will
%                   simply be ignored). 
% bhv_name_codes  : Vector whose i_th element is the number code associated
%                   with the i_th element of bhvs_to_exclude. As with
%                   bhvs_to_exclude, extra codes will be ignored.
%
% RETURNS
% -------
% mean_activity : n_bhvs x n_regions matrix whose i_th, j_th element is the
%                 mean activity (across all epochs) of the j_th region
%                 during the i_th behavior.
%
% Author: Jonathan Chien. 7/19/22.


% Number of times each behavior occurs (unrelated to region).
unique_bhvs = unique(bhv_labels);
if isrow(bhv_labels), bhv_labels = bhv_labels'; end
if iscolumn(unique_bhvs), unique_bhvs = unique_bhvs'; end
n_occur_by_bhv = sum(bhv_labels == unique_bhvs);

% Accumarray function below will return a vector of length =
% max(unique_bhvs, [], 1). Thus, if there are any behavior codes that were
% skipped, the array will be too long. Here we need to identify behaviors
% that were skipped, so that their corresponding element in the vector
% output of accumarray can be removed below.
missing_ind = setdiff(1:max(unique_bhvs), unique_bhvs);

% Number of regions.
n_regions = size(timeseries, 1);

% For each region, calculate mean activity across behaviors. This results
% in an n_bhvs x n_regions matrix, which is what is returned.
mean_activity = NaN(length(n_occur_by_bhv), n_regions);
for i_region = 1:n_regions
    activity_sum_by_bhv = accumarray(bhv_labels, timeseries(i_region,:));
    activity_sum_by_bhv(missing_ind) = [];
    mean_activity(:,i_region) ...
        = activity_sum_by_bhv ./ n_occur_by_bhv';
end

% Optionally remove a behavior.
remove_ind = bhv_name_to_remove_ind(bhvs_to_exclude, ...
                                    bhv_name_codes, ...
                                    unique_bhvs);
mean_activity(remove_ind,:) = [];

end
