function mean_durations ...
    = calc_mean_duration(bhv_labels, bhvs_to_exclude, bhv_name_codes)
% Calculate the mean duration (in units of timepoints), i.e., the mean
% epoch length, for each behavior. bhv_labels corresponds to the bhv_labels
% output of labeled_timeseries and is an n x 1 vector where n is the number
% of timepoints.
%
% PARAMETERS
% ----------
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
% mean_durations : n_bhvs x 1 vector whose i_th element is the mean
%                  duration (in units of timepoints) of the i_th behavior.
%
% Author: Jonathan Chien. 7/19/22.

unique_bhvs = unique(bhv_labels);
n_bhvs = length(unique_bhvs);

% n_timepoints x n_bhvs array of boolean indices.
if isrow(bhv_labels), bhv_labels = bhv_labels'; end
if iscolumn(unique_bhvs), unique_bhvs = unique_bhvs'; end
bhv_ind = bhv_labels == unique_bhvs;

% Prepend zeros so that the behavior that occurs first can be detected (it
% will have a 1 in the first timepoint rather than all behaviors having
% zero at the first timepoint). Also notice that the -1 ending marker
% occurs at the timepoint after the last annotated timepoint for a
% behavior. However, since we would have otherwise had to have added + 1
% (since an epoch from timepoints 118 to 120 really consists of three
% timepoints in the epoch, even though the difference is 2), we can just
% directly use the timepoint marked by -1. 
markers = diff([zeros(1, n_bhvs); bhv_ind], 1, 1);

mean_durations = NaN(n_bhvs, 1);
for i_bhv = 1:n_bhvs
    % Get start and endpoints of each epoch.
    startpoints = find(markers(:,i_bhv) == 1);
    endpoints = find(markers(:,i_bhv) == -1);

    % Handle edge case where behavior is last in session (thus no -1
    % flagging end of final epoch).
    if length(startpoints) > length(endpoints)
        assert(bhv_ind(end,i_bhv), ...
               "It appears that there are more startpoints than endpoints " + ...
               "for behavior code %d, but this behavior was not the last in " + ...
               "session either.", unique_bhvs(i_bhv))
        
        % The plus 1 is for the same reason as using the -1 as explained in
        % the above verbose comment.
        endpoints = [endpoints; size(bhv_ind, 1) + 1];
    end

    % Calculate mean duration acrosss epochs.
    mean_durations(i_bhv) = mean(endpoints - startpoints);   
end

% Optionally remove a behavior.
remove_ind = bhv_name_to_remove_ind(bhvs_to_exclude, ...
                                    bhv_name_codes, ...
                                    unique_bhvs);
mean_durations(remove_ind) = [];

end
