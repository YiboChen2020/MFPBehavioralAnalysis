function rmv_ind = remove_other_ind(labels, other_label, n_buffer)
% For a vector of labels for a timeseries, identify all timepoints
% corresponding to extended sequence of "Other" labels and mark these for
% removal, leaving a specified "buffer" of "Other" timepoints.
%
% PARAMETERS
% ----------
% labels      : k-vector whose i_th element is the behavioral label for the
%               i_th timepoint.
% other_label : Scalar, one of the behavioral labels, specifies which label
%               corresponds to "Other".
% n_buffer    : Non-negative scalar, specifies how many "Other" timepoints
%               from the ends of the "Other" epoch to leave. E.g., if
%               n_buffer = 5, and a given "Other" epoch is 70 timepoints
%               long, the first 2 and last 3 timepoints of the "Other"
%               epoch will be left, and the indices of the middle 65
%               timepoints (from the 3rd to 67th timepoint) of the "Other"
%               epoch will be returned.
%
% RETURNS
% -------
% rmv_ind : Vector of indices marking "Other" timepoints for removal.
%
% Author: Jonathan Chien

% Get logical indices of all timepoints marked Other.
other_ind = labels == other_label;

% Total number of timpoints.
n_timepoints = length(other_ind);

% Assuming that there are n periods of time (epochs) marked "Other", we
% want to get the indices of each of these timepoints, but stored in a
% separate cell (so that we can then trim from the ends of these indices
% vectors to leave some Other timepoints, in order to preserve transition
% statistics).
other_epochs = {};
i_timepoint = 1;
in_series = false;
curr_series = [];
while i_timepoint <= n_timepoints
    % If current timepoint is marked Other, and we are not already in the
    % middle of a series of Other timepoints, begin a new series (a new
    % Other epoch).
    if ~in_series && other_ind(i_timepoint) == 1 
        in_series = true;
        curr_series = i_timepoint;

        % This handles the (very unlikely) event that the last timepoint is
        % Other but the second to last timepoint is not.
        if i_timepoint == n_timepoints
            in_series = false;
            other_epochs = cat(1, other_epochs, {curr_series});
            curr_series = [];
        end

    % If currently in series, check if current timepoint is other, and add
    % if so.
    elseif in_series && other_ind(i_timepoint) == 1 
        curr_series = [curr_series; i_timepoint];

        % This handles the (quite likely) event that the last timepoint is
        % part of an Other epoch (that doesn't consist just of that last
        % timepoint).
        if i_timepoint == n_timepoints
            in_series = false;
            other_epochs = cat(1, other_epochs, {curr_series});
            curr_series = [];
        end

    % If current timepoint is not in series, set in_series to false and add
    % series to cell array other_epochs.
    elseif in_series && other_ind(i_timepoint) == 0 
        in_series = false;
        other_epochs = cat(1, other_epochs, {curr_series});
        curr_series = [];
    end

    i_timepoint = i_timepoint + 1;
end

% Remove end points of remove_ind for each Other epoch so that there are
% n_buffer Other timepoints left.
n_other_epochs = length(other_epochs);
rmv_ind = [];
half_buffer = n_buffer/2;
for i_epoch = 1:n_other_epochs
    curr_epoch = other_epochs{i_epoch};
    n_timepoints = length(curr_epoch);

    if n_timepoints < n_buffer
        exception ...
            = MException('remove_other_ind:n_buffer_too_long', ...
                         ['n_buffer = %d requested, but Other epoch number' ...
                          ' %d was only %d timepoints long.'], n_buffer, i_epoch, n_timepoints);
        throwAsCaller(exception);
    else
        % If n_timepoints == n_buffer, this is returned empty, else,
        % indices of middle n_timepoints - n_buffer timepoints (within
        % current Other epoch) are returned.
        curr_rmv_ind = setdiff(1:n_timepoints, ...
                       [1:floor(half_buffer), (n_timepoints - ceil(half_buffer) + 1):n_timepoints]);
    end
           
    rmv_ind = [rmv_ind; curr_epoch(curr_rmv_ind)];
end

end
