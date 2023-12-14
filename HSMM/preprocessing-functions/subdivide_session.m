function [subsections, other] = subdivide_session(timeseries, labels, other_label, n_subsections)
%
% Find Other periods that are particularly long. This may indicate a long
% period of separation betweens subjects possibly acting as a natural
% boundary between interactions (rather than a brief separation in the
% middle of an interaction). These Other periods are identified and removed
% to segment the session into the desired number of shorter subsections.
% 
% PARAMETERS
% ----------
% timeseries    : m x n array (e.g. timeseries).
% labels        : n-vector of labels whose i_th element corresponds to the 
%                 i_th column of "timeseries".
% other_label   : For A = unique(labels), other_label is a single element
%                 of A (it is the label for all "Other" timepoints).
% n_subsections : Positive integer specifying the number of desired
%                 subsections into which "timeseries" will be segmented.
%
% RETURNS
% -------
% subsections : n_subsections x 1 cell array, where the i_th cell contains 
%               a struct where the i_th cell contains a struct whose fields
%               have a segment of the overall timeseries, along with the
%               corresponding labels for the timepoints in that segement,
%               and the indices of the first and last timepoints of the
%               segment with respect to timeseries. 
% other       : k x n matrix, where k is the number of occurences of an
%               "Other" epoch (a consecutive sequence of Other labels). The
%               i_th row corresponds to the i_th such epoch. For the i_th
%               row, the 1st column element is the index (column index of
%               timeseries) of the first timepoint in this "Other" epoch,
%               the 2nd element is the index of the last timepoint in this
%               epoch, and the 3rd column element is the duration of the
%               epoch.
%
% Author: Jonathan Chien


% Check input.
assert(n_subsections > 0 && floor(n_subsections) == n_subsections, ...
       "n_subsections must be a positive integer.")

n_timepoints = length(labels);

% Will add rows to other. Each row corresponds to another Other period.
% First row element is head index of that period, second element is tail
% index, third element is the duration of that period.
other = [];
i_timepoint = 1;
while i_timepoint <= n_timepoints
    if labels(i_timepoint) == other_label
        % Record current index as head. 
        head = i_timepoint;

        % Set in_other flag to true to begin while loop.
        in_other = true;

        while in_other && i_timepoint 
            % Increment timepoint and check whether still Other.
            i_timepoint = i_timepoint + 1;

            % If no longer Other (or reached end of timeseries), record
            % tail index and duration.
            if i_timepoint > n_timepoints || labels(i_timepoint) ~= other_label
                tail = i_timepoint - 1;
                duration = tail - head + 1;
                in_other = false;
            % If still Other, increment time index again.       
            elseif labels(i_timepoint) == other_label
                i_timepoint = i_timepoint + 1;
            end      
        end

        % Add head index, tail index, and duration from the Other period
        % that was just completed.
        other = [other; [head tail duration]];
    end

    i_timepoint = i_timepoint + 1;
end


%% Subdivide by longest Other period

% Get the n_subsections largest Other periods (what we are getting are
% their indices among all of the Other periods).
[~, other_ind_all] = sort(other(:,3), 'descend');
if length(other_ind_all) >= n_subsections
    other_ind = other_ind_all(1:n_subsections);
else
    exception = MException('subdivide_session:n_subsections_gr_than_n_other', ...
                           ['The number of subsections requested (%d) is ' ...
                            'greather than the number of Other periods ' ...
                            '(%d).'], n_subsections, length(other_ind_all));
    throw(exception);
end


% Sort these Other indices by their order of occurence (using head
% indices), since they were originally ordered by length.
[~, occurence_ind] = sort(other(other_ind, 1));
other_ind_by_occurence = other_ind(occurence_ind);

subsections = cell(n_subsections, 1);
for i_subsection = 1:n_subsections
    % Handle head index of first subsection, since head will be start of
    % timeseries.
    if i_subsection == 1
        head = 1;
    else
        % Use the first point after the tail of the last Other period as
        % the head index.
        head = other(other_ind_by_occurence(i_subsection), 2) + 1;
    end

    % Handle tail index of last subsection, since this will be the end of
    % the timeseries. 
    if i_subsection == n_subsections
        tail = size(timeseries, 2);
    else
        % If not the last subsection, this current period will run up to
        % right before the head index of the next Other period.       
        tail = other(other_ind_by_occurence(i_subsection + 1), 1) - 1;       
    end

    % Get delineated portion of data and labels.
    subsections{i_subsection}.timeseries = timeseries(:,head:tail);
    subsections{i_subsection}.labels = labels(head:tail);
    subsections{i_subsection}.head = head;
    subsections{i_subsection}.tail = tail;

    % Compute how much of the subsection is not Other.
    subsections{i_subsection}.richness ...
        = sum(subsections{i_subsection}.labels ~= other_label) / (tail - head + 1);

    % Compute how many behaviors occured (not including Other). 
    subsections{i_subsection}.n_bhvs ...
        = length(unique(subsections{i_subsection}.labels(subsections{i_subsection}.labels ~= other_label)));
end

end
