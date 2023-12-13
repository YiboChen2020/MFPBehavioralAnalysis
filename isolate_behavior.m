function activity = isolate_behavior(timeseries, session_raw, behavior, varargin)
% Pull out sections of timeseries, for each region, that correspond to the
% specified behavior of interest. 
% 
% PARAMETERS
% ----------
% timeseries  : m x n timeseries array, where m = number of variables, and
%               n = number of timepoints/bins.
% session_raw : Scalar struct output (possibly within cell array) from
%               extract_single_subj. Needed here are the fields .behavior
%               (list of all behaviors annotated within session) and
%               .Fstart and .Fstop (containing respectively, the start and
%               stop timepoints of the annotated behavior).
% behavior    : Behavior of interest which we would like to pull out.
%               String value, watch capitalization.
% Name-Value Pairs
%   'min_length' : (Scalar) minimum number of timepoints required for
%                  specified behavior; default = 1. If the session had
%                  fewer timepoints for this behavior than the minimum
%                  number, an exception will be thrown. Note that different
%                  bin and step sizes may affect the number of total
%                  timepoints retained (due to trimming).
%   'bin_width'  : (Scalar) number of timepoints in a bin.
%   'step_size'  : (Scalar) step size by which each bin will be advanced.
%                  Units are in timepoints.
%
% RETURNS
% -------
% activity : m x p timeseries array, where the second array dim consists
%            of the concatenation of all m x p_i arrays, where i is in {1,
%            2, ... s}, with s = number of epochs with behavior of
%            interest, and m x p_i being the i_th timeseries array
%            annotated as behavior of interest.
%
% Author: Jonathan Chien. Refactored from local fx within xcorr_by_bhv on
% 5/31/22. Switched from bin_data to bin_nd_data on 12/22/22. 
% Added option for Investigate_M vs Investigate_F on 1/3/23.


%% Parse inputs

p = inputParser;
addRequired(p, 'timeseries');
addRequired(p, 'session_raw');
addRequired(p, 'behavior');
addParameter(p, 'min_length', 1);
addParameter(p, 'deriv', false);
addParameter(p, 'bin_width', 1);
addParameter(p, 'step_size', 1);
parse(p, timeseries, session_raw, behavior, varargin{:});
nv = p.Results;


%% Find timeperiods, bin, concatenate

% Handle Investigate_M and Investigate_F.
if ismember(behavior, {'Investigate_M', 'Investigate_F'})
    if strcmp(behavior, 'Investigate_M'), flag = 'M'; end
    if strcmp(behavior, 'Investigate_F'), flag = 'F'; end
    behavior = 'Investigate';
else
    flag = nan;
end

% Ensure current session has requested behavior. If not throw exception
% with identifier in invoking function; it seems throwAsCaller will
% continue to pass the MException object up through calling functions until
% it reaches a try/catch block (or perhaps the highest level invoking
% function?).
if ~ismember(behavior, session_raw.behaviors)
    exception = MException("isolate_behavior:behavior_not_found", ...
                           "The specified behavior did not occur " + ...
                           "during the passed-in session.");
    throwAsCaller(exception);
end

% If behavior was one of Investigate_M or Investigate_F, first find valid
% epochs.
if ~isnan(flag)
    intro_ind = find(ismember(session_raw.behaviors, sprintf('Intro_%s', flag)));
    remove_ind = find(ismember(session_raw.behaviors, sprintf('Rmv_%s', flag)));
    if length(intro_ind) ~= length(remove_ind)
        shorter = min(length(intro_ind), length(remove_ind));
        intro_ind = intro_ind(1:shorter);
        remove_ind = remove_ind(1:shorter);
    end

    valid_ind = [];
    for i_epoch = 1:length(intro_ind)
        valid_ind = [valid_ind; (intro_ind(i_epoch):remove_ind(i_epoch))'];
    end
else
    valid_ind = 1:length(session_raw.behaviors);
end

% Get m x n matrix woi (windows of interest) where m = the number of epochs
% corresponding to the behavior of interest (the total number of epochs =
% length(session_raw.behaviors) and n = 2, where the 1st column is the
% start timepoint, and 2nd column is the stop timepoint.
all_windows = [session_raw.Fstart session_raw.Fstop];
woi = all_windows(intersect(find(ismember(session_raw.behaviors, behavior)), valid_ind),:);

% For each of the epochs corresponding to behavior, first bin, then
% concatenate. Keep track of number of timepoints retained from each epoch,
% and check at end to ensure that there was a sufficient number.
activity = [];
n_timepts_total = 0;
for i_epoch = 1:size(woi, 1)
    % Get new epoch.
    new_period = timeseries(:, woi(i_epoch,1):woi(i_epoch,2));

    % Optionally take first derivative.
    if nv.deriv, new_period = [zeros(size(new_period, 1), 1) diff(new_period, 1, 2)]; end

    % Bin data. Add number of timepoints retained from current epoch to
    % tally across all epochs.
    [new_period_binned, n_timepts] ...
        = bin_nd_data(new_period, 2, nv.bin_width, nv.step_size);
    n_timepts_total = n_timepts_total + n_timepts;
    
    % Add binned data via concatenation.
    activity = [activity, new_period_binned];
end

% Ensure that number of timepoints retained across all epochs meets the
% minimum criterion.
if n_timepts_total < nv.min_length
   exception = MException("isolate_behavior:timeseries_of_insufficient_length", ...
                          "The isolated timeseries for %s for current " + ...
                          "subject and session is shorter than the " + ...
                          "specified minimum length of %d (see " + ...
                          "'min_length' name-value pair).", ...
                          behavior, nv.min_length);
   throwAsCaller(exception);
end

% TODO: add resampling here.

end
