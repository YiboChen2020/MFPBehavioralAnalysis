function activity = isolate_baseline_2(timeseries, behaviors, Fstart, varargin)
% Pull out sections of Lfold, for each region, that correspond to baseline.
% This is done in a separate function from the other behaviors because the
% baseline is best defined as the timeperiod between two labels
% ("Intro_Baseline" and "Rmv_Baseline"). This is the same as the original
% isolate_baseline function, but the argument session_raw has been replaced
% by the behaviors and Fstart arguments, since these fields of session_raw
% were all that were needed.
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
% activity : m x p timeseries array, where p = number of timepoints in the
%            baseline period.
%
% Author: Jonathan Chien


%% Parse inputs

p = inputParser;
addRequired(p, 'timeseries');
addRequired(p, 'behaviors');
addRequired(p, 'Fstart');
addParameter(p, 'min_length', 1);
addParameter(p, 'bin_width', 1);
addParameter(p, 'step_size', 1);
parse(p, timeseries, behaviors, Fstart, varargin{:});
nv = p.Results;


%% Find timeperiods

% If session contains both an 'Intro_Baseline' and 'Rmv_Baseline' label,
% use the intervening timeperiod as the baseline.
if ismember('Intro_Baseline', behaviors) ...
   && ismember('Rmv_Baseline', behaviors)
    
    % Get index of first and last timepoints in baseline.
    start_point = Fstart(find(ismember(behaviors, 'Intro_Baseline'), 1));
    end_point = Fstart(find(ismember(behaviors, 'Rmv_Baseline'), 1));
                                   
    activity = timeseries(:,start_point:end_point);

% If session lacks either 'Intro_Baseline', 'Rmv_Baseline', or both, search
% for "Intro" and use period before as the baseline. This assumes that the
% "Intro" is part of a label such as "Intro_Toy" or "Intro_M", etc., i.e.,
% the introduction of a stimulus. Warn the user if this route is taken. If
% this fails, throw exception passed to try/catch in an invoking function
% for further exception handling.
elseif any(~cellfun(@isempty, strfind(behaviors, 'Intro')))
    i_bin_with_intro = behaviors( ...
        find(~cellfun(@isempty, strfind(behaviors, 'Intro')), 1) ...
                                             );
    warning("'Intro_Baseline' and 'Rmv_Baseline' labels were not detected " ...
            + "in the current session. All activity prior to the occurence " ...
            + "of the label %s will be used as the baseline period", ...
            i_bin_with_intro);
    activity = timeseries(:,1:i_bin_with_intro);
    
else
    % As of 6/1/22, this exception is caught in a try/catch block in
    % snr_all_subjects_sessions.
    exception = MException("isolate_baseline:unable_to_isolate_baseline", ...
                           "Failed to identify a baseline period.");
    throwAsCaller(exception);  
end


%% Check against minimum length and bin

% Throw exception as caller (for exception handling at wrapper level across
% subjects/sessions if minimum timeseries length not met.
if size(activity, 2) < nv.min_length
    exception = MException("isolate_baseline:timeseries_of_insufficient_length", ...
                           "The isolated baseline timeseries for current " + ...
                           "subject and session is shorter than the " + ...
                           "specified minimum length of %d (see " + ...
                           "'min_length' name-value pair).", ...
                           nv.min_length);
    throwAsCaller(exception);
end

% Bin data.
[activity, ~, ~, ~] = bin_nd_data(activity, 2, nv.bin_width, nv.step_size, ...
                                  'trim_from', 'back', 'bin_excess', 'excess_only');

end
