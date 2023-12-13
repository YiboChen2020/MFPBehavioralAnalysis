function [extracted, n_pts_rmv] ...
    = isolate_annotated_interactions(timeseries, session_raw, varargin)
% Pull out and concatenate sections of timeseries that have human
% annotations. Any binning is applied to each section individually first
% before final concatenation.
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
% Name-Value Pairs
%   'exclude'    : cell array where each cell contains string pattern. If
%                  any part of a behavioral annotation matches any of these
%                  patters, it will be excluded (i.e., the correponding
%                  portion of the timeseries will not be used).
%   'min_length' : (Scalar) minimum number of timepoints required for
%                  specified behavior; default = 1. If the session had
%                  fewer timepoints for this behavior than the minimum
%                  number, an exception will be thrown. Note that different
%                  bin and step sizes may affect the number of total
%                  timepoints retained (due to trimming).
%   'bin_width'  : (Scalar) number of timepoints in a bin.
%   'step_size'  : (Scalar) step size by which each bin will be advanced.
%                  Units are in timepoints.
%   'smooth'     : (5 (default) | scalar | 0), specify the window size (in
%                  original, unbinned timepoints) for gaussian smoothing,
%                  or specify a value between 0 and 1 (e.g., 0.2) which
%                  will for each annotated section set the smoothing window
%                  size to this fraction multiplied by the length of the
%                  given section, or set to 0/false to suppress smoothing.
%
% RETURNS
% -------
% extracted : m x p timeseries array, where the second array dim
%             consists of the concatenation of all m x p_i arrays, where i
%             is in {1, 2, ... s}, with s = number of epochs with human
%             annotations (that were not excluded), and m x p_i being the
%             i_th annotated timeseries array.
% n_pts_rmv : s x 2 (s = number of epoch with human annotation not
%             excluded) array where the 1st and 2nd column elements of the
%             i_th row are, respectively, the number of leading and
%             trailing timepoints removed during binning of that epoch.
%
% Author: Jonathan Chien. 10/3/22.


%% Parse inputs

p = inputParser;
addRequired(p, 'timeseries');
addRequired(p, 'session_raw');
addParameter(p, 'exclude', {'Other', 'Intro_', 'Rmv_'})
addParameter(p, 'min_length', 1);
addParameter(p, 'bin_width', 1);
addParameter(p, 'step_size', 1);
parse(p, timeseries, session_raw, behavior, varargin{:});
nv = p.Results;


%% Find timeperiods, bin, concatenate

% Get all_windows, an a x b matrix where a = the total number of epochs =
% length(session_raw.behaviors) and b = 2, where the 1st column is the
% start timepoint, and 2nd column is the stop timepoint.
all_windows = [session_raw.Fstart session_raw.Fstop];
n_annotations_tot = length(session_raw.behaviors);
assert(n_annotations_tot == size(all_windows, 1));

% For each of the epochs corresponding to an annotated interaction behavior, first bin, then
% concatenate to growing array. Keep track of number of timepoints retained
% from each epoch, and check at end to ensure that there was a sufficient
% number. Also keep track of the number of timepoints discarded for
% binning.
extracted = [];
n_leading_pts_rmv = [];
n_trailing_pts_rmv = [];

for i_annotation = 1:n_annotations_tot
    % Skip if current annotation is to be excluded.
    if ismember(session_raw.behvariors{i_annotation}, nv.exclude), continue; end

    % Else, get new annotation.
    new_annotation = timeseries(:, all_windows(i_annotation,1):all_windows(i_annotation,2));

    % Optionally smooth data.
    if (nv.smooth > 0) && (nv.smooth < 1)
        % Set window width for smoothing as supplied fraction of length of
        % current annotated section (this results in different timescales
        % of smoothing across the annotated sections).
        new_annotation = smoothdata(new_annotation, 2, 'gaussian', ...
                                    size(new_annotation, 2) * nv.smooth);
    elseif nv.smooth
        % Smooth using supplied window width.
        new_annotation = smoothdata(new_annotation, 2, 'gaussian', nv.smooth);
    end

    % Bin data. Add number of timepoints retained from current epoch to
    % tally across all epochs.
    [new_annotation_binned, n_timepts_used, n_leading, n_trailing] ...
        = bin_data(new_annotation, nv.bin_width, nv.step_size);

    % Add number of leading and trailing timepoints removed to struct.
    n_leading_pts_rmv = [n_leading_pts_rmv; n_leading];
    n_trailing_pts_rmv = [n_trailing_pts_rmv; n_trailing];

    % Check number of timepoints from current epoch against min length.
    if n_timepts_used < nv.min_length
       exception ...
           = MException("isolate_annotated_interactions:epoch_of_insufficient_length", ...
                        "The isolated epoch timeseries for %s (index into " + ...
                        "session_raw.behaviors = %d) for current " + ...
                        "subject and session is shorter than the " + ...
                        "specified minimum length of %d (see " + ...
                        "'min_length' name-value pair).", ...
                        session_raw{i_annotation}, i_annotation, nv.min_length);
       throwAsCaller(exception);
    end
    
    % Add binned data via concatenation.
    extracted = [extracted, new_annotation_binned];
end

% Horizontally concatenate.
n_pts_rmv = [n_leading_pts_rmv n_trailing_pts_rmv];

end
