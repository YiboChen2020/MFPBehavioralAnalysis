function [binned_timeseries, n_timepts_to_use, n_leading_pts_rmv, n_trailing_pts_rmv] ...
    = bin_data(timeseries, bin_width, step_size)
% Accepts an n_vars x n_timepoints array, and bins according to specified
% bin width and step size.
% 
% PARAMETERS
% ----------
% timeseries : n_vars x n_timepoints array of data to be binned.
% bin_width : Desired bin width (as number of timepoints). Set both 
%             bin_width and step_size to 1 to eschew binning.
% step_size : Desired step size (as number of timepoints). Set both 
%             bin_width and step_size to 1 to eschew binning. Set step_size
%             to be less than bin_width to achieve overlap of bins.
% 
% RETURNS
% -------
% binned_timeseries  : n_vars x n_bins array of data formed by averaging
%                      input timeseries within bins after trimming.
% n_timepts_to_use   : Scalar number of timepoints from the original
%                      timeseries retained for binning (i.e., max number of
%                      timepoints determined to be compatible for desired
%                      bin and step size.
% n_leading_pts_rmv  : Scalar number of timepoints trimmed from head of
%                      timeseries to accomodate bin width and step size.
% n_trailing_pts_rmv : Scalar number of timepoints trimmed from tail of
%                      timeseries to accomodate bin width and step size.
%
% Author: Jonathan Chien. Refactored from local fx within xcorr_by_bhv on
% 5/31/22. Last edit 10/3/22.

% Determine number of timepoints/bins compatible with desired bin width and
% step size.
n_timepts_tot = size(timeseries, 2);
[n_timepts_to_use, n_bins] ...
    = determine_n_bins(n_timepts_tot, bin_width, step_size);

% Trim ends of timeseries in balanced fashion to preserve middle (which
% is likely to be more accurate).
trimmed_timeseries = trim_ends(timeseries, n_timepts_to_use);

% Get number of leading and trailing timepoints removed. If the number of
% points removed is odd, the extra point will be removed from the end.
n_leading_pts_rmv = floor((n_timepts_tot - n_timepts_to_use) / 2);
n_trailing_pts_rmv = ceil((n_timepts_tot - n_timepts_to_use) / 2);

% Prepare indices of timepoints at beginning of each bin.
bin_starts = 1 : step_size : (n_bins-1) * step_size + 1;

% Bin data by averaging.
binned_timeseries = NaN(size(trimmed_timeseries, 1), n_bins);
for i_bin = 1:n_bins
    binned_timeseries(:,i_bin) = mean( ...
        trimmed_timeseries(:, bin_starts(i_bin):bin_starts(i_bin) + bin_width -1), 2 ...
                                      );        
end

end


% --------------------------------------------------
function trimmed_timeseries = trim_ends(timeseries, n_timepts_to_use)
% For a timeseries of length m, trim ends of timeseries so that it is of
% length n. Trimming is performed equally from both ends if m - n is even,
% and as evenly as possible if m - n is odd, with the extra timepoint
% trimmed from the end of the timeseries. 

n_timepts_to_remove = size(timeseries, 2) - n_timepts_to_use;

% Copy input.
trimmed_timeseries = timeseries;

% Handle parity of number of points to be removed with floor/ceil.
trimmed_timeseries(:, 1:floor(n_timepts_to_remove/2)) = [];
trimmed_timeseries(:, end - (ceil(n_timepts_to_remove/2) - 1):end) = [];

assert(size(trimmed_timeseries, 2) == n_timepts_to_use);

end
