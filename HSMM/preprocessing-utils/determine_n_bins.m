function [n_timepts_to_use, n_bins] ...
    = determine_n_bins(n_timepts, bin_width, step_size)
% If we wish to bin the data, the length of the timeseries passed in might
% not be compatible with the bin/step sizes. This function takes in the
% number of points in the timeseries, as well as the desired bin width and
% step size and decreases the number of timepoints if necessary to fit with
% the specified bin width and step size.
% 
% PARAMETERS
% ----------
% n_timepts : Number of timepoints in full original timeseries.
% bin_width : Desired bin width (as number of timepoints). Set both 
%             bin_width and step_size to 1 to eschew binning.
% step_size : Desired step size (as number of timepoints). Set both 
%             bin_width and step_size to 1 to eschew binning. Set step_size
%             to be less than bin_width to achieve overlap of bins.
% 
% RETURNS
% -------
% n_timepts_to_use : Number of timepoints in trimmed timeseries that can
%                    be binned.
% n_bins           : Number of bins that can be fitted to trimmed
%                    timeseries.
%
% Author: Jonathan Chien. Refactored from local fx within xcorr_by_bhv on
% 5/31/22.


% alpha is the numnber of timepoints minus the bin width, divided by the
% step size. If alpha is an integer, it is the number of bins. 
calc_alpha = @(x, y, z) ((x - y) / z) + 1;
alpha = calc_alpha(n_timepts, bin_width, step_size);

if floor(alpha) == alpha
    % alpha is an integer and is thus the number of bins.
    n_bins = alpha;
    n_timepts_to_use = n_timepts;
else
    % Set n_timepoints as initial number of timepoints and decrease until
    % alpha becomes an integer.
    n_timepts_to_use = n_timepts;
    counter = 0;
    while floor(alpha) ~= alpha
        n_timepts_to_use = n_timepts_to_use - 1;
        alpha = calc_alpha(n_timepts_to_use, bin_width, step_size);

        % Graceful exit if timeseries was too short/incompatible with
        % bin/step size.
        counter = counter + 1;
        if counter == n_timepts
            exception = MException('determine_n_bins:no_timepoints_left', ...
                                   "All timepoints were removed before " + ...
                                   "an acceptable value for n_bins was " + ...
                                   "determined.");
            throwAsCaller(exception);
        end
    end

    n_bins = alpha;
end

end
