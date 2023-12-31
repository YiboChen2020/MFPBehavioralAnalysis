function [X_binned, n_timepts_to_use, n_excess_leading_pts, n_excess_trailing_pts] ...
    = bin_nd_data(X, dim, bin_width, step_size, nv)
% Accepts an n-dim array, where the d_th dimension (given by dim) is e.g.
% time, and bins along the d_th dimension according to the specified bin
% width and step size. Both non-overlapping and overlapping bins are
% supported (see bin_width and step_size under PARAMETERS).
% 
% PARAMETERS
% ----------
% X         : n-d array of data to be binned along the d_th dimension
%             (specified by dim).
% dim       : Positive integer 1 <= dim <= n denoting the dimension of X
%             along which we wish to bin (e.g., the dimension that is
%             time).
% bin_width : Desired bin width (as number of timepoints). Set both 
%             bin_width and step_size to 1 to eschew binning (e.g., if this
%             is function is called in a larger pipeline). 
% step_size : Desired step size (as number of timepoints). Set both 
%             bin_width and step_size to 1 to eschew binning (e.g., if this
%             is function is called in a larger pipeline). Set step_size
%             to be less than bin_width to achieve overlap of bins.
% Name-Value Pairs (nv)
%   'excess_at'  : ('back' | 'front' | 'both_ends'). If the number of
%                  timepoints is not evenly divisble by the bin size, the
%                  remainder/excess needs to be handled. If 'excess_at' =
%                  'back', the leading edge of the first bin is aligned
%                  with the start of the timeseries such that the excess
%                  timepoints are at the end or "back" of the timeseries.
%                  If 'excess_at' = 'front', the trailing edge of the last
%                  bin is aligned with the end of the timeseries, and the
%                  excess timepoints are at the start of the timeseries. If
%                  'excess_at' = 'both_ends', for n excess timepoints,
%                  floor(n/2) and ceil(n/2) timepoints are removed from the
%                  front and back of the timeseries, respectively. If
%                  'bin_excess' (below) is false, the excess timepoints
%                  (under any of the three schemes described above) are
%                  dropped, changing the original timeseries length. If the
%                  points are to be kept, the exact method is specified by
%                  'bin_excess'.
%   'bin_excess' : ('excess_only' | 'same_bin_size' | false). If
%                  false, any excess timepoints are removed. If
%                  'excess_only', the excess timepoints are averaged to
%                  create an additional bin, appended to the end of the
%                  timeseries if 'excess_at' = 'back', to the start if
%                  'excess_at' = 'front', and to both ends if 'excess_at' =
%                  'both_ends'. If 'bin_excess' = 'same_bin_size', an
%                  additional bin is added as with 'excess_only', but the
%                  additional bin is of the same size as all others, thus
%                  resulting in overlap. For example, if:
%
%                      X is an n-D array featuring 14 timepoints, 
%                      bin_width = step_size = 5,
%                      'excess_at' = 'back', 
%
%                      then the two main bins will have {start_index,
%                      stop_index} as {1, 5} and {6, 10}. Timepoints 11-14
%                      are excess (and at the "back").
%                         
%                      If 'bin_excess' = 'excess_only', a third bin
%                      consisting of timepoints from {11, 14} will be
%                      appended to the back. If 'bin_excess' =
%                      'same_bin_size', the third bin will consist of
%                      timepoints from {10, 14}, thus preserving the bin
%                      size (5) but introducing overlap between the last
%                      two bins.
% 
% RETURNS
% -------
% X_binned           : n-d array of data formed by averaging
%                      input along the d_th dimension. All dim sizes are
%                      the same, except for the d_th dimension size, which
%                      now equals the number of bins.
% n_timepts_to_use   : Scalar number of timepoints from the original
%                      timeseries retained for binning. If 'bin_excess' =
%                      false, this will be the max number of timepoints
%                      determined to be compatible for desired bin and step
%                      size, as the excess timepoints were trimmed off. If
%                      'bin_excess' is not false, n_timepts_to_use equals
%                      the original number of timepoints.
% n_leading_pts_rmv  : Scalar number of timepoints deemed as excess at the
%                      start of timeseries to accomodate bin width and step
%                      size (see 'excess_at' and 'bin_excess' Name-Value
%                      pairs above).
% n_trailing_pts_rmv : Scalar number of timepoints deemed as excess at the
%                      end of timeseries to accomodate bin width and step
%                      size (see 'excess_at' and 'bin_excess' Name-Value
%                      pairs above).
%
% Author: Jonathan Chien. 10/3/22.

arguments
    X
    dim
    bin_width
    step_size
    nv.excess_at = 'back' % 'back', 'front', or 'both_ends'
    nv.bin_excess = 'excess_only' % 'excess_only', 'same_bin_size', or false
end

% Check input data.
% assert(~any(isnan(X), 'all'), "Input data should not contain NaN values.")

% Determine number of timepoints/bins compatible with desired bin width and
% step size.
n_timepts_tot = size(X, dim);
[n_timepts_to_use, n_bins] ...
    = determine_n_bins(n_timepts_tot, bin_width, step_size);

% Place temporal dimension at front and flatten remaining dims.
[X_perm_2d, X_perm_dim_sizes] = to_2d_array(X, dim);

% Trim ends of timeseries as specified.
[X_perm_2d_trimmed, n_excess_leading_pts, n_excess_trailing_pts] ...
    = trim_ends(X_perm_2d, n_timepts_to_use, nv.excess_at);

% Prepare indices of timepoints at beginning of each bin.
bin_starts = 1 : step_size : (n_bins-1) * step_size + 1;

% Bin data by averaging. This is just the part that is evenly divisible by
% the number of requested bins.
X_perm_2d_trimmed_binned = NaN(n_bins, size(X_perm_2d_trimmed, 2));
for i_bin = 1:n_bins
    X_perm_2d_trimmed_binned(i_bin,:) = mean( ...
        X_perm_2d_trimmed(bin_starts(i_bin):bin_starts(i_bin) + bin_width - 1, :), ...
        1 ...
                                             );        
end

% Optionally create leading and trailing bins as necessary to handle extra
% timepoints.
if nv.bin_excess
    % If the leading and trailing timepoints are to be included, we can
    % either put them in their own bin (only excess), or we can create bins
    % of the same size as the others, but with the first timepoint of front
    % bin as the first of the leading timepoints and the last timepoint of
    % the back bin as the last of the trailing timepoints.
    if strcmp(nv.bin_excess, 'same_bin_size')
        % Check that there were excess points trimmed off the front. If so,
        % create a bin of the same size, beginning with the first of these
        % points. If no excess points at the front, set to 0.
        if n_excess_leading_pts ~= 0
            front_bin_size = bin_width; 
        else
            front_bin_size = 0;
        end

        % Check that there were excess points trimmed off the back. If so,
        % create a bin of the same size, where the last timpoint in the bin
        % is the last of the trimmed off timepoints. If no excess points at
        % the back, set to 0.
        if n_excess_trailing_pts ~= 0
            back_bin_size = bin_width;
        else
            back_bin_size = 0;
        end
    elseif strcmp(nv.bin_excess, 'excess_only')
        front_bin_size = n_excess_leading_pts;
        back_bin_size = n_excess_trailing_pts;
    else
        error("Unrecognized value for 'bin_excess'. Must be " + ...
              "'same_bin_size' | 'excess_only'.")
    end

    % Bin excess points as determined above.
    front_bin = X_perm_2d(1:front_bin_size, :);
    if ~isempty(front_bin), front_bin = mean(front_bin, 1); end % Do this on separate line or will generate NaNs
    back_bin = X_perm_2d((end-back_bin_size+1):end, :);
    if ~isempty(back_bin), back_bin = mean(back_bin, 1); end % Do this on separate line or will generate NaNs
    
    % Add any new bins of excess points.
    X_perm_2d_trimmed_binned = [front_bin; X_perm_2d_trimmed_binned; back_bin];
    
    % Update number of bins, dim sizes, timepoints used (all are now used).
    n_bins = n_bins + ~isempty(front_bin) + ~isempty(back_bin);
    X_perm_dim_sizes(1) = n_bins;
    n_timepts_to_use = n_timepts_tot;
end

% Reshape permuted array back to original array shape/dim order.
X_binned = to_nd_array(X_perm_2d_trimmed_binned, dim, X_perm_dim_sizes);

end


%% Local functions

% --------------------------------------------------
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
% Author: Jonathan Chien. Refactored from local fx within xcorr_by_bhv to
% separate function on 5/31/22. Added as local function here on 12/21/22.


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


% --------------------------------------------------
function [X_perm_2d, perm_array_sizes] = to_2d_array(X, d)
% X is an n-D array. d is a positive integer 1 <= d <= n, denoting the
% dimension along which we will average (bin). The d_th dimension will be
% moved to the 1st dimension. The permuted array will then be flattened to
% yield a 2D array where the i_th row of the flattened array is the i_th
% slice along the d_th dimension of the original array.

% Get shape and dimensionality of array.
array_sizes = size(X);
n = length(array_sizes); 

% Check for valid d.
assert((1 <= d) && (d <= n), "The following must hold: 1 <= d <= n.");

% Place temporal dimension first.
perm_dim_order = [d setdiff(1:n, d)];
X_perm = permute(X, perm_dim_order); 
perm_array_sizes = array_sizes(perm_dim_order);

% Flatten to 2D array.
X_perm_2d = reshape(X_perm, size(X_perm, 1), []);

end


% --------------------------------------------------
function X = to_nd_array(X_perm_2d_trimmed_binned, d, perm_array_sizes)
% perm_array_sizes here are the dim sizes of the permuted n-D array (with
% the d_th dimension first).

% Get dimensionality of array.
n = length(perm_array_sizes);

% Check for valid d.
assert((1 <= d) && (d <= n), "The following must hold: 1 <= d <= n.");

% Reshape into n-D array. We use the first dim size of X_perm_2d, since it
% may have been trimmed and binned. The other dim sizes should not change.
X_perm = reshape(X_perm_2d_trimmed_binned, ...
                 [size(X_perm_2d_trimmed_binned, 1) perm_array_sizes(2:end)]);

% Place temporal dimension back to its original position.
X = permute(X_perm, [2:d 1 d+1:n]);

end


% --------------------------------------------------
function [X_trimmed, n_leading_timepts_removed, n_trailing_timepts_removed] ...
    = trim_ends(X, n_timepts_to_use, excess_at)
% For a timeseries of length m (m x p array X), trim ends of timeseries so
% that it is of length n. Trimming is performed equally from both ends if m
% - n is even, and as evenly as possible if m - n is odd, with the extra
% timepoint trimmed from the end of the timeseries.

n_timepts_to_remove = size(X, 1) - n_timepts_to_use;

% Copy input.
X_trimmed = X;

% Trim excess timepoints.
if strcmp(excess_at, 'front')
    X_trimmed(1:n_timepts_to_remove,:) = [];
    n_leading_timepts_removed = n_timepts_to_remove;
    n_trailing_timepts_removed = 0;
elseif strcmp(excess_at, 'back')
    X_trimmed((end-n_timepts_to_remove+1):end,:) = [];
    n_leading_timepts_removed = 0;
    n_trailing_timepts_removed = n_timepts_to_remove;
elseif strcmp(excess_at, 'both_ends')
    % Handle parity of number of points to be removed with floor/ceil.
    X_trimmed(1:floor(n_timepts_to_remove/2), :) = [];
    X_trimmed(end - (ceil(n_timepts_to_remove/2) - 1):end, :) = [];  
    n_leading_timepts_removed = floor(n_timepts_to_remove/2);
    n_trailing_timepts_removed = ceil(n_timepts_to_remove/2);
else
    error('Invalid value for excess_at.')
end

assert(size(X_trimmed, 1) == n_timepts_to_use);

end
