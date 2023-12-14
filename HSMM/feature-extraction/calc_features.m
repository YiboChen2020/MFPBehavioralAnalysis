function [F, feat_names, features] = calc_features(T1, T2, tracked_var_names, varargin)
% 
% Generates suite of behavioral features based on 2D tracking data of two
% subjects.
%
% PARAMETERS
% ----------
% T1                : m x n x p array, where m = number of tracked body
%                     parts, n = number of spatial dimensions (2 because
%                     video is top-down hence 2D), and p = number of
%                     timepoints. The (i_th,1,k_th) element of T1 is the X
%                     position of the i_th body part of subject 1 at the
%                     k_th timepoint and the (i_th,2,k_th) element of T1 is
%                     the Y axis of the same body part at the same point in
%                     time.
% T2                : Same as T1 but for subject 2.
% tracked_var_names : Cell array of length c where each cell contains the
%                     name of a body part. The order of the names must
%                     match the order in which they appear along the first
%                     dimension of T1 and T2; it is through this argument
%                     that the function knows which slice of T1 (i, :, :)
%                     corresponds to which body part.
% Name-Value Pairs (nvp)
%   'frame_rate'       : Sampling frequency (Hz). Used to convert duration
%                        of time into number of timepoints.
%   'features_to_keep' : {'all' (default) | 'both' | 'subj1' | 'subj2'}. If
%                        'all', all features are retained. If 'both' only
%                        features computed based on the 2 subjects are
%                        retained. If 'subj1' or 'subj2', only features
%                        computed soley based on that respective subject
%                        are retained.
%
% RETURNS
% -------
% F          : f x p array, where f = number of features and p = number of
%              timepoints. The i_th,j_th element is the value of the i_th
%              feature at the j_th timepoint.
% feat_names : f x 1 cell array whose i_th element is a char array
%              consisting of the name of the i_th feature.
% features   : Scalar struct where each field consists of a single feature
%              (fieldnames match those in feat_names).
%
% Author: Jonathan Chien 2022


p = inputParser;
addRequired(p, 'T1');
addRequired(p, 'T2');
addRequired(p, 'tracked_var_names');
addParameter(p, 'frame_rate', 25);
addParameter(p, 'features_to_keep', 'all')
parse(p, T1, T2, tracked_var_names, varargin{:});
nv = p.Results;


%% Prepare for feature calculation

% Check inputs.
assert(all(size(T1) == size(T2)), ...
       "Array shapes for the two subjects do not match.")
n_tracked_vars = size(T1, 1);
assert(n_tracked_vars == length(tracked_var_names), ...
       "The number of tracked variable names does not match the number " + ...
       "of tracked variables supplied.")

% Create dictionary mapping tracked variable names to indices along first
% array dimension.
dict = create_dict(tracked_var_names);

% Determine n values for "last n" ms features (e.g. tail movement in the
% last 500 ms at each timepoint).
window_size(1) = ms_to_timepts(166, nv.frame_rate);
window_size(2) = ms_to_timepts(500, nv.frame_rate);

% Initialize.
features = struct();


%% Calculate features based on euclidean dist between two body parts
% Body parts may come from same subject or one body part from each subject.
% Note that features where both body parts come from the same subject are
% listed with that subject as part of the feature name (e.g.,
% subj1_ear_to_ear is the distance between ears of subject 1). Otherwise,
% the first body part listed comes from subject 1, while the second body
% part comes from subject 2 (e.g., nose_to_lat_left indicates distance from
% subject 1 nose to subject 2 lateral left, while lat_left_to_nose
% indicates the distance from subject 1 lateral left to subject 2 nose).

% Mounting 3.
features.both_nose_to_tail_base_dist ...
    = two_point_dist(T1, T2, dict, 'nose', 'tail base'); 

% Anogenital sniff 4, lateral threat 3.
features.both_nose_to_nose_dist ...
    = two_point_dist(T1, T2, dict, 'nose', 'nose'); 

% Lateral threat 2, Flee 1, Anogenital sniff 2, Allogrooming vigorous 2,
% Allogrooming normal 2.
features.both_nose_to_lat_left ...
    = two_point_dist(T1, T2, dict, 'nose', 'lateral left'); 

% Lateral threat 1, Flee 2, Anogenital sniff 3, Pursuit 3 (presumably nose
% to left is also good for pursuit too), Allogrooming vigorous 1,
% Allogrooming normal 1.
features.both_nose_to_lat_right...
    = two_point_dist(T1, T2, dict, 'nose', 'lateral right'); 

% This is from Bing's function and not in the paper.
features.both_lat_left_to_nose ...
    = two_point_dist(T1, T2, dict, 'lateral left', 'nose');

% This is from Bing's function and not in the paper.
features.both_lat_right_to_nose ...
    = two_point_dist(T1, T2, dict, 'lateral right', 'nose');

% Lateral threat 4.
features.subj1_ear_to_ear ...
    = two_point_dist(T1, T1, dict, 'left ear', 'right ear');

% Sum of subject 2 width in last 500 ms (Upright submissive 3).
features.subj2_sum_last_width ...
    = last_n_moving(two_point_dist(T2, T2, dict, 'lateral left', 'lateral right'), ...
                    window_size(1), 'sum');

% Sum of subject 1 width in last 500 ms (Scramble 3).
features.subj1_sum_last_width ...
    = last_n_moving(two_point_dist(T1, T1, dict, 'lateral left', 'lateral right'), ...
                    window_size(1), 'sum');


%% Convex hull features

% Calculate convex hull and min/max distances.
conv_hull_1 = convex_hull(T1);
conv_hull_2 = convex_hull(T2);

% Median of longest distance between any two points in convex hull in the
% last 166 ms (Pursuit 1).
features.subj1_median_last_longest_hull ...
    = last_n_moving(conv_hull_1.max_dist, window_size(1), 'median');

% Same as below but with bin size of 166 ms instead of 500 ms.
% (Standardized) deviation from subject 2 mean shortest hull distance in
% moving 166 ms bins centered at each timepoint (Mounting 4).
features.subj2_meandev_bin1_shortest_hull ...
    = moving_bin_stats(conv_hull_2.min_dist, window_size(1), 'mean');

% (Standardized) deviation from subject 2 mean shortest hull distance in
% moving 500 ms bins centered at each timepoint (Mounting 2).
[~, features.subj2_meandev_bin2_shortest_hull] ...
    = moving_bin_stats(conv_hull_2.min_dist, window_size(2), 'mean');

% Percentile rank of subject 2 shortest hull dsitance in 500 ms bins
% (Mounting 1).
features.subj2_prctile_rank_mean_bin_shortest_hull ...
    = calc_prctile_rank(moving_bin_stats(conv_hull_2.min_dist, window_size(2), 'min'));


%% Body part distance and movement features

% Calculate body part metrics.
bp_1 = calc_body_parts(T1);
bp_2 = calc_body_parts(T2);

% Subject 1 centroid X position.
features.subj1_centroid_X_pos = bp_1.centroid(1,:);

% Subject 1 centroid Y position.
features.subj1_centroid_Y_pos = bp_1.centroid(2,:);

% Subject 2 centroid X position.
features.subj2_centroid_X_pos = bp_2.centroid(1,:);

% Subject 2 centroid Y position.
features.subj2_centroid_Y_pos = bp_2.centroid(2,:);

% Percentile rank of distance between subject centroids (Flee 4).
features.both_prctile_rank_centroid_dist ...
    = calc_prctile_rank(euclidean_distance(bp_1.centroid, bp_2.centroid));

% Percentile rank of sum of body part distances (Flee 3).
features.both_prctile_rank_sum_bp_dist ...
    = calc_prctile_rank(sum(euclidean_distance(T1, T2)));

% Median centroid movement in last 500 ms (Pursuit 4).
features.subj1_median_last_centroid_mvmt ...
    = last_n_moving(calc_mvmt(bp_1.centroid), window_size(2), 'median');

% Tail end was not tracked.
% % Median tail end movement in the last 500 ms (Attack 1).
% features.subj1_median_last_tail_base_mvmt ...
%     = last_n_moving(bp_1.mvmt.mvmt(dict('tail end'),:), last_n(2), 'median');

% Median tail base movement in the last 500 ms (Attack 4).
features.subj1_median_last_tail_base_mvmt ...
    = last_n_moving(bp_1.mvmt.mvmt(dict('tail base'),:), window_size(2), 'median');

% Sum of shortest distances in subject 2 within the last 166 ms (Upright
% submissive 4).
features.subj2_sum_last_shortest_bp_dist ...
    = last_n_moving(bp_2.dist.min, window_size(1), 'sum');

% Mean of shortest distances in subject 2 within the last 166 ms (Upright
% submissive 1).
features.subj2_mean_last_shortest_bp_dist ...
    = last_n_moving(bp_2.dist.min, window_size(1), 'mean');
 
% Mean movement of all body parts (assuming among both subjects?) in last
% 166 ms (Scramble 4).
features.both_mean_last_all_bp_mvmt ...
    = last_n_moving(mean(calc_mvmt([T1; T2])), window_size(1), 'mean');

% Percentile rank of mean (over last 500 ms) of sum of distances (Scramble
% 2; note that this feature is similar to that of Flee 3).
features.both_prctile_rank_mean_last_sum_bp_dist ...
    = calc_prctile_rank( ...
        last_n_moving(sum(euclidean_distance(T1, T2)), window_size(2), 'mean') ...
                        );                      

% Shortest value (within sliding 500 ms window) for mean body part dist for
% subject 2 (Allogrooming vigorous 4).
features.subj2_shortest_last_mean_bp_dist ...
    = last_n_moving(bp_2.dist.mean, window_size(2), 'min');

% Shortest value (within sliding 500 ms window) for median body part dist
% for subject 2 (Allogrooming vigorous 3).
features.subj2_shortest_last_median_bp_dist ...
    = last_n_moving(bp_2.dist.median, window_size(2), 'min');

% (Standardized) deviation from mean of mean movement of all body parts
% (both subjects) in 500 ms bins (Attack 3).
[~, features.both_devmean_bin_all_body_part_mvmt] ...
    = moving_bin_stats(mean([bp_1.mvmt.mean; bp_2.mvmt.mean]), ...
                       window_size(2), ...
                       'mean');

% Percentile rank of the sum of centroid movements of subject 1 over the
% last 500 ms (Attack 2).
features.subj1_prctile_rank_last_sum_centroid_mvmt ...
    = calc_prctile_rank(last_n_moving(calc_mvmt(bp_1.centroid), window_size(2), 'sum'));


%% Convert to array

[F, feat_names] = to_array(features);


%% Take specified subset of features

% Can retain subset of features (e.g., only those features that were
% calculated solely based on subject 1).
if ~strcmp(nv.features_to_keep, 'all')
    assert(ismember(nv.features_to_keep, {'both', 'subj1', 'subj2'}), ...
           "'features_to_keep' must be 'all' | 'both' | 'subj1' | 'subj2'.")
    del_ind = cellfun(@isempty, regexp(feat_names, nv.features_to_keep));
    features = rmfield(features, feat_names(del_ind));
    feat_names(del_ind) = [];
    F(del_ind,:) = [];
end


end


% --------------------------------------------------
function n_timepoints = ms_to_timepts(duration, frequency)
% Convert a time duration (in ms), given a sampling frequency, into the
% nearest number of timepoints.

n_timepoints = round(duration / frequency);

end


% --------------------------------------------------
function dict = create_dict(tracked_var_names)
% Create dictionary mapping passed in names of tracked variables to their
% corresponding first dimension index in the passed in array of tracked
% variables in the time domain. Arrays can be passed in with any order of
% the tracked variables, and subroutines called, so long as the string name
% of the tracked variables are preserved.

n_vars = length(tracked_var_names);
dict = containers.Map(tracked_var_names, 1:n_vars);

end


% --------------------------------------------------
function [S1, S2] = select_var(T1, T2, dict, var1, var2)
% Use dictionary to get 2D arrays corresponding to X, Y coords of specified
% tracked variable, for both subjects. T1 and T2 are n_vars x 2 x
% n_timepoints arrays, and S1 and S2 are 2 x n_timepoints arrays
% (corresponding to a slice across dims 2 and 3 of T1 and T2).

S1 = squeeze(T1(dict(var1),:,:));
S2 = squeeze(T2(dict(var2),:,:));

end


% --------------------------------------------------
function D = euclidean_distance(S1, S2)
% If S1 and S2 are 2 x n_timepoints (first row corresponds to X coords and
% second row to Y coords),  D is a 1 x n_timepoints vector of distances. If
% S1 and S2 are 2 x m x n_timepoints arrays, (where first array dim
% contains x and y coordinates), D is an m x n_timepoints matrix of
% distances.

if length(size(S1)) == 3 && length(size(S2)) == 3
    S1 = permute(S1, [2 1 3]);
    S2 = permute(S2, [2 1 3]);
end 

D = sqrt(sum((S2 - S1) .^ 2));
D = squeeze(D);

end


% --------------------------------------------------
function D = two_point_dist(T1, T2, dict, var1, var2)
% D is a 1 x n_timepoints vector.

[S1, S2] = select_var(T1, T2, dict, var1, var2);
D = euclidean_distance(S1, S2);

end


% --------------------------------------------------
function M = calc_mvmt(A)
% A is either an m x 2 x n_timepoints array, where m corresponds to an
% arbitrary number of variables, or a 2 x n_timepoints array (for just one
% variable). For 3D inputs, the i_th,j_th element of M is the movement of
% the i_th variable from the j-1_th timepoint to the j_th timepoint; for 2D
% inputs, the i_th element of M (M is a row vector here) is the movement
% from the i-1_th timepoint to the i_th timepoint. The first slice along
% the first array dimension of M in either case contains all NaNs as
% padding.

if length(size(A)) == 2
    A = reshape(A, [1 size(A)]);
end
n_vars = size(A, 1);
M = sqrt(sum(diff(A, 1, 3).^2, 2));
M = squeeze(M); if iscolumn(M), M = M'; end
M = [NaN(n_vars, 1), M];

end


% --------------------------------------------------
function M = last_n_moving(A, last_n, metric)
% Use sliding window over last n points to calculate moving mean, median,
% or sum along second dimension of A, an m x n_timepoints array, where m is
% an arbitary number of variables. 

switch metric
    case 'mean'
        M = movmean(A, [last_n 0], 2);
    case 'median'
        M = movmedian(A, [last_n 0], 2);
    case 'sum'
        M = movsum(A, [last_n 0], 2);
    case 'min'
        M = movmin(A, [last_n 0], 2);
end

% Disqualify edge cases.
M(:,1:last_n) = NaN;

end

% --------------------------------------------------
function [M, D] = moving_bin_stats(A, n, metric)
% A is an m x n_timepoints array, where m is an arbitrary number of
% varaibles. A sliding window of size n (units are in timepoints) is passed
% along the second dimension of A; the window is stepped at 1 timepoint and
% centers on each timepoint if n_timepoints is odd and on the current and
% previous timepoint if n_timepoints is even. For each timepoint, first and
% second order sample statistics are calculated for the window current
% timepoint. Edge cases are handled by averaging over available elemnts
% (NaN padding + mean with omitnan option).

% Calculate statistic.
switch metric
    case 'mean'
        M = movmean(A, n, 2);
    case 'median'
        M = movmedian(A, n, 2);
    case 'min'
        M = movmin(A, n, 2);
end

% Disqualify edge cases (handle based on parity of window size).
if mod(n, 2) == 0
    M(:,1:n/2+1) = NaN;
    M(:,end-(n/2-1):end) = NaN;
elseif mod(n, 2) ~= 0
    M(:,1:ceil(n/2)) = NaN;
    M(:,end-(ceil(n/2)-1)) = NaN;
end

% Calculate standardized deviation.
D = (A - M) ./ M;

end


% --------------------------------------------------
function C = calc_prctile_rank(A)
% A is an m x n_timepoints array, where m corresponds to an arbitary number
% of variables. The i_th, j_th element of C is the proportion of elements
% of the i_th row (i_th variable) that fall below the value of the i_th
% variable at the j_th timepoint.

[n_vars, n_timepoints] = size(A);
C = NaN(n_vars, n_timepoints);
parfor i_timept = 1:n_timepoints
    C(:,i_timept) = sum(A(:,i_timept) > A, 2) / n_timepoints;
end

end


% --------------------------------------------------
function conv_hull = convex_hull(T)

% Get number of timepoints.
n_timepoints = size(T, 3);

% Preallocate. Note that convex hull may be defined by a different number
% of points at different timepoints.
conv_hull.conv_hull = cell(1, n_timepoints);
conv_hull.max_dist = NaN(1, n_timepoints);
conv_hull.min_dist = NaN(1, n_timepoints);

% At each timepoint, find convex hull and max/min distance.
for i_timept = 1:n_timepoints
    % Slice corresponding to current timepoint.
    curr_T = T(:,:,i_timept);

    % Convex hull at current timepoint.
    conv_hull.conv_hull{i_timept} = curr_T(unique(convhull(curr_T)),:);
    
    % Calculate pairwise distances.
    D = pdist(conv_hull.conv_hull{i_timept}, 'euclidean');

    % Distance between furthest points in convex hull.
    conv_hull.max_dist(i_timept) = max(D);

    % Distance between closest points in convex hull.
    conv_hull.min_dist(i_timept) = min(D);
end

end


% --------------------------------------------------
function body_parts = calc_body_parts(T)

[n_body_parts, ~, n_timepoints] = size(T);

% Calculate pairwise distances.
D = NaN(n_body_parts*(n_body_parts-1)/2, n_timepoints);
for i_timept = 1:n_timepoints
    D(:,i_timept) = pdist(T(:,:,i_timept));
end

% Mean, sum, max, min of body part distances at each timepoint.
body_parts.dist.mean = mean(D, 1);
body_parts.dist.median = median(D, 1);
body_parts.dist.sum = sum(D, 1);
body_parts.dist.sum_prctile_rank = calc_prctile_rank(body_parts.dist.sum);
body_parts.dist.max = max(D);
body_parts.dist.max_prctile_rank = calc_prctile_rank(body_parts.dist.max);
body_parts.dist.min = min(D);
body_parts.dist.min_prctile_rank = calc_prctile_rank(body_parts.dist.min);

% Centroid of body parts.
body_parts.centroid = squeeze(mean(T, 1));
assert(size(body_parts.centroid, 1) == 2);

% Movement statistics.
body_parts.mvmt.mvmt = calc_mvmt(T); % Movement of each body part
body_parts.mvmt.mean = mean(body_parts.mvmt.mvmt); % Mean across body parts
body_parts.mvmt.std = std(body_parts.mvmt.mvmt, 0, 1); % Std across body parts

end


% --------------------------------------------------
function [F, feat_names] = to_array(features)
% Takes in struct with each feature in a field and returns feature matrix
% (each feature as a row) as well as an accompanying ordered cell array of
% feature names.

% Get feature names.
feat_names = fieldnames(features);

% Get feature matrix.
n_features = length(feat_names);
n_timepoints = size(features.(feat_names{1}), 2);
F = NaN(n_features, n_timepoints);
for i_feature = 1:n_features
    F(i_feature,:) = features.(feat_names{i_feature});   
end

end
