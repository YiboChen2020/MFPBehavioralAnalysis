function mpp = marked_point_process(timeseries, varargin)
% Accepts as input an n_regions x n_timepoints array of neural activity and
% finds all peaks above specified threshold (as percent of max activity).
% Optionally convolve filter across resulting peaks.
%
% PARAMETERS
% timeseries : m x n array of neural data where m is the number of regions
%              and n the number of timepoints.
% Name-Value Pairs (nv)
%   'min_threshold' : (Scalar, default = 0.25) Threshold (proportion) of
%                     signal max (by region) below which all values will be
%                     set to 0.
%   'zscore'        : (1 (default) | 0), specify whether or not to zscore
%                     data.
%   'convolve'      : (vector | false), pass filter in as a vector or set
%                     to false to skip convolution.
%   'plot'          : (1 | 0 (default)), specify whether or not to
%                     visualize resulting point process.
%   'window'        : 2-vector indexing the beginning and ending timepoints
%                     of the period of the timeseries that we would like to
%                     visualize.
%   'figure'        : (1 | 0 (default)), specify whether or not to plot in 
%                     new figure.
%
% RETURNS
% -------
% mpp : m x n array matching shape of input argument "timeseries". Each row
%       is a (filtered) marked point process.


%% Parse inputs.

p = inputParser;
addRequired(p, 'timeseries');
addParameter(p, 'min_threshold', 0.25, @(x) (isscalar(x) && 0 <= x && x <= 1) || isinf(x));
addParameter(p, 'zscore', true, @islogical); % Optionally zscore first.
addParameter(p, 'convolve', [0.1353	0.3679	1	0	0]); % [0.1 0.2 0.4 0 0]/0.7)
addParameter(p, 'plot', false, @islogical);
addParameter(p, 'window', [3000 7000]);
addParameter(p, 'figure', false, @islogical);
parse(p, 'timeseries', varargin{:});
nv = p.Results;


%% Peak detection

% Option to z-score data first.
if nv.zscore, timeseries = zscore(timeseries, 0, 2); end

n_regions = size(timeseries, 1);
mpp = zeros(size(timeseries));
for i_region = 1:n_regions
    % Find threshold for current region.
    threshold = nv.min_threshold * max(timeseries(i_region,:));

    % Find peaks.
    [peaks, ind] = findpeaks(timeseries(i_region,:), 'MinPeakHeight', threshold);
    mpp(i_region,ind) = peaks;
end

% Optionally convolve to create filtered marked point process.
if any(nv.convolve)
    for i_region = 1:n_regions
        mpp(i_region,:) = conv(mpp(i_region,:), nv.convolve, 'same');
    end
end


%% Optional plotting

if nv.plot
    if nv.figure, figure; end

    for i_region = 1:n_regions
        subplot(n_regions, 1, i_region)
        stem(mpp(i_region,nv.window(1):nv.window(2)), 'Marker', 'none')
    end
end

end