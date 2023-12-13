function neural_data = preprocess_timeseries_arch1(session_raw, data_type, params, varargin)
% Archived on 8/4/22 upon removal of redundancy between params and 

p = inputParser;
addRequired(p, 'session_raw');
addRequired(p, 'data_type');
addRequired(p, 'params');
addParameter(p, 'resample', true);
addParameter(p, 'resample_params', {5, 'scaled'});
addParameter(p, 'resample_seed', 'default');
parse(p, session_raw, data_type, params, varargin{:});
nv = p.Results;


%% Apply preprocessing routine

switch data_type
    % Deconvolution.
    case 'deconv'
        neural_data ...
            = apply_deconvolution(session_raw, ...
                                  params, ...
                                  false); % Plotting
    % Conversion to filtered Markov point process.
    case 'mpp'
        neural_data ...
            = markov_point_process(session_raw.Lfold, ...
                                   'min_threshold', params.min_threshold, ...
                                   'zscore', params.zscore, ...
                                   'convolve', params.convolve, ...
                                   'plot', false);
    case 'raw'
        neural_data = session_raw.Lfold;
    otherwise
        error("Please provide a valid value for 'data_type'.")
end


%% Optional resampling

n_regions = size(neural_data, 1);

% Get variances for Poisson distribution.
if strcmp(nv.resample_params{2}, 'scaled')
    lambdas = get_scaled_lambdas(neural_data, nv.resample_params{1});
elseif strcmp(nv.resample_params{2}, 'flat')
    lambdas = ones(n_regions, 1) * params.resampling.mean_var;
else
    assert(isvector(nv.resample_params{2}) ...
           && length(nv.resample_params) == n_regions, ...
           "If the second element of 'resample_params' is not a string " + ...
           "'scaled' or 'flat', it must be a vector matching the number " + ...
           "of regions in length.")
    lambdas = nv.resample_params{2};
end

% Resample for each region.
for i_region = 1:n_regions
    neural_data(i_region,:) ...
        = resample2poisson(neural_data(i_region,:), ...
                           lambdas(i_region), ...
                           nv.resample_seed);
end

end


% --------------------------------------------------
function deconv = apply_deconvolution(session_raw, params, plot_flag)
% Wrapper for deconvolution. 

% Preallocate.
deconv = NaN(size(session_raw.Lfold));
n_regions = size(deconv, 1);
        
for i_region = 1:n_regions
    % Deconvolve.
    [~, deconv(i_region,:), ~] ...
        = deconvolveCa(session_raw.Lfold(i_region,:), ...
                       params.kernel_model, params.kernel_params, ...
                       params.method, params.extra_params, ...
                       'window', params.window, ...
                       'sn', params.noise_std, ...
                       'b', params.baseline, ...
                       'lambda', params.lambda, ...
                       'optimize_b', params.optimize_b, ...
                       'optimize_pars', params.optimize_pars, ...
                       'optimize_smin', params.optimize_smin);       
end

% Optional plotting. Helpful for debugging/pausing for visualization.
if plot_flag
    for i_region = 1:n_regions
        % Prepare window.
        if true
            plot_window = [1 size(deconv, 2)];
        else
            plot_window = [2000 5000]; % Set desired window
        end
        
        % Plot deconvolved.
        subplot(n_regions, 1, i_region)
        hold on
        stem(deconv(i_region,:), 'Marker', 'none')
        xlim(plot_window)
    end
end

end