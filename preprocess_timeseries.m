function neural_data = preprocess_timeseries(session_raw, data_type, params)

% p = inputParser;
% addRequired(p, 'session_raw');
% addRequired(p, 'data_type');
% addRequired(p, 'params');
% parse(p, session_raw, data_type, params, varargin{:});
% nv = p.Results;


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
            = marked_point_process(session_raw.Lfold, ...
                                   'min_threshold', params.min_threshold, ...
                                   'zscore', params.zscore, ...
                                   'convolve', params.convolve, ...
                                   'plot', false);
    case 'raw'
        neural_data = session_raw.Lfold;
    otherwise
        error("Please provide a valid value for 'data_type'.")
end


%% Optional smoothing

if params.smooth
    neural_data = smoothdata(neural_data, 2, 'movmean', params.smooth);
end


% %% Optional resampling
% 
% if params.resample_params.resample
%     % TODO: Move this routine into the isolate_behavior and isolate_baseline
%     % functions.
%     
%     n_regions = size(neural_data, 1);
%     
%     % Get variances for Poisson distribution.
%     if strcmp(params.resample_params{2}, 'scaled')
%         lambdas = get_scaled_lambdas(neural_data, params.resample_params{1});
%     elseif strcmp(params.resample_params{2}, 'flat')
%         lambdas = ones(n_regions, 1) * params.resample_params{1};
%     elseif strcmp(params.resample_params{2}, 'specified')
%         assert(isvector(params.resample_params{1}) ...
%                && length(params.resample_params{1}) == n_regions, ...
%                "If the second element of 'resample_params' is 'specified' " + ...
%                "the first element must be a vector matching the number " + ...
%                "of regions in length.")
%         lambdas = params.resample_params{1};
%     end
%     
%     % Resample for each region.
%     for i_region = 1:n_regions
%         neural_data(i_region,:) ...
%             = resample2poisson(neural_data(i_region,:), ...
%                                lambdas(i_region), ...
%                                params.resample_seed);
%     end
% end

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