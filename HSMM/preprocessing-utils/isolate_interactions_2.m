function [interactions, labels, n_pts_rmv, intro_epoch_ind, rmv_epoch_ind]...
    = isolate_interactions_2(timeseries, session_raw, dict, varargin)
% timeseries can be 2D or 3D. Last dimension is always time.


%% Parse inputs

p = inputParser;
addRequired(p, 'timeseries');
addRequired(p, 'session_raw');
addRequired(p, 'dict');
addParameter(p, 'second_subject', {'M', 'F'}); % Must be {'M', 'F', 'Pup', 'Toy'} or some subset thereof
addParameter(p, 'human_annotations', true);
addParameter(p, 'smooth', false)
addParameter(p, 'bin_width', 1)
addParameter(p, 'step_size', 1)
parse(p, timeseries, session_raw, dict, varargin{:});
nv = p.Results;


%% Find epochs corresponding to two subjects (as defined by 'second_subject')

[intro_epoch_ind, rmv_epoch_ind] ...
    = get_interaction_epoch_ind(session_raw, nv.second_subject);


%% Extract epochs from timeseries

% Dimensionality of input array.
n_array_dims = length(size(timeseries));

% Number of second subjects.
n_interactions = length(intro_epoch_ind);

% Extract portion of timeseries labeled as having two subjects present
% (optionally including toy as a subject).
interactions = cell(n_interactions, 1);
labels = cell(n_interactions, 1);
n_pts_rmv = cell(n_interactions, 1);
for i_interaction = 1:n_interactions
    % Intialize for concatenating across epochs (within interaction).
    interactions{i_interaction} = [];
    labels{i_interaction} = [];
    n_pts_rmv{i_interaction} = [];

    % Here, within the current interaction, we want to extract portions of
    % timeseries with human labels (not 'Other'). Each such portion is
    % termed an epoch.
    if nv.human_annotations
        % A vector of indices into the session_raw.behavior cell array.
        interaction_ind ...
            = intro_epoch_ind(i_interaction):rmv_epoch_ind(i_interaction);

        % For each annotated period (not 'Other'), apply any potential
        % binning or smoothing and return the result. Note that this
        % includes the Intro and Remove epochs.
        for i_epoch = 1:length(interaction_ind)
            % Skip if behavior is not one of the annotated and included ones.
            if ~ismember(session_raw.behaviors(interaction_ind(i_epoch)), keys(dict))
                continue
            end

            % Get start and end timepoints of current epoch.
            curr_epoch_ind ...
                = session_raw.Fstart(interaction_ind(i_epoch)) ...
                  :session_raw.Fstop(interaction_ind(i_epoch));

            % Get portion of timeseries corresponding to current epoch.
            if n_array_dims == 2
                curr_epoch = timeseries(:,curr_epoch_ind);
            elseif n_array_dims == 3
                curr_epoch = timeseries(:,:,curr_epoch_ind);
            end
            
            % Optionally preprocess portion of timeseries
            % associated with current epoch.
            [binned_epoch, n_pts_rmv_for_binning] ...
                = smooth_and_bin(curr_epoch, ...
                                 nv.smooth, nv.bin_width, nv.step_size);

            % Concatenate.
            interactions{i_interaction} ...
                = cat(n_array_dims, interactions{i_interaction}, binned_epoch);
            labels{i_interaction} = [labels{i_interaction}; ...
                               repelem(dict(session_raw.behaviors{interaction_ind(i_epoch)}), ...
                                       size(binned_epoch, n_array_dims), ...
                                       1)];
            n_pts_rmv{i_interaction} ...
                = [n_pts_rmv{i_interaction}; n_pts_rmv_for_binning];
        end
    else
        % Treat entire interaction period as one epoch and get frame
        % indices corresponding to beginning and end. NOTE: here, we use
        % the start point of the epoch following Intro as the start point,
        % and the endpoint of the epoch preceding Rmv as the endpoint.
        curr_epoch_ind = session_raw.Fstart(intro_epoch_ind(i_interaction)+1) ...
                     :session_raw.Fstop(rmv_epoch_ind(i_interaction)-1);

        % Get portion of timeseries corresponding to current epoch.
        if n_array_dims == 2
            curr_epoch = timeseries(:,curr_epoch_ind);
        elseif n_array_dims == 3
            curr_epoch = timeseries(:,:,curr_epoch_ind);
        end

        % Optionally preprocess entire interaction as one epoch.
        [interactions{i_interaction}, n_pts_rmv_for_binning] ...
            = smooth_and_bin(curr_epoch, ...
                             nv.smooth, nv.bin_width, nv.step_size);

        % Labels cannot be generated here since bins may overlap between
        % differen annotations.

        % Numbers of points (leading and trailing) dropped for binning.
        n_pts_rmv{i_interaction} = n_pts_rmv_for_binning;
    end
end

end


% --------------------------------------------------
function [epoch_timeseries_binned, n_pts_rmv] ...
    = smooth_and_bin(epoch_timeseries, smooth, bin_width, step_size)

% Dimensionality of array.
n = length(size(epoch_timeseries));

% Optionally smooth data.
if (smooth > 0) && (smooth < 1)
    % Set window width for smoothing as supplied fraction of length of
    % current annotated section (this results in different timescales
    % of smoothing across the annotated sections).
    epoch_timeseries = smoothdata(epoch_timeseries, n, 'gaussian', ...
                                  size(epoch_timeseries, n) * smooth);
elseif smooth
    % Smooth using supplied window width.
    epoch_timeseries = smoothdata(epoch_timeseries, n, 'gaussian', smooth);
end

% Bin data. Add number of timepoints retained from current epoch to
    % tally across all epochs.
    [epoch_timeseries_binned, ~, n_leading_pts_rmv, n_trailing_pts_rmv] ...
        = bin_nd_data(epoch_timeseries, n, bin_width, step_size);

% Arrange number of leading and trailing points dropped as elements of
% length 2 row vector.
n_pts_rmv = [n_leading_pts_rmv n_trailing_pts_rmv];

end
