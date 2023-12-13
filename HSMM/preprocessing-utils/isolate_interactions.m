function [interactions, intro_epoch_ind, rmv_epoch_ind]...
    = isolate_interactions(timeseries, session_raw, varargin)
% timeseries can be 2D or 3D. Third dimension is always time, and the
% epochs are indexed along this dimension.


%% Parse inputs

p = inputParser;
addRequired(p, 'timeseries');
addRequired(p, 'session_raw');
addParameter(p, 'second_subject', {'M', 'F'}); % Must be {'M', 'F', 'Pup', 'Toy'} or some subset thereof
addParameter(p, 'legacy', false); % If true, indexes from start of intro to end of remove. If false, indexes from end of intro to start of remove.
parse(p, timeseries, session_raw, varargin{:});
nv = p.Results;


%% Find epochs corresponding to two subjects (as defined by 'second_subject')

[intro_epoch_ind, rmv_epoch_ind] ...
    = get_interaction_epoch_ind(session_raw, nv.second_subject);


%% Extract epochs from timeseries

% Dimensionality of input array.
n_array_dims = length(size(timeseries));

% Extract portion of timeseries labeled as having two subjects present
% (optionally including toy as a subject).
n_epochs = length(intro_epoch_ind);
interactions = cell(n_epochs, 1);
for i_epoch = 1:n_epochs
    % Vector indexing timepoints of current epoch.
    if nv.legacy
        % This was the older (and original) way of indexing. We take the
        % beginning of the intro epoch to the end of the remove epoch.
        curr_epoch = session_raw.Fstart(intro_epoch_ind(i_epoch)) ...
                     :session_raw.Fstop(rmv_epoch_ind(i_epoch));
    else
        % Newer way of indexing taking from end of intro epoch to beginning
        % of remove epoch.
        curr_epoch = session_raw.Fstop(intro_epoch_ind(i_epoch)) ...
                     :session_raw.Fstart(rmv_epoch_ind(i_epoch));
    end

    % If input is 2D array (usually if this is neural data : n_regions x
    % n_timepoints).
    if n_array_dims == 2
        interactions{i_epoch} = timeseries(:,curr_epoch);
    % If input array is 3D (usually if this is tracking data:
    % n_body_parts_tracked x 2 x n_timepoints, with 2 bc data are in 2D (X
    % and Y) coordinates).
    elseif n_array_dims == 3
        interactions{i_epoch} = timeseries(:,:,curr_epoch);
    end
end


end
