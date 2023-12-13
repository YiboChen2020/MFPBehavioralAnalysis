function [interactions, labels, second_subj, intro_epoch_ind, rmv_epoch_ind]...
    = isolate_interactions_3(timeseries, session_raw, dict, varargin)
% timeseries can be 2D or 3D. Last dimension is always time. Pass in dict
% as empty to skip generating labels. Otherwise pass in a map/container
% structure.
% interactions is the timeseries intercepted in the shape of (n_interaction, timeseries 2D or 3D )
% labels is the mark of behavior according to dict, in the shape of (n_interacion, length(interactions{i}))
% second_subj is the subject of an ineraction, in the shape of (n_interaction, 1)

%% Parse inputs

p = inputParser;
addRequired(p, 'timeseries');
addRequired(p, 'session_raw');
addRequired(p, 'dict');
addParameter(p, 'allowed_second_subjects', {'M', 'F'}); % Must be {'M', 'F', 'Pup', 'Toy'} or some subset thereof
addParameter(p, 'include_intro_rmv', true);
parse(p, timeseries, session_raw, dict, varargin{:});
nv = p.Results;


%% Find epochs corresponding to two subjects (as defined by 'second_subject')

[intro_epoch_ind, rmv_epoch_ind] ...
    = get_interaction_epoch_ind(session_raw, nv.allowed_second_subjects);


%% Extract epochs from timeseries

% Dimensionality of input array.
n_array_dims = length(size(timeseries));

% Number of second subjects.
n_interactions = length(intro_epoch_ind);

% Extract portion of timeseries labeled as having two subjects present
% (optionally including toy as a subject).
interactions = cell(n_interactions, 1);
labels = cell(n_interactions, 1);
second_subj = cell(n_interactions, 1);
for i_interaction = 1:n_interactions

    % Set whether to use the start point of the epoch
    % (annotated behavior) following Intro as the start point, and the
    % endpoint of the epoch preceding Rmv as the endpoint, OR the start
    % point of the Intro and endpoint of the Rmv epochs as the start and
    % end points, respectively.
    if nv.include_intro_rmv, offset = 0; else, offset = 1; end

    % Get frame indices corresponding to beginning and end of current
    % interaction. 
    curr_epoch_ind = session_raw.Fstart(intro_epoch_ind(i_interaction)+offset) ...
                 :session_raw.Fstop(rmv_epoch_ind(i_interaction)-offset);

    % Get portion of timeseries corresponding to current epoch.
    if n_array_dims == 2
        curr_interaction = timeseries(:,curr_epoch_ind);
    elseif n_array_dims == 3
        curr_interaction = timeseries(:,:,curr_epoch_ind);
    end

    % Assign into cell.
    interactions{i_interaction} = curr_interaction;

    % Determine identity of current 2nd subject.
    second_subj{i_interaction} = session_raw.behaviors{intro_epoch_ind(i_interaction)}(7:end);

    % ---------------------------------------------------------------------
    % Optionally generate labels on all timepoints.
    if ~isempty(dict)
         % A vector of indices into the session_raw.behavior cell array.
        interaction_ind ...
            = intro_epoch_ind(i_interaction)+offset:rmv_epoch_ind(i_interaction)-offset;
    
        % Number of timepoints in current interaction.
        n_timepoints_curr = size(curr_interaction, 2);
        labels{i_interaction} = nan(n_timepoints_curr, 1);
    
        % Behavior (annotation) list from current interaction.
    %     curr_bhvs = session_raw.behaviors(intro_epoch_ind(i_interaction)+offset ...
    %                                       : rmv_epoch_ind(i_interaction)-offset);
        curr_bhvs = session_raw.behaviors(interaction_ind);
    
        for i_epoch = 1:length(curr_bhvs)
            if i_epoch == 1
                i_start = 1;
            else
                i_start = i_end + 1;
            end
    
            % Get index for end of current epoch by adding to the start index
            % the length of the current epoch.
            i_end = i_start + (session_raw.Fstop(interaction_ind(i_epoch)) ...
                               - session_raw.Fstart(interaction_ind(i_epoch)));
    
            % Assign into labels vector.
            labels{i_interaction}(i_start:i_end) = dict(curr_bhvs{i_epoch});
    
            if i_epoch == length(curr_bhvs), assert(i_end == length(labels{i_interaction})); end
        end
    end
end


