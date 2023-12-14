function [binned_timeseries, bhv_labels, aux] ...
    = labeled_timeseries(timeseries, sessions_raw, varargin)
% Accepts a cell array of sessions for one subject and bins data with
% accompnaying behavioral labels. Note: as of 12/2023, this function is
% mostly unused except to get behavioral labels in order to estimate things
% like mean duration, transition frequencies etc.

% PARAMETERS
% timeseries   : n_sessions x 1 cell array whose i_th cell contains
%                an m x n numeric array whose j_th, l_th element is the
%                value of the j_th variable (e.g. brain region) at the n_th
%                timepoints. 
% sessions_raw : n_sessions x 1 cell array whose i_th cell contains session
%                info for the i_th session.
% Name-Value Pairs 
%
% RETURNS
% -------
% binned_timeseries  : n_sessions x 1 cell array whose i_th cell contains
%                      an m x k numeric array whose j_th, l_th element is
%                      the value of the j_th variable (e.g. brain region)
%                      in the l_th bin. See 'combine_sessions' name-value
%                      pair.
% bhv_labels         : n_sessions x 1 cell array whose i_th cell contains
%                      an m-vector whose j_th element is the behavioral
%                      label of the j_th timepoint in the i_th session. See
%                      'combine sessions' name-value pair.
% aux.bhv_list       : list of retained behaviors; these are in
%                      alphabetical and their associated number codes are
%                      also in ascending order (though note that some
%                      bhvs/number codes can be skipped, namely if the
%                      epoch is too short to support the specified bin
%                      width and step size).
% aux.bhv_name_codes : Note that some values may be skipped if all
%                      instances of that behavior were too short to support
%                      the specified bin width and step size.
% aux.remove_epochs  : first column lists behaviors where (at least once)
%                      more consecutive timebins with all zeros (than the
%                      threshold value for number of bins) were identified
%                      and removed. Second column lists the number of times
%                      this occurs for each behavior (not the number of
%                      timebins but the number of times that a consecutive
%                      group of timebins exceeding threshold length was
%                      identified).
%
% Author: Jonathan Chien

% Parse inputs.
p = inputParser;
addRequired(p, 'timeseries');
addRequired(p, 'sessions_raw');
addParameter(p, 'smooth', {'movmean', 5}, @(x) iscell(x) || (islogical(x) && ~x));
addParameter(p, 'bin_width', 10);
addParameter(p, 'step_size', 5);
addParameter(p, 'remove_thresh', false);
addParameter(p, 'suppress_thresh_warning', true);
addParameter(p, 'combine_sessions', false);
parse(p, timeseries, sessions_raw, varargin{:});
nv = p.Results;

% Preallocate.
n_sessions = length(sessions_raw);
binned_timeseries = cell(n_sessions, 1);
bhv_labels = cell(n_sessions, 1);
bhv_list = cell(n_sessions, 1);
remove_epochs = cell(n_sessions, 1);

% Create dictionary for all behaviors across all sessions.
[dict, r_dict, bhv_name_codes] = gen_dict(sessions_raw);

for i_sess = 1:n_sessions
    % Get current session data.
    sess = sessions_raw{i_sess};

    % Optionally smooth entire timeseries first.
    if iscell(nv.smooth)
        timeseries{i_sess} = smoothdata(timeseries{i_sess}, 2, nv.smooth{1}, nv.smooth{2});
    end

    % Number of epochs (different observed behavior) and timepoints for
    % current session.
    n_epochs = length(sess.behaviors);
    n_timepoints = size(timeseries{i_sess}, 2);

    % For each epoch, bin data and store associated label. 
    for i_epoch = 1:n_epochs
        % Some frame labels occur after the end of the timeseries. Need to
        % check with Dayu about this.
        if all([sess.Fstart(i_epoch) sess.Fstop(i_epoch)] <= n_timepoints)
            % Get epoch and bin the data. Concatenate.
            epoch = timeseries{i_sess}(:, sess.Fstart(i_epoch) : sess.Fstop(i_epoch));
            curr_bhv = sessions_raw{i_sess}.behaviors{i_epoch};

            % Attempt to bin current epoch. This may fail if the epoch is
            % too short; in that case, skip the epoch and continue to the
            % next; else, rethrow any other exceptions.
            try
                % binned_epoch = bin_data(epoch, nv.bin_width, nv.step_size);
                binned_epoch = bin_nd_data(epoch, 2, nv.bin_width, nv.step_size, 'excess_at', 'both_ends', 'bin_excess', false); % Name-value pairs are for backwards compatability with bin_data
            catch exception
                if strcmp(exception.identifier, 'determine_n_bins:no_timepoints_left')
                    continue
                else
                    rethrow(exception);
                end
            end

            % Optionally handle timepoints that consist of all zeros across
            % all regions.
            if nv.remove_thresh
                % Remove timepoints with zeros across all regions. If this
                % results in removal of all timepoints in epoch, skip this
                % behavior and move onto the next.
                remove_ind = sum(binned_epoch) == 0;
                binned_epoch(:,remove_ind) = [];
                if isempty(binned_epoch), continue; end
    
                % If zeros across regions are present among some (but not
                % all) of the timepoints, check if there are consecutive
                % timepoints with all zero values where the number of
                % consecutive such points is greater than the specified
                % threshold. If so, issue warning and mark down these
                % epochs.
                mov_sum = movsum(remove_ind, nv.remove_thresh + 1);                
                if any(mov_sum > nv.remove_thresh)
                    % Mark down current behavior as having had some time bins
                    % removed.
                    remove_epochs{i_sess} = vertcat(remove_epochs{i_sess}, {curr_bhv});
                    if ~nv.suppress_thresh_warning
                        warning("More timebins than the specified threshold " + ...
                                "value of %d were removed during a(n) %s epoch.", ...
                                nv.remove_thresh, curr_bhv);
                    end        
                end
            end

            % Number of bins in current epoch.
            epoch_length = size(binned_epoch, 2);

            % Add current epoch timeseries and behavioral labels.
            binned_timeseries{i_sess} = [binned_timeseries{i_sess} binned_epoch];
            
            bhv_labels{i_sess} = [bhv_labels{i_sess}; ...
                                  repelem(dict(curr_bhv), epoch_length)'];  
        end
    end

    % Use values to get keys corresonding to behaviors that were included.
    % We do this here rather than just grabbing the behavior list stored in
    % sessions_raw because upon attempting to bin, it may come about that
    % all instances of a behavior were to0 short to support the desired bin
    % width and step. In that case, that behavior (and its corresponding
    % timepoints) will be skipped. Thus, after binning, we see what
    % behaviors remain, and save those as the behaviors for the current
    % session.
    codes_retained = unique(bhv_labels{i_sess});
    for i_code = 1:length(codes_retained)
        bhv_list{i_sess} ...
            = vertcat(bhv_list{i_sess}, {r_dict(codes_retained(i_code))});
    end
end

% Optionally concatenate/combine sessions.
if nv.combine_sessions
    % Concatenate timeseries and behavioral labels across sessions.
    binned_timeseries = {cell2mat(binned_timeseries')};
    bhv_labels = {cell2mat(bhv_labels)};
    
    % Pool lists of behaviors, as well as list of behaviors that had
    % timepoints removed.
    bhv_list_ = {};
    remove_epochs_ = {};
    for i_sess = 1:n_sessions
        bhv_list_ = vertcat(bhv_list_, bhv_list{i_sess});
        remove_epochs_ = vertcat(remove_epochs_, remove_epochs{i_sess});
    end

    % Remove redundant items in list of behaviors.
    bhv_list = unique(bhv_list_);

    % Combine lists of behaviors removed due to presence of zeros and get
    % counts.
    remove_epochs = remove_epochs_;
    [remove_epochs, ~, ind] = unique(remove_epochs); 
    count = histc(ind, 1:numel(remove_epochs));
    remove_epochs = horzcat(remove_epochs, ...
                                       mat2cell(count, repelem(1, length(count))));
else 
    for i_sess = 1:n_sessions
        % Remove any redundancies in list of behaviors that had timepoints
        % removed. Also, get count of how many times a group of bins greater
        % than threshold was removed, for each behavior.
        [remove_epochs{i_sess}, ~, ind] = unique(remove_epochs{i_sess}); 
        count = histc(ind, 1:numel(remove_epochs{i_sess}));
        remove_epochs{i_sess} = horzcat(remove_epochs{i_sess}, ...
                                       mat2cell(count, repelem(1, length(count))));
    end  
end

% Gather auxiliary outputs.
aux.bhv_list = bhv_list;
aux.bhv_name_codes = bhv_name_codes;
aux.remove_epochs = remove_epochs;

end


% --------------------------------------------------
function [dict, r_dict, key_value_pairs] = gen_dict(sessions_raw)
% Create a dictionary with all behaviors across sessions as keys and a
% number code as value.

% Gather unique labels from each session.
n_sessions = length(sessions_raw);
all_bhvs = {};
for i_sess = 1:n_sessions
    all_bhvs = vertcat(all_bhvs, unique(sessions_raw{i_sess}.behaviors));
end

% Get unique labels across all sessions.
all_bhvs = unique(all_bhvs);
n_bhvs = length(all_bhvs);

% Set up dictionary.
dict = containers.Map(all_bhvs, 1:n_bhvs);

% Also set up reverse dictionary (where number codes are keys and behaviors
% are values) for easier retrieval of behaviors keys from code values.
r_dict = containers.Map(1:n_bhvs, all_bhvs);

% Get and return keys and values in n_bhvs x 2 cell array (1st column
% behaviors/keys, 2nd column codes/values).
key_set = keys(dict);
value_set = values(dict);
key_value_pairs = horzcat(key_set', value_set');

end
