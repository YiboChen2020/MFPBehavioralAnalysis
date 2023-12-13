function timeseries = isolate_all_subjects_sessions_bhvs_slist(path_to_data, slist, bhv_list, nv)

arguments
    path_to_data
    slist
    bhv_list
    nv.remove_nan = true
    nv.nan_warning = true
    nv.zscore = false
    nv.smooth = false
    nv.relabel_behaviors = false
    nv.bin_width = 1
    nv.step_size = 1
    nv.min_length = 1
end

% Get indices of all subjects.
subj_ind = find(cellfun(@isempty, slist) == 0);
if iscolumn(subj_ind), subj_ind = subj_ind'; end

% Get highest subject index, max number of sessions (for any subject) and
% max number of behaviors.
i_last_subj = subj_ind(end);
n_sessions_by_subj = nan(i_last_subj, 1);
for i_subj = 1:i_last_subj
    if ~isempty(slist{i_subj})
        n_sessions_by_subj(i_subj) = length(slist{i_subj});
    end
end
n_sessions_max = max(n_sessions_by_subj, [], 'omitnan');
n_bhvs = length(bhv_list);


timeseries = cell(i_last_subj, n_sessions_max, n_bhvs);
for i_subj = subj_ind
    
    
    % Iterate over all sessions for current subject.
    session_ids = slist{i_subj};
    n_sessions = length(session_ids);
    for i_session = 1:n_sessions
        % Current session ID.
        curr_session_id = session_ids{i_session};

        % Get all sessions for current subject.
        sessions_raw = extract_single_subj(path_to_data, i_subj, ...
                                           'session', session_ids, ...
                                           'combine_sessions', false, ...
                                           'region_names', true, ...
                                           'remove_nan', true, ...
                                           'nan_warning', true, ...
                                           'zscore', nv.zscore, ...
                                           'smooth', nv.smooth, ...
                                           'plot', false);

        % Optionally relabel behaviors.
        if nv.relabel_behaviors
            assert(iscell(nv.relabel_behavriors));
            assert(size(nv.relabel_behaviors, 2) == 2);
            for i_relabel = 1:size(nv.relabel_behaviors, 1)
                sessions_raw{i_session} ...
                    = relabel_behavior(sessions_raw{i_session}, ...
                                       nv.relabel_behaviors{i_relabel,1}, ...
                                       nv.relabel_behaviors{i_relabel,2});
            end
        end

        % Try all behaviors. Will not examine baseline presently as Dayu's
        % group did not seem to look at this in their analyses.
        for i_bhv = 1:n_bhvs
            % Try current behavior.
            try 
                timeseries{i_subj,i_session,i_bhv} ...
                    = isolate_behavior(sessions_raw{i_session}.Lfold, ...
                                       sessions_raw{i_session}, ...
                                       bhv_list{i_bhv}, ...
                                       'min_length', nv.min_length, ...
                                       'bin_width', nv.bin_width, ...
                                       'step_size', nv.step_size);
            catch e
                % Legitimate exceptions, continue to next behavior.
                if ismember(e.identifier, ...
                            ["isolate_behavior:behavior_not_found", ...
                             "isolate_behavior:timeseries_of_insufficient_length", ...
                             "determine_n_bins:no_timepoints_left"])
                    
                    continue
                % Rethrow for debugging.
                else   
                    rethrow(e)
                end
            end
        end
    end
end

end
