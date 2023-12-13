function [timeseries_bhv, timeseries_base] = isolate_all_subjects_sessions_bhvs_alist(path_to_data, area_list, bhv_list, nv)

arguments
    path_to_data
    area_list
    bhv_list
    nv.remove_nan = true
    nv.nan_warning = true
    nv.zscore = false
    nv.smooth = false
    nv.relabel_behaviors = false
    nv.deriv = false
    nv.bin_width = 1
    nv.step_size = 1
    nv.min_length = 1
end

% Get indices of all subjects.
subj_ind = find(cellfun(@isempty, area_list) == 0);
if iscolumn(subj_ind), subj_ind = subj_ind'; end

% Get highest subject index, max number of sessions (for any subject) and
% max number of behaviors.
i_last_subj = subj_ind(end);
N_SESSIONS_MAX = 100;
n_bhvs = length(bhv_list);

% Extract portions of each session.
timeseries_bhv = cell(i_last_subj, N_SESSIONS_MAX, n_bhvs);
timeseries_base = cell(i_last_subj, N_SESSIONS_MAX, n_bhvs);
for i_subj = subj_ind
    % Get all sessions for current subject.
    sessions_raw = extract_single_subj(path_to_data, i_subj, ...
                                       'session', 'all', ...
                                       'combine_sessions', false, ...
                                       'region_names', true, ...
                                       'remove_nan', false, ...
                                       'nan_warning', true, ...
                                       'zscore', nv.zscore, ...
                                       'smooth', nv.smooth, ...
                                       'plot', false);

    % Iterate over all sessions for cu n_vars x n_timepoints dataset, compute varaince for each variable
% then rescale so that the variances across regions is distributed as
% Poi(mean_var).rrent subject.
    n_sessions = length(sessions_raw);
    for i_session = 1:n_sessions
        n_regions = size(sessions_raw{i_session}.Lfold, 1);
        % Set areas not on list to NaN.
        rmv_ind = setdiff(1:n_regions, area_list{i_subj});
        sessions_raw{i_session}.Lfold(rmv_ind,:) = nan;

        % Optionally relabel behaviors.
        if iscell(nv.relabel_behaviors)
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
                % Extract behavior.
                timeseries_bhv{i_subj,i_session,i_bhv} ...
                    = isolate_behavior(sessions_raw{i_session}.Lfold, ...
                                       sessions_raw{i_session}, ...
                                       bhv_list{i_bhv}, ...
                                       'deriv', nv.deriv, ...
                                       'min_length', nv.min_length, ...
                                       'bin_width', nv.bin_width, ...
                                       'step_size', nv.step_size);

                % If behavior sucessfully extracted, extract baseline as
                % well.
                timeseries_base{i_subj,i_session,i_bhv} ...
                    = isolate_baseline(sessions_raw{i_session}.Lfold, ...
                                       sessions_raw{i_session}, ...
                                       'deriv', nv.deriv, ...
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
