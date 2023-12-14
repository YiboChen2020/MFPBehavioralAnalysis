
clearvars
close all
clc


%% Add to path; run setup for CNMF-E

if true
    cnmfe_setup
end


%% Load data and preprocess

subj_sessions ...
    = {1, '190912', 8; ...  % 1
       6, '200311', 8; ...  % 2 
       6, '200312', 8; ...  % 3 
       8, '200527', 8; ...  % 4 
       8, '200611', 8; ...  % 5 
       8, '200615', 8; ...  % 6
       8, '200617', 8; ...  % 7
       12, '200801', 8; ... % 8 
       12, '200803', 8; ... % 9 
       17, '200708', 8; ... % 10 
       17, '200709', 7; ... % 11 
       17, '200719', 8; ... % 12 
       7,  '200527', 8; ... % 13 
       7,  '200530', 7; ... % 14 
       7,  '200604', 8; ... % 15
       7,  '200617', 8};    % 16

% Choose a subject of the above sessions to run.
selection_vec = 1:16;
subj_sessions = subj_sessions(selection_vec,:);
n_subj_sess = size(subj_sessions, 1);

% Set whether or not to save empirical results.
params.save_emp = false;

% Set whether to simulate data (this will automatically be saved).
params.save_sim = false;

% Set directory for saving results.
params.save_dir_emp = '';
params.save_dir_sim = '';

% Set whether to preprocess the entire session or individual interactions.
params.portion = 'interaction'; % {'wholesession' | 'interaction'}

for i_subj_sess = 1:n_subj_sess
    % Specify data location and subject to be analyzed. 
    params.path_to_data = '';
    
    % Set whether or not to z-score and/or smooth raw data after extraction.
    params.zscore_raw_extracted = true;
    params.smooth_raw_extracted = false;
    
    % Subject and session ID.
    params.i_subject = subj_sessions{i_subj_sess,1};
    params.session_id = subj_sessions{i_subj_sess,2};
    
    % Extract all sessions (or subset of sessions) data for one subject.
    [session_raw, ~] ...
        = extract_single_subj(params.path_to_data, params.i_subject, ...
                              'session', {params.session_id}, ...
                              'combine_sessions', false, ...
                              'bhvs_to_exclude', {}, ...
                              'region_names', true, ...
                              'zscore', params.zscore_raw_extracted, ...
                              'smooth', params.smooth_raw_extracted, ...
                              'plot', false);
    
    % Get list of behaviors that occurred across all sessions for current
    % subject. This will be used to set up a dictionary assigning a unique
    % number to each behavior (across all sessions for a single subject).
    all_subj_behaviors = get_all_behaviors(params.path_to_data, params.i_subject);
    
    % Number of regions retained.
    n_regions = session_raw{:}.n_regions_retained;
    
    
    %% Optionally relabel some behaviors
    
    % If true, this will cause certain human annotations to be replaced by
    % others. Note that this does not affect that all_subj_behaviors list,
    % which is used to generate number codes, as we want that list (and the
    % number codes, as well as the dicts linking bhvs to number codes) to
    % be as consistent across sessions/preprocessing as possible.
    params.relabeling.relabel = true;
    params.relabeling.old_new_labels ...
        = {'Attempted_Attack', 'Attack'};
    
    if params.relabeling.relabel 
        for i_relabel = 1:size(params.relabeling.old_new_labels, 1)
            session_raw = relabel_behavior(session_raw{:}, ...
                                           params.relabeling.old_new_labels{i_relabel,1}, ...
                                           params.relabeling.old_new_labels{i_relabel,2});
            session_raw = {session_raw};
        end
    end
    
    
    %% Calculate mean durations
    % This can inform the selection of bin size. The labels here should not be
    % used, as they will not be valid after any binning.
    
    % Get behavior labels.
    [~, exploratory.bhv_labels, exploratory.aux] ...
        = labeled_timeseries({session_raw{:}.Lfold}, session_raw, 'bin_width', 1, 'step_size', 1);
    
    % Calculate mean durations for each behavior in each session.
    exploratory.mean_durations ...
        = calc_mean_duration(exploratory.bhv_labels{:}, [], []);
    shortest_mean_duration = min(exploratory.mean_durations);
    suggested_max_bin_size = floor(0.5 * shortest_mean_duration);
    
     
    %% Set which types of preprocessing to apply
    
    params.data_type = 'fmpp';
    
     
    %% Optional deconvolution
    % Need to run setup for CNMF-E.
    
    if strcmp(params.data_type, 'deconv')
    
    % Set parameters related to kernel.
    params.deconv.kernel_model = 'ar1'; 
    params.deconv.kernel_params = []; % Corresponds to 'pars', default = []
    params.deconv.method = 'constrained'; 
    params.deconv.extra_params = []; % Only if method = 'mcmc'
    params.deconv.window = 5; % Width of kernel
    
    % Other parameters.
    params.deconv.noise_std = []; % Default is [], in which case this will be estimated by function
    params.deconv.baseline = 0; % Default = 0
    params.deconv.lambda = 0; % Default = 0
    params.deconv.optimize_b = true;
    params.deconv.optimize_pars = true;
    params.deconv.optimize_smin = true;
    
    % As of 5/19/22, any parameters not specified should be assumed to
    % match the defaults at the time the repo was cloned on 5/18/22.
    
    % Preallocate.
    neural_data = NaN(size(session_raw{:}.Lfold));
    
    % Deconvolve.
    parfor i_region = 1:n_regions
        [~, neural_data(i_region,:), ~] ...
            = deconvolveCa(session_raw{:}.Lfold(i_region,:), ...
                           params.deconv.kernel_model, params.deconv.kernel_params, ...
                           params.deconv.method, params.deconv.extra_params, ...
                           'window', params.deconv.window, ...
                           'sn', params.deconv.noise_std, ...
                           'b', params.deconv.baseline, ...
                           'lambda', params.deconv.lambda, ...
                           'optimize_b', params.deconv.optimize_b, ...
                           'optimize_pars', params.deconv.optimize_pars, ...
                           'optimize_smin', params.deconv.optimize_smin);       
    end
    
    end
    
    
    %% Option to visualize deconvolved data (if deconvolution performed).
    
    params.plotting.deconv.plot = true;
    params.plotting.deconv.window = [5500 65000];
    
    if strcmp(params.data_type, 'deconv') && params.plotting.deconv.plot
        figure
        i_subplot = 1;
        
        for i_region = 1:n_regions
            % Prepare window.
            if strcmp(params.plotting.deconv.window, 'entire')
                plot_window = [1 size(neural_data, 2)];
            else
                plot_window = params.plotting.deconv.window;
            end
            
            % Plot current session.
            subplot(n_regions, 1, i_subplot)
            hold on
            stem(neural_data(i_region,:), 'Marker', 'none')
            xlim(plot_window)
    %         xlabel(sprintf('Time at 25 Hz from %dth to %dth timepoint', plot_window(1), plot_window(2)))
    %         title(sprintf('Session %s', selected_session_ids{i_session}))
            clear plot_window
            
            i_subplot = i_subplot + 1;
        end  
    end
    
    
    %% Optional peak detection + convolution to create filtered Markov point process
    
    if strcmp(params.data_type, 'fmpp')
    
    % Set threshold (as percentage of max activity for each region) below which
    % all activity will be set to zero).
    params.fmpp.threshold_percent = -Inf;
    params.fmpp.t_0 = 0;
    params.fmpp.t = -2:2;
    params.fmpp.alpha = 1;
    params.fmpp.filter ...
        = create_filter(params.fmpp.t_0, params.fmpp.t, params.fmpp.alpha, 'renorm', false); % [0.1 0.2 0.4 0 0] / 0.7;
    params.fmpp.zscore = false;
    
    % Get filtered MPP.
    neural_data ...
        = marked_point_process(session_raw{:}.Lfold, ...
                               'min_threshold', params.fmpp.threshold_percent, ...
                               'zscore', params.fmpp.zscore, ...
                               'convolve', params.fmpp.filter, ...
                               'plot', false, ...
                               'figure', true, ...
                               'window', [3000 7500]);
    end
    
    
    %% If using raw data, assign to neural data
    
    if strcmp(params.data_type, 'raw'), neural_data = session_raw{:}.Lfold; end
    
    
    %% Extract, bin, and label annotated periods during interaction

    % Create dictionaries mapping behaviors to number codes and back.
    params.dict = containers.Map(all_subj_behaviors, 1:length(all_subj_behaviors));
    params.r_dict = containers.Map(1:length(all_subj_behaviors), all_subj_behaviors);

    if strcmp(params.portion, 'wholesession')
        % Generate labels on all timepoints at original 25 Hz.
        n_timepoints_total = size(neural_data, 2);
        labels = nan(n_timepoints_total, 1);
        for i_epoch = 1:length(session_raw{:}.behaviors)
            i_start = session_raw{:}.Fstart(i_epoch);
            i_end = session_raw{:}.Fstop(i_epoch);
        
            labels(i_start:i_end) = params.dict(session_raw{:}.behaviors{i_epoch});
        end
        
        % Sometimes, there may be a mismatch between last Fstop and the
        % number of timepoints in the neural data. If Fstop exceeds the
        % number of timepoints, this should result in an error here.
        % However, if the last Fstop is a number smaller than the number of
        % timepoints, this will result in labels with nan, which will
        % propagate forward past binning (an exception will be thrown
        % during computation of mean activity, which uses the labels, due
        % to nan values in the labels, which originate from the
        % preallocated nan values for the var labels here). Thus, we will
        % trim off any excess timepoints from the end of the neural data.
        if any(isnan(labels))
            % First check for NaNs.
            if sum(isnan(labels)) == 1
                warning("There is one NaN value in the var labels.")
            else
                warning("There are %d NaN values in the var labels.", sum(isnan(labels)))
            end
        
            % Remove any trailing timepoints past those indexed by Fstop.
            neural_data(:,i_end+1:end) = [];
            labels(i_end+1:end) = [];
        
            % Check if any NaNs remain. If so, warn; optionally, throw error.
            if any(isnan(labels))
                if true
                    error("NaNs in data not due to data beyond final Fstop.")
                else
                    warning("NaNs in data not due to data beyond final Fstop.")
                end
            end
        end
    
        % Optionally remove data at the beginning of the array that is all
        % zeros. This is done before binning so that the transition indices are
        % preserved etc.
        params.remove_leading_zeros = false;
        if params.remove_leading_zeros
            % Use while loop to ensure only indices of all zero bins in the
            % front portion are identified.
            remove_ind = [];
            i_timepoint = 1;
            while true
                if sum(neural_data(:,i_timepoint)) == 0
                    remove_ind = [remove_ind; i_timepoint];
                    i_timepoint = i_timepoint + 1;
                else
                    break
                end
            end
    
            % Remove data.
            params.leading_zeros_remove_ind = remove_ind;
            neural_data(:,remove_ind) = [];
            labels(remove_ind) = [];
        end

        % Place in cell array to match the interaction format.
        neural_data_to_bin = {neural_data};
        labels = {labels};
        session_raw{:}.n_interactions = 1;
        
    elseif strcmp(params.portion, 'interaction')
        % Set possible identities of second subject. Also set whether or not to
        % include the Intro and Rmv epochs.
        params.allowed_second_subjects = {'M', 'F', 'Toy'};
        params.include_intro_rmv = false;
        
        % Extract only portions of timeseries that featured two subjects
        % (optionally counting a toy as a subject).
        [neural_data_to_bin, labels, second_subj, ~, ~] ...
            = isolate_interactions_3(neural_data, session_raw{:}, params.dict, ...
                                     'allowed_second_subjects', params.allowed_second_subjects, ...
                                     'include_intro_rmv', params.include_intro_rmv);

        % Save identities of all second subjects from current session.
        session_raw{:}.all_second_subjects = second_subj;
        session_raw{:}.n_interactions = length(session_raw{:}.all_second_subjects);

        % Note that prior to binning in the case of preprocessing the full
        % timeseries (hmm_preprocessing_pipeline_2) we first remove any
        % leading zeros, but this is not necessary here, as we are
        % isolating the portion of the timeseries that contains an
        % interaction between two subjects.
    else
        error("Unrecognized value for params.portion.")
    end
    
    % Bin, resample, and save for each interaction (or whole session).
    for i_interaction = 1:session_raw{:}.n_interactions
        % Get behaviors that were kept.
        params.bhv_codes = unique(labels{i_interaction});
        params.bhv_names = cell(length(params.bhv_codes), 1);
        for i_bhv = 1:length(params.bhv_codes)
            params.bhv_names{i_bhv} = params.r_dict(params.bhv_codes(i_bhv));
        end

        if strcmp(params.portion, 'interaction')
            % If was interaction only, store identity of second subject.
            params.second_subject = second_subj{i_interaction}; % Same number of letters as M and F
            assert(ismember(params.second_subject, params.allowed_second_subjects));
            params.second_subject = strrep(second_subj{i_interaction}, 'Toy', 'T'); % Same number of letters as M and F
            
            % Print progress to console.
            disp(" ")
            disp("-------------------------------------------------------------")
            if i_interaction == 1
                disp("-------------------------------------------------------------")
            end
            fprintf('\nSubject: %d, Session: %s, 2nd %s. Interaction %d. Number of behaviors: %d.\n', ...
                    params.i_subject, ...
                    params.session_id, ...
                    params.second_subject, ...
                    i_interaction, ...
                    length(params.bhv_codes))
        else
            % Print progress to console.
            fprintf('\nSubject: %d, Session: %s, Whole session. Number of behaviors: %d.\n', ...
                    params.i_subject, ...
                    params.session_id, ...
                    length(params.bhv_codes))

        end

        
        % Set bin width and step size.
        params.bin_width = subj_sessions{i_subj_sess, 3};
        params.step_size = params.bin_width;
        
        % Bin data. If the attempted bin size is too large, attempt smaller
        % bin sizes until no bin contains more than one transition.
        while true
            try
                % Try current bin size.
                [neural_data_binned, bin_labels, transitions, mixed_transitions] ...
                    = bin_with_transitions(neural_data_to_bin{i_interaction}, ...
                                           labels{i_interaction}', ...
                                           params.bin_width, ...
                                           params.step_size);
                break
            catch e
                % If current bin size results in multiple transitions
                % within one bin, warn, decrease bin size by 1 and try
                % again.
                if strcmp(e.identifier, 'bin_with_transitions:multiple_transitions_within_bin')
                    warning(['Bin size of %d resulted in more than one transition' ...
                             ' within a bin. Trying bin_width = %d.'], ...
                             params.bin_width, params.bin_width - 1);
                    params.bin_width = params.bin_width - 1;
                    params.step_size = params.bin_width;
                    if params.bin_width == 0
                        error(['Something went wrong during binning ' ...
                               '(a bin of size 1 still contained multiple ' ...
                               'transitions.'])
                    end
                else
                    rethrow(e);
                end
            end
        end


        %% Optional poisson resampling
        
        % Set to false/0 to skip resampling. Else, set to lambda, mean/var
        % of the Poisson distribution to which to resample.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        params.resampling.resample = 5;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if params.resampling.resample
        
        % Set other resampling parameters.
        params.resampling.resample_regions_together = false;
        params.resampling.seed = 'default';
        
        if params.resampling.resample_regions_together
            % Resample all regions together.
            neural_data_binned = resample2poisson(neural_data_binned, ...
                                                  params.resampling.resample, ...
                                                  params.resampling.seed);
        else
            % Resample for each region separately.
            for i_region = 1:n_regions
                neural_data_binned(i_region,:) ...
                    = resample2poisson(neural_data_binned(i_region,:), ...
                                       params.resampling.resample, ...
                                       params.resampling.seed);
            end
        end
        
        end


        %% Optionally divide into smaller subsections by cutting at Other
        % See choose_opt_subdivision in the proprocessing-functions
        % directory about automated segmentation.
        
        % Try a range of number of subsections and print the results of
        % each to the console ---------------------------------------------

        % Set the values for n_subsectsions to try. SET TO 1 to skip
        % subdivision.
        if ~strcmp(params.portion, 'wholesession')
            n_subsections_vals_to_try = 1:5;
        else
            n_subsections_vals_to_try = 1;
        end

        % Try each number of subsections (each value for n_subsections).
        subsections_trial = cell(length(n_subsections_vals_to_try), 1);
        other_trial = cell(length(n_subsections_vals_to_try), 1);
        for i_n_subsections = 1:length(n_subsections_vals_to_try)
            % Try current value for n_subsections.
            try
                [subsections_trial{i_n_subsections}, other_trial{i_n_subsections}] ...
                    = subdivide_session(neural_data_binned, ...
                                        bin_labels, ...
                                        params.dict('Other'), ...
                                        n_subsections_vals_to_try(i_n_subsections));
            catch e
                if strcmp(e.identifier, 'subdivide_session:n_subsections_gr_than_n_other')
                    break
                else
                    rethrow(e);
                end
            end

            % Print results to console.
            fprintf('\nn_subsections = %d \n', n_subsections_vals_to_try(i_n_subsections))
            for i_subsection = 1:n_subsections_vals_to_try(i_n_subsections)
                % Get names of behaviors in current subsection.
                labels_list = unique(subsections_trial{i_n_subsections}{i_subsection}.labels);
                bhvs_in_subsection = cell(length(labels_list), 1);
                for i_label = 1:length(labels_list)
                    bhvs_in_subsection{i_label} = params.r_dict(labels_list(i_label));
                end

                % Print.
                fprintf('duration = %d. richness = %f. diversity = %d of %d non-Other bhvs: %s out of %s. \n', ...
                        length(subsections_trial{i_n_subsections}{i_subsection}.labels), ...
                        subsections_trial{i_n_subsections}{i_subsection}.richness, ...
                        subsections_trial{i_n_subsections}{i_subsection}.n_bhvs, ...
                        length(params.bhv_codes) - 1, ...
                        strjoin(setdiff(bhvs_in_subsection, 'Other'), ', '), ...
                        strjoin(setdiff(params.bhv_names, 'Other'), ', '));

            end
        end
        
        % Stop to allow user to check subdivisions before setting a final
        % value and saving.
        if ~all(n_subsections_vals_to_try == 1)
            params.n_subsections = input('\nEnter a final value for n_subsections: ');
            params.subsections_to_keep = input('\nEnter indices of the subsections to keep: ');
            assert(all(ismember(params.subsections_to_keep, 1:params.n_subsections)), ...
                   sprintf(['One or more of the indices of the subsections to keep is out of the' ...
                            'valid range from 1 to %d'], params.n_subsections));
        else
            params.n_subsections = 1;
            params.subsections_to_keep = 1;
        end

        % Subdivide.
        [subsections, other] ...
            = subdivide_session(neural_data_binned, ...
                                bin_labels, ...
                                params.dict('Other'), ...
                                params.n_subsections);

        % Store info about subdivisions. 
        params.other_subdivision = other;

        for i_subsection = 1:params.n_subsections
            % Skip if not one of the subsections to be kept.
            if ~ismember(i_subsection, params.subsections_to_keep), continue; end

            % Get first subsection and rename it.
            neural_data_binned = subsections{i_subsection}.timeseries;
            bin_labels = subsections{i_subsection}.labels;
            subsection_info.head = subsections{i_subsection}.head;
            subsection_info.tail = subsections{i_subsection}.tail;

            % Get behaviors in subsection.
            labels_retained = unique(bin_labels);
            subsection_info.bhvs = cell(length(labels_retained), 1);
            for i_label = 1:length(labels_retained)
                subsection_info.bhvs{i_label} = params.r_dict(labels_retained(i_label));
            end
        
            %% Optionally visualize neural data and labels
            
            if false
                figure
                imagesc(neural_data_binned)
                colorbar
                title(sprintf('Subject: %d, Session: %s, Interaction %d, 2nd %s', ...
                              params.i_subject, ...
                              params.session_id, ...
                              i_interaction, ...
                              params.second_subject))
                
                figure
                imagesc(bin_labels')
                colorbar
                title(sprintf('Subject: %d, Session: %s, Interaction %d, 2nd %s', ...
                              params.i_subject, ...
                              params.session_id, ...
                              i_interaction, ...
                              params.second_subject))
            end
            
            
            %% Calculate Markov statistics after any binning/resampling
            
            hmm_stats.mean_duration = calc_mean_duration(bin_labels, [], []);
            hmm_stats.mean_activity = calc_mean_activity(neural_data_binned, bin_labels, [], []);
            hmm_stats.trans_mat = get_transition_mat(bin_labels, [], []);
            hmm_stats.bhvs = params.bhv_names;
            
            
            %% Save empirical data

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            run_id = 'B';            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Create string for either wholesession or interaction to place
            % in filename.
            if strcmp(params.portion, 'wholesession')
                portion_name = 'wholesession';
            elseif strcmp(params.portion, 'interaction')
                portion_name = sprintf("2nd%s%d", params.second_subject, i_interaction);
            else
                error("Unrecognized value for params.portion.")
            end

            % If n_subsections == 1, this means that we used the entire
            % interaction (nothing is cut out either). In that case, we
            % want to create a special designiation, as subsection1 would
            % cause the file to be overwritten when subdividing later.
            if params.n_subsections == 1
                subsection_name = 0;
            else
                subsection_name = i_subsection;
            end

            % Build save filename.
            save_filename = sprintf("%d-%s-%s-zscore%d-bin%d-relabeled%d-resampled%d-%s-subsection%d-%s", ...
                             params.i_subject, strjoin(session_raw{:}.session_id, '.'), ...
                             params.data_type, ...
                             params.zscore_raw_extracted, ...
                             params.bin_width, ...
                             params.relabeling.relabel, ...
                             params.resampling.resample, ...
                             portion_name, ...
                             subsection_name, ...
                             run_id);
            
            if params.save_emp
                % Switch to save directory then back.
                curr_dir = pwd;
                cd(params.save_dir_emp)
                if params.resampling.resample
                    cd('resampled')
                else
                    cd('not-resampled')
                end
            
                % Replace Lfold with its shape to save space.
                session_raw{:}.Lfold = size(session_raw{:}.Lfold);
    
                save(save_filename, ...
                     "neural_data_binned", ...
                     "bin_labels", ...
                     "params", ...
                     "session_raw", ...
                     "hmm_stats", ...
                     "exploratory", ...
                     "subsection_info", ...
                     "other");
                
                cd(curr_dir);
            end

            %% Optionally generate and save simulation data
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sim_id = "A";
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
            % Skip if there's only one behavior.
            if size(hmm_stats.mean_activity, 1) == 1, continue; end

            % Generate simulated data.
            [neural_data_binned, bin_labels, degenerate, seed] ...
                = gen_poiss_data_trans_mat(hmm_stats.mean_activity, ...
                                           hmm_stats.trans_mat, ...
                                           2000, ...
                                           'starting_bhv_idx', [], ...
                                           'remove_degenerate', true, ...
                                           'seed', 'shuffle');
            degenerate = find(degenerate);

            % Switch to save directory then back.
            curr_dir = pwd;
            cd(params.save_dir_sim)
            
            % If Lfold was not replaced with its size above, do so
            % here.
            if size(session_raw{:}.Lfold, 1) > 1
                session_raw{:}.Lfold = size(session_raw{:}.Lfold);
            end

            % Save.
            if params.save_sim
                save(save_filename + "-" + sim_id, ...
                     "neural_data_binned", ...
                     "bin_labels", ...
                     "degenerate", ...
                     "seed", ...
                     "params", ...
                     "session_raw", ...
                     "hmm_stats", ...
                     "exploratory", ...
                     "subsection_info", ...
                     "other");

                cd(curr_dir);
            end


        end
    end
end

