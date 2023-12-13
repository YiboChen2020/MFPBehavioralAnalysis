% This script is used to generate the HSMM annotation.

% Author: yibochen, 11.29.2023

clc
clearvars
close all

subj_sessions ...
    = {1, '190912', 8;
       1, '190922', 8;
       1, '190925', 8;
       1, '191010', 8;
       6, '200311', 8; %5
       6, '200312', 8;
       6, '200313', 8;
       6, '200314', 8;
       6, '200315', 8;
       7, '200527', 8; %10
       7, '200530', 7;
       7, '200604', 8;
       7, '200617', 8;
       8, '200527', 8;
       8, '200611', 8; %15
       8, '200615', 8;
       8, '200617', 8;
       12, '200801', 8;
       12, '200803', 8;
       17, '200708', 8; % 20
       17, '200709', 7;
       17, '200719', 8;
       };

selection_vec = [1,5,6];
subj_sessions = subj_sessions(selection_vec,:);
n_subj_sess = size(subj_sessions, 1);
Interaction_Type = 'Female';% Female, Male, Wholesession

pre_path = '/home/jmc/Dropbox (NYU Langone Health)/bhv-mfp/preprocessing/results/preprocessed-data/poisson-hmm/emp/resampled/';
if strcmp(Interaction_Type,'Wholesession')
    HSMM_path = '/home/jmc/Dropbox (NYU Langone Health)/py-hmm-2/results/2023-06-14-A/';
else
    HSMM_path = '/home/jmc/Dropbox (NYU Langone Health)/py-hmm-2/results/2023-06-07-C/';
end

for i_subj_sess = 1:n_subj_sess
    Animal_ID = num2str(subj_sessions{i_subj_sess,1});
    Session_ID = subj_sessions{i_subj_sess,2};
    bin = subj_sessions{i_subj_sess,3};
    
    if strcmp(Interaction_Type,'Wholesession')
        % load HSMM data
        cur_dir = pwd;
        cd(HSMM_path)
        file = dir([Animal_ID,'-',Session_ID,'*zscore0*.mat']);
        load([HSMM_path,file.name])
        cd(cur_dir)
        % load preprocess data
        cur_dir = pwd;
        cd(pre_path)
        file = dir([Animal_ID,'-',Session_ID,'*zscore0*wholesession*.mat']);
        load([pre_path,file.name])
        cd(cur_dir)
        % Annotation start and end time point
        RawS = session_raw{1,1};
        start_idx = 1;
        stop_idx = RawS.Lfold(2);
    elseif strcmp(Interaction_Type,'Male')
        % load HSMM data
        cur_dir = pwd;
        cd(HSMM_path)
        file = dir([Animal_ID,'-',Session_ID,'*zscore0*2ndM*subsection0*B*.mat']);
        load([HSMM_path,file.name])
        cd(cur_dir)
        % load preprocess data
        cur_dir = pwd;
        cd(pre_path)
        file = dir([Animal_ID,'-',Session_ID,'*zscore0*2ndM*subsection0*B*.mat']);
        load([pre_path,file.name])
        cd(cur_dir)
        RawS = session_raw{1,1};
        idx_intro = find(strcmp(RawS.behaviors,{'Intro_M'}));
        idx_rmv = find(strcmp(RawS.behaviors,{'Rmv_M'}));
        start_idx = idx_intro(end)+1;
        stop_idx = idx_rmv(1)-1;
    else % Female
        % load HSMM data
        cur_dir = pwd;
        cd(HSMM_path)
        file = dir([Animal_ID,'-',Session_ID,'*zscore0*2ndF*subsection0*B*.mat']);
        load([HSMM_path,file.name])
        cd(cur_dir)
        % load preprocess data
        cur_dir = pwd;
        cd(pre_path)
        file = dir([Animal_ID,'-',Session_ID,'*zscore0*2ndF*subsection0*B*.mat']);
        load([pre_path,file.name])
        cd(cur_dir)
        RawS = session_raw{1,1};
        idx_intro = find(strcmp(RawS.behaviors,{'Intro_F'}));
        idx_rmv = find(strcmp(RawS.behaviors,{'Rmv_F'}));
        start_idx = idx_intro(end)+1;
        stop_idx = idx_rmv(1)-1;
    end
    
    true_label = bin_labels;
    true_data= true_label;
    hmm_data = hidden_states + 1;
    
   
    if strcmp(Interaction_Type,'Female') == 1
        bhvrs = keys(params.dict);
        bhv_names = cell(1,6);
        for i = 1:length(bhvrs)
            switch bhvrs{i}
                case 'Investigate'
                    idx_i = find(true_label == params.dict('Investigate'));
                    true_data(idx_i) = 5;
                    if ~isempty(idx_i)
                        bhv_names{5} = 'Investigate';
                    end
                case 'Attempted_Mount'
                    idx = find(true_label == params.dict('Attempted_Mount'));
                    true_data(idx) = 1;
                    if ~isempty(idx)
                        bhv_names{1} = 'Attempted_Mount';
                    end 
                case 'Mount'
                    idx = find(true_label == params.dict('Mount'));
                    true_data(idx) = 2;
                    if ~isempty(idx)
                        bhv_names{2} = 'Mount';
                    end 
                case 'Thrust'
                    idx = find(true_label == params.dict('Thrust'));
                    true_data(idx) = 3;
                    if ~isempty(idx)
                        bhv_names{3} = 'Thrust';
                    end 
                case 'Ejaculate'
                    idx = find(true_label == params.dict('Ejaculate'));
                    true_data(idx) = 4;
                    if ~isempty(idx)
                        bhv_names{4} = 'Ejaculate';
                    end 
                otherwise
                    true_data(true_label == params.dict(bhvrs{i})) = 6;
                    bhv_names{6} = 'Other';
            end
        end
    elseif strcmp(Interaction_Type,'Male') == 1
        bhvrs = keys(params.dict);
        bhv_names = cell(1,3);
        for i = 1:length(bhvrs)
            switch bhvrs{i}
                case 'Attack'
                    idx = find(true_label == params.dict('Attack'));
                    true_data(idx) = 1;
                    if ~isempty(idx)
                        bhv_names{1} = 'Attack';
                    end
                case 'Investigate'
                    idx_i = find(true_label == params.dict('Investigate'));
                    true_data(idx_i) = 2;
                    if ~isempty(idx_i)
                        bhv_names{2} = 'Investigate';
                    end
                otherwise
                    true_data(true_label == params.dict(bhvrs{i})) = 3;
                    bhv_names{3} = 'Other';
            end
        end
    else
        bhvrs = keys(params.dict);
        bhv_names = cell(1,9);
        for i = 1:length(bhvrs)
            switch bhvrs{i}
                case 'Attack'
                    idx = find(true_label == params.dict('Attack'));
                    true_data(idx) = 2;
                    if ~isempty(idx)
                        bhv_names{2} = 'Attack';
                    end
                case 'Investigate'
                    idx_i = find(true_label == params.dict('Investigate'));
                    if isKey(params.dict, 'Intro_M')
                        idx1 = find(true_label == params.dict('Intro_M'));
                        idx2 = find(true_label == params.dict('Rmv_M'));
                        idx_m = [idx1(end), idx2(1)];
                        idx_MI = idx_i( idx_i > idx_m(1));
                        idx_MI = idx_MI(idx_MI < idx_m(end));
                        true_data(idx_MI) = 1;
                        idx_i = setdiff(idx_i,idx_MI);
                        if ~isempty(idx_MI)
                            bhv_names{1} = 'Male_Investigate';
                        end
                    end
                    if isKey(params.dict, 'Intro_F')
                        idx1 = find(true_label == params.dict('Intro_F'));
                        idx2 = find(true_label == params.dict('Rmv_F'));
                        idx_f = [idx1(end), idx2(1)];
                        idx_FI = idx_i( idx_i > idx_f(1));
                        idx_FI = idx_FI(idx_FI < idx_f(end));
                        true_data(idx_FI) = 3;
                        idx_i = setdiff(idx_i,idx_FI);
                        if ~isempty(idx_FI)
                            bhv_names{3} = 'Female_Investigate';
                        end
                    end
                    if isKey(params.dict, 'Intro_Toy')
                        idx1 = find(true_label == params.dict('Intro_Toy'));
                        if ~isempty(idx1)
                            idx2 = find(true_label == params.dict('Rmv_Toy'));
                            idx_t = [idx1(end), idx2(1)];
                            idx_TI = idx_i( idx_i > idx_t(1)); 
                            idx_TI = idx_TI(idx_TI < idx_t(end));
                            true_data(idx_TI) = 8;
                            idx_i = setdiff(idx_i,idx_TI);
                            if ~isempty(idx_FI)
                                bhv_names{8} = 'Toy_Investigate';
                            end
                        end
                    end
                    true_data(idx_i) = 9;% there is a pup invest in one session, set it as other
                case 'Attempted_Mount'
                    idx = find(true_label == params.dict('Attempted_Mount'));
                    true_data(idx) = 4;
                    if ~isempty(idx)
                        bhv_names{4} = 'Attempted_Mount';
                    end 
                case 'Mount'
                    idx = find(true_label == params.dict('Mount'));
                    true_data(idx) = 5;
                    if ~isempty(idx)
                        bhv_names{5} = 'Mount';
                    end 
                case 'Thrust'
                    idx = find(true_label == params.dict('Thrust'));
                    true_data(idx) = 6;
                    if ~isempty(idx)
                        bhv_names{6} = 'Thrust';
                    end 
                case 'Ejaculate'
                    idx = find(true_label == params.dict('Ejaculate'));
                    true_data(idx) = 7;
                    if ~isempty(idx)
                        bhv_names{7} = 'Ejaculate';
                    end 
                otherwise
                    true_data(true_label == params.dict(bhvrs{i})) = 9;
                    bhv_names{9} = 'Other';
            end
        end
    end
    unique(true_data);
    bhv_names = bhv_names(~cellfun(@isempty, bhv_names));
    
    % Get all_state assignment for each state.
    SKIP_OTHER_PURITY = false;
    [reorder_ind, ...
     row_membership, ...
     row_membership_unsorted, ...
     max_overlaps, ...
     occupancy] ...
        = all_state_match(true_data, hmm_data);
    hidden_states ...
        = reorder_hidden_states(hmm_data, reorder_ind);
    hidden_states_list ...
        = unique(hidden_states);
    occupancy ...
            = get_occupancy_mat(true_data, ...  
                            hidden_states);
    [purity, purity_by_bhv] ...
        = compute_all_state_purity(row_membership, ...
                                   occupancy, ...
                                   SKIP_OTHER_PURITY);
    %% save the HSMM annotation
    annotation = hidden_states;
    annotation = reshape(repmat(annotation,bin,1),[],1);
    Confusion_Matrix = occupancy;
    Behavior_Names = bhv_names;
    HSMM_Annotation = annotation;
    cur_dir = pwd;
    if strcmp(Interaction_Type,'Wholesession')
        cd("/home/jmc/Dropbox (NYU Langone Health)/bhv-mfp/hmm-results/HSMMannotation/resultss/whole session/")
    else
        cd("/home/jmc/Dropbox (NYU Langone Health)/bhv-mfp/hmm-results/HSMMannotation/resultss/without_baseline/")
    end
%     save([Animal_ID,'-',Session_ID,'-',Interaction_Type,'.mat'], ...
%         "HSMM_Annotation","stop_idx","start_idx","Confusion_Matrix", ...
%         "Animal_ID","Session_ID","Interaction_Type","Behavior_Names")
    cd(cur_dir)
    
    %%
    figure;
    h1 = heatmap(Confusion_Matrix);
    h1.YDisplayLabels = strrep(Behavior_Names,'_',' ');
    title([Animal_ID '-' Session_ID]);
end