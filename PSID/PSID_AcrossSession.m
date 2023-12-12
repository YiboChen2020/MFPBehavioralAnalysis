%% Notes:
% This is a script calling functions to use PSID to find the connection
% between the neural signal and behavior for different session of one subject.
% Using one session's male interaction neural and tracking data as train
% data, test on another session's male interaction neural data.Plot an
% imagesc.

% Author: yibochen, 10/02/23
clearvars
close all
clc

% Set path
scriptFullPath = mfilename('fullpath');
[scriptPath,~,~] = fileparts(scriptFullPath);
cd(scriptPath)
% Add path
allPaths = genpath("./source");
addpath(allPaths);

%% Params
sessions = {
    6, '200311';
    6, '200312';
    6, '200313';
    6, '200314';
    6, '200315';
    };
select_session = 1:5;
select_feature = [10,20,21,24,25,27,29];

%% initialize
male_neural_data = cell(length(select_session),1);
male_tracking_data = cell(length(select_session),1);
corr_matrix_neural = zeros(length(select_session));
corr_matrix_behavior = zeros(length(select_session));

for ii = 1:length(select_session)
    %% load data
    load(['../Data/MFP_raw/MFP-ERa-GC', num2str(sessions{select_session(ii),1}), '.mat'])
    RawS = Raw.(['r',sessions{select_session(ii),2}]);
    neural_data = zscore(RawS.Lfold,0,2);
    behaviors = RawS.behaviors;
    FrameS = RawS.Fstart;
    FrameE = RawS.Fstop;
    idx_IntroM = find(strcmp(behaviors,{'Intro_M'}));
    idx_RmvM = find(strcmp(behaviors,{'Rmv_M'}));
    male_idx = FrameS(idx_IntroM):FrameE(idx_RmvM);
    male_neural_data{ii} = neural_data(:,male_idx);


    %% load tracking data
    curr_dir = pwd;
    cd('../Data/tracking_preprocessed');
    all_files = dir();
    tracking_data_path = {};
    for i = 1:length(all_files)
        file_name = all_files(i).name;
        if (contains(file_name, strcat(num2str(sessions{select_session(ii),1}),'-',sessions{select_session(ii),2}))) && (contains(file_name, 'bin1'))
            tracking_data_path = file_name;
        end
    end
    load(tracking_data_path);
    tracking_data = F_binned(select_feature,:);
    tracking_data(isnan(tracking_data)) = 0;
    tracking_data = zscore(tracking_data,0,2);
    % normalize to [0.1, 0.9]
    tracking_min_per_row = min(tracking_data, [], 2);
    tracking_max_per_row = max(tracking_data, [], 2);
    tracking_norm = (tracking_data - 0.99*tracking_min_per_row) ./ (1.01*tracking_max_per_row - tracking_min_per_row);
    tracking_data = 0.8 * tracking_norm + 0.1;
    male_tracking_data{ii} = tracking_data(:,male_idx);
    cd(curr_dir)
end

%% PSID
% params
nx = 13;
n1 = 6;
i_para = 10;
k_fold = 10;
%k-fold validation for get best nx and n1
for i = 1:length(select_session)
    for j = 1:length(select_session)
        %k-fold validation for i-i prediction
        if i == j 
            scaler_temp1 = 0;
            scaler_temp2 = 0;
            for k  = 1:k_fold
                % divide train and test data segment
                testlength = floor(length(male_neural_data{i})/k_fold);
                test_idx = 1 + (k-1)*testlength: k*testlength;
                train_idx = setdiff(1:length(male_neural_data{i}),test_idx);
                % apply PSID
                Sys_male = PSID(male_neural_data{i}(:,train_idx)',male_tracking_data{i}(:,train_idx)',nx,n1,i_para);
                [zTestPred, yTestPred, xTestPred] = PSIDPredict(Sys_male, male_neural_data{i}(:,test_idx)');
                zTestPred = zTestPred';
                yTestPred = yTestPred';
                xTestPred = xTestPred';
                %calculate behavior prediction correlation
                for ii = 1:size(zTestPred,1)
                    scaler_temp1 = scaler_temp1 + corr(zTestPred(ii,:)',male_tracking_data{i}(ii,test_idx)');
                end
                for ij = 1:size(yTestPred,1)
                    scaler_temp2 = scaler_temp2 + corr(yTestPred(ij,:)',male_neural_data{i}(ij,test_idx)');
                end
                %self prediction plot
                if k == 1
                    figure;
                    hold on;
                    ymulti = 10;
                    for ii=1:size(zTestPred,1)
                        zTestPred(ii,:) = zscore(zTestPred(ii,:),0,2);
                        plot(ii*ymulti+zTestPred(ii,:),'linewidth',0.5,'Color','r');
                        plot(ii*ymulti+zscore(male_tracking_data{j}(ii,test_idx),0,2),'linewidth',0.5,'Color','k');
                        text(min(xlim),ii*ymulti,strrep(feat_names{select_feature(ii)},'_',' '),'HorizontalAlignment','right');
                        temp_corr = corr(zTestPred(ii,:)',zscore(male_tracking_data{j}(ii,test_idx)',0,1));
                        text(max(xlim),ii*ymulti,num2str(temp_corr),'HorizontalAlignment','right');
                    end
                    xlabel('time');
                    title("self prediction");
                    yticks([]);
                    hold off;
                end
            end
            corr_matrix_behavior(i,j) = scaler_temp1/k/ii;
            corr_matrix_neural(i,j) = scaler_temp2/k/ij;
            continue
        end
        % apply PSID
        Sys_male = PSID(male_neural_data{i}',male_tracking_data{i}',nx,n1,i_para);
        [zTestPred, yTestPred, xTestPred] = PSIDPredict(Sys_male, male_neural_data{j}');
        zTestPred = zTestPred';
        yTestPred = yTestPred';
        xTestPred = xTestPred';
        %calculate behavior prediction correlation
        temp = 0;
        for k = 1:size(yTestPred,1)
            temp = temp + corr(yTestPred(k,:)',male_neural_data{j}(k,:)');
        end
        corr_matrix_neural(i,j) = temp/k;
        temp = 0;
        for k = 1:size(zTestPred,1)
            temp = temp + corr(zTestPred(k,:)',male_tracking_data{j}(k,:)');
        end
        corr_matrix_behavior(i,j) = temp/k;

        %make behavior label
        load(['../Data/MFP_raw/MFP-ERa-GC', num2str(sessions{select_session(j),1}), '.mat'])
        RawS = Raw.(['r',sessions{select_session(j),2}]);
        behaviors = RawS.behaviors;
        FrameS = RawS.Fstart;
        FrameE = RawS.Fstop;
        idx_IntroM = find(strcmp(behaviors,{'Intro_M'}));
        idx_RmvM = find(strcmp(behaviors,{'Rmv_M'}));
        idx_attack = find(strcmp(behaviors,{'Attack'}));
        tS_idx_attack = FrameS(idx_attack) - FrameS(idx_IntroM);
        tE_idx_attack = FrameE(idx_attack) - FrameE(idx_IntroM);
        idx_Minv = find(strcmp(behaviors,{'Investigate'}));
        idx_Minv = idx_Minv( idx_Minv > idx_IntroM & idx_Minv < idx_RmvM);
        tS_idx_Minv = FrameS(idx_Minv) - FrameS(idx_IntroM);
        tE_idx_Minv = FrameE(idx_Minv) - FrameE(idx_IntroM);
        %behavior prediction plot
        figure;
        hold on;
        ymulti = 10;
        for ii=1:size(zTestPred,1)
            zTestPred(ii,:) = zscore(zTestPred(ii,:),0,2);
            plot(ii*ymulti+zTestPred(ii,:),'linewidth',0.5,'Color','r');
            plot(ii*ymulti+zscore(male_tracking_data{j}(ii,:),0,2),'linewidth',0.5,'Color','k');
            text(min(xlim),ii*ymulti,strrep(feat_names{select_feature(ii)},'_',' '),'HorizontalAlignment','right');
            temp_corr = corr(zTestPred(ii,:)',zscore(male_tracking_data{j}(ii,:)',0,1));
            text(max(xlim),ii*ymulti,num2str(temp_corr),'HorizontalAlignment','right');
        end
        y_limits = ylim;
        for ii = 1:length(idx_attack)
            patch([tS_idx_attack(ii) tE_idx_attack(ii) tE_idx_attack(ii) tS_idx_attack(ii)], ...
                [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], 'r','EdgeColor','none','FaceAlpha',0.3);
        end
        if ~isempty(idx_Minv)
            for ii = 1:length(idx_Minv)
                patch([tS_idx_Minv(ii) tE_idx_Minv(ii) tE_idx_Minv(ii) tS_idx_Minv(ii)], ...
                    [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], 'b','EdgeColor','none','FaceAlpha',0.3);
            end
        end
        xlabel('time');
        title(strcat(num2str(i), '-',  num2str(j), "-behavior"));
        yticks([]);
        hold off;

        % plot extracted hidden state
        figure;
        hold on;
        ymulti = 10;
        for ii=1:size(xTestPred,1)
            xTestPred(ii,:) = zscore(xTestPred(ii,:),0,2);
            plot(ii*ymulti+xTestPred(ii,:),'linewidth',0.5,'Color','k');
        end
        y_limits = ylim;
        for ii = 1:length(idx_attack)
            patch([tS_idx_attack(ii) tE_idx_attack(ii) tE_idx_attack(ii) tS_idx_attack(ii)], ...
                [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], 'r','EdgeColor','none','FaceAlpha',0.3);
        end
        if ~isempty(idx_Minv)
            for ii = 1:length(idx_Minv)
                patch([tS_idx_Minv(ii) tE_idx_Minv(ii) tE_idx_Minv(ii) tS_idx_Minv(ii)], ...
                    [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], 'b','EdgeColor','none','FaceAlpha',0.3);
            end
        end
        xlabel('time');
        title(strcat(num2str(i), '-',  num2str(j), "-latent state"));
        yticks([]);
        hold off;

    end
end

%%plot correlation matrix
figure;
imagesc(corr_matrix_neural);
colormap('gray');
colorbar;
% caxis([0,1]);
title('GC6-neural prediction');

figure;
imagesc(corr_matrix_behavior);
colormap('gray');
colorbar;
% caxis([0,1]);
title('GC6-behavior prediction');