%% Notes:
% This is a script using GPFA to visualize the partern PSID predict X, and
% to find the pattern of trajectures' angle feature
%author: yibochen, 10/10/23

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
    6, '200313';
    6, '200314';
    6, '200315';
    };

% params
choose_behavior = 'Thrust'; % 'Mount' or 'Thrust'
if strcmp(choose_behavior,'Mount')%mount 1:4,thrust 3:4
    select_session = 1:4;
else
    select_session = 3:4;
end
select_feature = [10,20,21,24,25,27,29];

%% initialize
female_neural_data = cell(length(select_session),1);
female_tracking_data = cell(length(select_session),1);
angle_sessions = cell(length(select_session),1);
attack_durations = cell(length(select_session),1);
for ii = 1:length(select_session)
    %load MFP data
    load(['../Data/MFP_raw/MFP-ERa-GC', num2str(sessions{select_session(ii),1}), '.mat'])
    RawS = Raw.(['r',sessions{select_session(ii),2}]);
    neural_data = zscore(RawS.Lfold,0,2);
    behaviors = RawS.behaviors;
    FrameS = RawS.Fstart;
    FrameE = RawS.Fstop;
    idx_IntroF = find(strcmp(behaviors,{'Intro_F'}));
    idx_RmvF = find(strcmp(behaviors,{'Rmv_F'}));
    female_idx = FrameS(idx_IntroF(end)):FrameE(idx_RmvF(end));
    female_neural_data{ii} = neural_data(:,female_idx);

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
    female_tracking_data{ii} = tracking_data(:,female_idx);
    cd(curr_dir)
end

%% PSID
% params
%initilization
nx = 13;
n1 = 6;
i_para = 10;
for avewndw = 10
    for lag = 3
        %k-fold validation for i-i prediction
        for i = 1:length(select_session)
            for j = 1:length(select_session)
                if i~=j
                    continue
                end
                Sys_female = PSID(female_neural_data{i}',female_tracking_data{i}',nx,n1,i_para);
                [zTestPred, yTestPred, xTestPred] = PSIDPredict(Sys_female, female_neural_data{j}');
                zTestPred = zTestPred';
                yTestPred = yTestPred';
                xTestPred = xTestPred';
        
                load(['../Data/MFP_raw/MFP-ERa-GC', num2str(sessions{select_session(j),1}), '.mat'])
                RawS = Raw.(['r',sessions{select_session(j),2}]);
                behaviors = RawS.behaviors;
                FrameS = RawS.Fstart;
                FrameE = RawS.Fstop;
                idx_IntroF = find(strcmp(behaviors,{'Intro_F'}));
                idx_RmvF = find(strcmp(behaviors,{'Rmv_F'}));
                idx_mount = find(strcmp(behaviors,{choose_behavior}));
                idx_mount(idx_mount>idx_RmvF(end)) = [];
                tS_idx_mount = FrameS(idx_mount) - FrameS(idx_IntroF(end));
                tE_idx_mount = FrameE(idx_mount) - FrameS(idx_IntroF(end));
        
                color = [1,0,0;
                        0,1,0;
                        0,0,1;
                        0,0,0;
                        1,0,1;
                        0,1,1;
                        ];
                max_lines = min(6,length(tS_idx_mount));
                angles = zeros(1,max_lines);
                duration = zeros(1,max_lines);
                legends = cell(1,max_lines);
                handles = zeros(1,max_lines);
                f=figure;
                hold on;
                select_trajectory = 1:max_lines;
                for ii = 1:max_lines
                    [~,Sz,Vz] = svd(Sys_female.Cz);
                    temp = xTestPred(:,tS_idx_mount(select_trajectory(ii))-avewndw*lag:min(size(xTestPred,2),tS_idx_mount(select_trajectory(ii))+avewndw*lag-1));
                    bin_pred = [];
                    for k = 1:size(xTestPred,1)
                       n = length(temp(k,:));
                       reshaped_matrix = reshape(temp(k,:),[avewndw, n/avewndw]);
                       averaged_vector = mean(reshaped_matrix);
                       bin_pred = [bin_pred;averaged_vector];
                    end
                    manifold_z = Sz*Vz'*bin_pred;
                    % only consider head and tail
                    manifold_z = manifold_z - manifold_z(:,lag+1);
                    v1 = manifold_z(:,1);
                    v2 = manifold_z(:,end);
                    cos_theta = (v1'*v2)/(sqrt(sum(v1.^2))*sqrt(sum(v2.^2)));
                    angles(1,ii) = acosd(cos_theta);
                    %   2D
                    handles(ii) = line(manifold_z(1,:),manifold_z(2,:),'color',color(mod(ii-1,size(color,1))+1,:));
                    plot(manifold_z(1,:),manifold_z(2,:),'Marker','+','Color',color(mod(ii-1,size(color,1))+1,:),'MarkerSize',3)
                    plot(manifold_z(1,1),manifold_z(2,1),'o','Color',color(mod(ii-1,size(color,1))+1,:),'MarkerFaceColor','auto','MarkerSize',10);
                    plot(manifold_z(1,end),manifold_z(2,end),'^','Color',color(mod(ii-1,size(color,1))+1,:),'MarkerFaceColor','auto','MarkerSize',10);
                    duration(1,ii) = -tS_idx_mount(select_trajectory(ii))+tE_idx_mount(select_trajectory(ii));
                    legends{ii} = num2str(angles(1,ii));
                end
                legend(handles,legends);
                axis equal;
                hold off;
                angle_sessions{i} = angles;
                attack_durations{i} = duration;
            end
        end
        %% Draw figures
        total_angle = 0;
        total_duration = 0;
        count = 0;
        for i = 1:length(angle_sessions)
            total_angle = [total_angle,angle_sessions{i}];
            total_duration = [total_duration,attack_durations{i}];
        end
        %T-test whether greater than 90
        [h, p, ci, stats] = ttest(total_angle, 90);
        if h == 1 && stats.tstat > 0
            disp('Significantly greater than 90');
            disp(['p-value: ', num2str(p)]);
        elseif h == 1 && stats.tstat < 0
            disp('Significantly less than 90');
            disp(['p-value: ', num2str(p)]);
        else
            disp('Not significantly different from 90');
            disp(['p-value: ', num2str(p)]);
        end
        
        data = total_angle;
        
        % histograph
        figure;
        histogram(data, 'FaceColor', 'blue','BinWidth',10);
        hold on;
        ylim([0,50])
        yL = ylim; 
        line([90 90], yL, 'Color', 'red', 'LineWidth', 2);
        title('Histogram with Reference Line');
        xlabel('Degree','FontSize',20,'FontWeight','bold');
        ylabel('Number of Trace','FontSize',20,'FontWeight','bold');
        
        % box graph
        axes('Position',[0.2 0.7 0.15 0.2]);
        boxplot(data,'Widths',1.0);
        hold on;
        ylim([-5,185])
        yticks([0,45,90,135,180]);
        yL = ylim; 
        xL = xlim;
        line([min(xL) max(xL)], [90 90], 'Color', 'red', 'LineWidth', 2);
        title('Boxplot with Reference Line');
        ylabel('Degree');
        set(gca,"XTick",[]);
        
        
        % 95% Confidence Interval for Mean
        figure;
        mu = mean(data);
        sigma = std(data);
        n = length(data);
        SEM = sigma/sqrt(n);
        ts = tinv([0.025  0.975],n-1);
        CI = mu + ts*SEM;
        errorbar(1, mu, mu-CI(1), CI(2)-mu, 'ob');
        hold on;
        line([0.5 1.5], [90 90], 'Color', 'red', 'LineWidth', 5);
        xlim([0.5 1.5]);
        title('95% Confidence Interval for Mean');
        ylabel('Degree');
        set(gca,"XTick",[]);
        
        [rho, pval]=corr(total_angle',total_duration','type','Spearman');
        disp(["Degree mean: ",num2str(mu),"median: ",num2str(median(data))]);
        disp(["lag: ",num2str(lag),"avewndw: ", num2str(avewndw),"corr: ",num2str(rho),"pval: ",num2str(pval)])
        %scatter
        figure;
        hold on;
        idx = total_duration > 250;
        total_duration(idx) = [];
        total_angle(idx) = [];
        disp(["point number: ",num2str(length(total_duration))]);
        scatter(total_duration, total_angle, 200, '.');
        [b,stats] = robustfit(total_duration,total_angle);
        x_fit = linspace(min(total_duration),max(total_duration));
        y_fit = b(1) + b(2) * x_fit;
        plot(x_fit,y_fit,'r','LineWidth',2);
        xlabel("Duration","FontSize",20);
        ylabel("Degree",'FontSize',20);
        xlim([-5,max(total_duration)]);
        ylim([-5,max(total_angle)]); 
        pbaspect([max(total_duration) max(total_angle) 1]);
        hold off;
    end
end