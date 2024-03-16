%% fieldtrip
addpath('C:\Users\siann\Downloads\fieldtrip-20231113\fieldtrip-20231113')
ft_defaults

% look at markers in dataset
cfg = [];
cfg.dataset             = 'biomed2024_sub03_ses1_SP.eeg';
cfg.trialdef.eventtype = '?';
dummy                   = ft_definetrial(cfg);

%% create trial stucture 
cfg = [];
cfg.datafile = 'biomed2024_sub03_ses1_SP.eeg';
cfg.continous = 'yes';
cfg.channel = {'FDIr'};
cfg.trialdef.prestim = 0.05; 
cfg.trialdef.poststim = 0.1; 
cfg.trialdef.eventtype = 'Out';
cfg.trialdef.eventvalue = 'B';
cfg = ft_definetrial(cfg);
trial_matrix = cfg.trl
% read-in data
% cfg = []; 
cfg.datafile = 'biomed2024_sub03_ses1_SP.eeg';
cfg.headerfile = 'biomed2024_sub03_ses1_SP.vhdr';
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
cfg.demean='yes';
cfg.baselinewindow = [-0.05 -0.03];
data_raw_sp_01= ft_preprocessing(cfg);
save('data_raw_sp_01','data_raw_sp_01','-v7.3');

%% create trial stucture 
cfg = [];
cfg.datafile = 'biomed2024_sub03_ses1_SP.eeg';
cfg.continous = 'yes';
cfg.channel = {'FDIr'};
cfg.trialdef.prestim = 0.05; 
cfg.trialdef.poststim = 0.1; 
cfg.trialdef.eventtype = 'Out';
cfg.trialdef.eventvalue = 'B';
cfg = ft_definetrial(cfg);

% read-in data
% cfg = []; 
cfg.datafile = 'biomed2024_sub03_ses2_SP.eeg';
cfg.headerfile = 'biomed2024_sub03_ses2_SP.vhdr';
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
cfg.demean='yes';
cfg.baselinewindow = [-0.05 -0.03];
data_raw_sp_02= ft_preprocessing(cfg);
save('data_raw_sp_02','data_raw_sp_02','-v7.3');

%% create trial stucture 
cfg = [];
cfg.datafile = 'biomed2024_sub03_ses1_SICI.eeg';
cfg.continous = 'yes';
cfg.channel = {'FDIr'};
cfg.trialdef.prestim = 0.05; 
cfg.trialdef.poststim = 0.1; 
cfg.trialdef.eventtype = 'Out';
cfg.trialdef.eventvalue = 'B';
cfg = ft_definetrial(cfg);

% read-in data
% cfg = []; 
cfg.datafile = 'biomed2024_sub03_ses1_SICI.eeg';
cfg.headerfile = 'biomed2024_sub03_ses1_SICI.vhdr';
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
cfg.demean='yes';
cfg.baselinewindow = [-0.05 -0.03];
data_raw_sici_01= ft_preprocessing(cfg);
save('data_raw_sici_01','data_raw_sici_01','-v7.3');

%% create trial stucture 
cfg = [];
cfg.datafile = 'biomed2024_sub03_ses2_SICI.eeg';
cfg.continous = 'yes';
cfg.channel = {'FDIr'};
cfg.trialdef.prestim = 0.05; 
cfg.trialdef.poststim = 0.1; 
cfg.trialdef.eventtype = 'Out';
cfg.trialdef.eventvalue = 'B';
cfg = ft_definetrial(cfg);

% read-in data
% cfg = []; 
cfg.datafile = 'biomed2024_sub03_ses2_SICI.eeg';
cfg.headerfile = 'biomed2024_sub03_ses2_SICI.vhdr';
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
cfg.demean='yes';
cfg.baselinewindow = [-0.05 -0.03];
data_raw_sici_02= ft_preprocessing(cfg);
save('data_raw_sici_02','data_raw_sici_02','-v7.3');

%% create trial stucture 
cfg = [];
cfg.datafile = 'biomed2024_sub03_ses1_SICF.eeg';
cfg.continous = 'yes';
cfg.channel = {'FDIr'};
cfg.trialdef.prestim = 0.05; 
cfg.trialdef.poststim = 0.1; 
cfg.trialdef.eventtype = 'Out';
cfg.trialdef.eventvalue = 'B';
cfg = ft_definetrial(cfg);

% read-in data
% cfg = []; 
cfg.datafile = 'biomed2024_sub03_ses1_SICF.eeg';
cfg.headerfile = 'biomed2024_sub03_ses1_SICF.vhdr';
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
cfg.demean='yes';
cfg.baselinewindow = [-0.05 -0.03];
data_raw_sicf_01= ft_preprocessing(cfg);
save('data_raw_sicf_01','data_raw_sicf_01','-v7.3');

%% create trial stucture 
cfg = [];
cfg.datafile = 'biomed2024_sub03_ses2_SICF.eeg';
cfg.continous = 'yes';
cfg.channel = {'FDIr'};
cfg.trialdef.prestim = 0.05; 
cfg.trialdef.poststim = 0.1; 
cfg.trialdef.eventtype = 'Out';
cfg.trialdef.eventvalue = 'B';
cfg = ft_definetrial(cfg);

% read-in data
% cfg = []; 
cfg.datafile = 'biomed2024_sub03_ses2_SICF.eeg';
cfg.headerfile = 'biomed2024_sub03_ses2_SICF.vhdr';
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
cfg.demean='yes';
cfg.baselinewindow = [-0.05 -0.03];
data_raw_sicf_02= ft_preprocessing(cfg);
save('data_raw_sicf_02','data_raw_sicf_02','-v7.3');

%% extract min and max points for for the time window from 0.015-0.05 after
%% TMS pulse
cfg=[];
cfg.latency=[0.015 0.05];
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
% cfg.trials=find(data.trialinfo==B); 
data_MEPs_sp_01=ft_selectdata(cfg,data_raw_sp_01);
for i=1:size(data_MEPs_sp_01.trial,2) % for each trial (dimension 2)
    min_val=min(data_MEPs_sp_01.trial{1,i}(1,:)); % only one channel (FDI), look through all timepoints in selected window
    max_val=max(data_MEPs_sp_01.trial{1,i}(1,:));
    data_MEPs_sp_01.mep(i,1)=abs(min_val)+abs(max_val); % create list of MEP amplitudes with new row for each trial
end
% data_MEPs_sp_01.mep(81:100,:) = zeros(20,1); % not enough trials here, battery died

%% extract min and max points for for the time window from 0.015-0.05 after TMS pulse
cfg=[];
cfg.latency=[0.015 0.05];
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
% cfg.trials=find(data.trialinfo==B); 
data_MEPs_sp_02=ft_selectdata(cfg,data_raw_sp_02);
for i=1:size(data_MEPs_sp_02.trial,2) % for each trial (dimension 2)
    min_val=min(data_MEPs_sp_02.trial{1,i}(1,:)); % only one channel (APBr), look through all timepoints in selected window
    max_val=max(data_MEPs_sp_02.trial{1,i}(1,:));
    data_MEPs_sp_02.mep(i,1)=abs(min_val)+abs(max_val); % create list of MEP amplitudes with new row for each trial
end

%% extract min and max points for for the time window from 0.015-0.05 after
% TMS pulse
cfg=[];
cfg.latency=[0.015 0.05];
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
% cfg.trials=find(data.trialinfo==B); 
data_MEPs_sici_01=ft_selectdata(cfg,data_raw_sici_01);
for i=1:size(data_MEPs_sici_01.trial,2) % for each trial (dimension 2)
    min_val=min(data_MEPs_sici_01.trial{1,i}(1,:)); % only one channel (APBr), look through all timepoints in selected window
    max_val=max(data_MEPs_sici_01.trial{1,i}(1,:));
    data_MEPs_sici_01.mep(i,1)=abs(min_val)+abs(max_val); % create list of MEP amplitudes with new row for each trial
end

%% extract min and max points for for the time window from 0.015-0.05 after
% TMS pulse
cfg=[];
cfg.latency=[0.015 0.05];
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
% cfg.trials=find(data.trialinfo==B); 
data_MEPs_sici_02=ft_selectdata(cfg,data_raw_sici_02);
for i=1:size(data_MEPs_sici_02.trial,2) % for each trial (dimension 2)
    min_val=min(data_MEPs_sici_02.trial{1,i}(1,:)); % only one channel (APBr), look through all timepoints in selected window
    max_val=max(data_MEPs_sici_02.trial{1,i}(1,:));
    data_MEPs_sici_02.mep(i,1)=abs(min_val)+abs(max_val); % create list of MEP amplitudes with new row for each trial
end

%% extract min and max points for for the time window from 0.015-0.05 after
% TMS pulse
cfg=[];
cfg.latency=[0.015 0.05];
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
% cfg.trials=find(data.trialinfo==B); 
data_MEPs_sicf_01=ft_selectdata(cfg,data_raw_sicf_01);
for i=1:size(data_MEPs_sicf_01.trial,2) % for each trial (dimension 2)
    min_val=min(data_MEPs_sicf_01.trial{1,i}(1,:)); % only one channel (APBr), look through all timepoints in selected window
    max_val=max(data_MEPs_sicf_01.trial{1,i}(1,:));
    data_MEPs_sicf_01.mep(i,1)=abs(min_val)+abs(max_val); % create list of MEP amplitudes with new row for each trial
end

%% extract min and max points for for the time window from 0.015-0.05 after
% TMS pulse
cfg=[];
cfg.latency=[0.015 0.05];
cfg.channel = {'FDIr'}; % indicate the channels we would like to read and/or exclude.
% cfg.trials=find(data.trialinfo==B); 
data_MEPs_sicf_02=ft_selectdata(cfg,data_raw_sicf_02);
for i=1:size(data_MEPs_sicf_02.trial,2) % for each trial (dimension 2)
    min_val=min(data_MEPs_sicf_02.trial{1,i}(1,:)); % only one channel (APBr), look through all timepoints in selected window
    max_val=max(data_MEPs_sicf_02.trial{1,i}(1,:));
    data_MEPs_sicf_02.mep(i,1)=abs(min_val)+abs(max_val); % create list of MEP amplitudes with new row for each trial
end
%% get conditions on each trial
conditions_SP = cell2mat(vertcat(SP_01.trialinfomatrix(:, 4), SP_02.trialinfomatrix(:, 4))); % for some reason SP_ses02 missing
conditions_SICI = cell2mat(vertcat(pp_SICI_01.trialinfomatrix(:, 4), pp_SICI_02.trialinfomatrix(:, 4)));
conditions_SICF = cell2mat(vertcat(pp_SICF_01.trialinfomatrix(:, 4), pp_SICF_02.trialinfomatrix(:, 4)));

%% get MEP amplitude for each condition
data_SP = zeros(200,2); 
data_SP(: ,1) = conditions_SP; 
data_SP(:,2) = vertcat(data_MEPs_sp_01.mep, data_MEPs_sp_02.mep); 

data_SICI = zeros(200,2);
data_SICI(:, 1) = conditions_SICI;
data_SICI(:, 2) = vertcat(data_MEPs_sici_01.mep, data_MEPs_sici_02.mep); 

data_SICF = zeros(200,2);
data_SICF(:, 1) = conditions_SICF; 
data_SICF(:, 2) = vertcat(data_MEPs_sicf_01.mep, data_MEPs_sicf_02.mep);

%% count how many times each ITI occurs and then get mean MEP amplitude for each ITI
[unique_ITI, ~, idx] = unique(data_SP(:,1));
ITI_counts = accumarray(idx, 1);

% Calculate average MEP amplitudes for each unique ITI
average_MEP_sp = accumarray(idx, data_SP(:,2), [], @mean);

% Create a structure to store results
MEP_avg_sp = struct('ITI', unique_ITI, 'Count', ITI_counts, 'Average_MEP', average_MEP_sp);

%% count how many times each ITI occurs and then get mean MEP amplitude for each ITI
[unique_ITI, ~, idx] = unique(data_SP(:,1));
ITI_counts = accumarray(idx, 1);

% Calculate average MEP amplitudes for each unique ITI
average_MEP_sp = accumarray(idx, data_SP(:,2), [], @mean);

% Create a structure to store results
MEP_avg_sp = struct('ITI', unique_ITI, 'Count', ITI_counts, 'Average_MEP', average_MEP_sp);

%% count how many times each ITI occurs and then get mean MEP amplitude for each ITI
[unique_ITI, ~, idx] = unique(data_SICI(:,1));
ITI_counts = accumarray(idx, 1);

% Calculate average MEP amplitudes for each unique ITI
average_MEP_sici = accumarray(idx, data_SICI(:,2), [], @mean);

% Create a structure to store results
MEP_avg_sici = struct('ITI', unique_ITI, 'Count', ITI_counts, 'Average_MEP', average_MEP_sici);

%% count how many times each ITI occurs and then get mean MEP amplitude for each ITI
[unique_ITI, ~, idx] = unique(data_SICF(:,1));
ITI_counts = accumarray(idx, 1);

% Calculate average MEP amplitudes for each unique ITI
average_MEP_sicf = accumarray(idx, data_SICF(:,2), [], @mean);

% Create a structure to store results
MEP_avg_sicf = struct('ITI', unique_ITI, 'Count', ITI_counts, 'Average_MEP', average_MEP_sicf);

%% plot
figure; 
plot(MEP_avg_sp.ITI, MEP_avg_sp.Average_MEP, '-ok', 'MarkerFaceColor','k')
title('single pulse MEP amplitudes')
xlabel('ITIs')
ylabel('MEP amplitude')

figure; 
plot(MEP_avg_sici.ITI, MEP_avg_sici.Average_MEP, '-ok', 'MarkerFaceColor','k')
title('SICI MEP amplitudes')
xlabel('ITIs')
ylabel('MEP amplitude')

figure; 
plot(MEP_avg_sicf.ITI, MEP_avg_sicf.Average_MEP, '-ok', 'MarkerFaceColor','k')
title('SICF MEP amplitudes')
xlabel('ITIs')
ylabel('MEP amplitude')

%% ANOVA
% [p,tbl,stats] = anova1(data_MEPs.mep, conditions);

