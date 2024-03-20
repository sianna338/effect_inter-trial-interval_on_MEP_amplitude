%% preliminaries 
addpath('C:\Users\siann\Downloads\fieldtrip-20231113\fieldtrip-20231113')
ft_defaults

datapath = 'C:\Users\siann\Data\biomed_2024\';
% subjects = {'sub01', 'sub02', 'sub03', 'sub04', 'sub05'};
subjects = {'sub03'};
sessions = {'ses1', 'ses2'};
conditions = {'SP', 'SICI', 'SICF'};


%% main analysis: loop over subjects
for subi = 1:numel(subjects)
    for sesi = 1:numel(sessions)
        for condi = 1:numel(conditions )
            % look at markers in dataset and segment based on markers
            cfg = [];
            cfg.datafile = strcat(datapath, subjects{subi}, filesep, 'biomed2024', '_', subjects{subi}, '_', sessions{sesi}, '_', conditions{condi}, '.eeg')
            cfg.continous = 'yes';
            cfg.channel = 'FDIr';
            cfg.trialdef.prestim = 0.05; 
            cfg.trialdef.poststim = 0.1; 
            cfg.trialdef.eventtype = 'Out';
            cfg.trialdef.eventvalue = 'B';
            cfg = ft_definetrial(cfg);
            trial_matrix = cfg.trl
            % read-in data
            % cfg = []; 
            cfg.headerfile = strcat(datapath, 'biomed2024_', subjects{subi}, '_', sessions{sesi}, '_', conditions{condi}, '.vhdr');
            cfg.channel = 'FDIr'; % indicate the channels we would like to read and/or exclude.
            data_raw = ft_preprocessing(cfg);
            save(strcat(datapath,subjects{subi}, filesep,'data_raw_', subjects{subi}, '_', sessions{sesi}, '_', conditions{condi}), '-v7.3');


            % extract min and max points for for the time window from 0.015-0.05 after
            % TMS pulse
            cfg=[];
            cfg.latency=[0.015 0.05];
            cfg.channel = {'FDIr'};
            data_MEPs = ft_selectdata(cfg, data_raw);
            for i = 1:numel(data_MEPs.trial)
                min_val = min(data_MEPs.trial{i}(1,:));
                max_val = max(data_MEPs.trial{i}(1,:));
                data_MEPs.mep(i,1) = abs(min_val) + abs(max_val);
            end
            save(strcat(datapath, subjects{subi}, filesep, 'data_raw_', subjects{subi}, '_', sessions{sesi}, '_', conditions{condi}, '_data_MEPs'), '-v7.3');
        end
    end 
end 
%% get ITIs and MEPs for each condition (concatenate sessions together)
% Initialize variables
data_all = [];
conditions_all = [];
data_MEPs_stored = []
for subi = 1:numel(subjects)
    for condi = 1:numel(conditions)
        for sesi = 1:numel(sessions)
            % Load data
            load(strcat(datapath, subjects{subi}, filesep,'data_raw_', subjects{subi}, '_', sessions{sesi}, '_', conditions{condi}, '_', 'data_MEPs.mat'));
            data_all = [data_all; data_MEPs.mep];
            % Get ITIs on each trial
            load(strcat(datapath, subjects{subi}, filesep, 'BEST_toolbox', filesep, subjects{subi}, filesep, 'Experiment', filesep, 'BESTData_biomed_', subjects{subi}, '_', sessions{sesi}, '_', conditions{condi}));
            conditions_all = [conditions_all; cell2mat(BESTData.trialinfomatrix(:,4))];
            data_MEPs_stored = [data_MEPs_stored; BESTData.Results.FDIr.MEPAmplitude]
        end
        data = [conditions_all, data_all, data_MEPs_stored];
        save(strcat(datapath, subjects{subi}, filesep, 'data_', subjects{subi}, '_', conditions{condi}, '_data_MEPs_ITI'), '-v7.3');
    end
end 
    %% Count how many times each ITI occurs and then get mean MEP amplitude for each ITI
    ITI_counts = [];
    average_MEP = [];
     for subi = 1:length(subjects)
        figure
        for condi = 1:length(conditions)
            load(strcat(datapath, subjects{subi}, filesep, 'data_', subjects{subi}, '_', conditions{condi}, '_data_MEPs_ITI'));
            [unique_ITI, ~, idx] = unique(data(:,1));
            ITI_counts = accumarray(idx, 1);
            average_MEP = accumarray(idx, data(:,2), [], @mean);
            average_MEP_stored = accumarray(idx, data(:,3), [], @mean)
            MEP_avg = struct('ITI', unique_ITI, 'Count', ITI_counts, 'Average_MEP', average_MEP, 'Average_MEP_stored', average_MEP_stored, 'diff', abs(average_MEP_stored-average_MEP));
            save(strcat(datapath, subjects{subi}, filesep, 'data_', subjects{subi}, '_', conditions{condi}, '_data_average_MEPs_ITI'), '-v7.3');

            % Plot
            subplot(1, numel(conditions), condi);
            plot(MEP_avg.ITI, MEP_avg.Average_MEP, '-ok', 'MarkerFaceColor','k');
            title([subjects{subi}, ' ', conditions{condi}, ' MEP amplitudes']);
            xlabel('ITIs');
            ylabel('MEP amplitude');
        end
    end 
