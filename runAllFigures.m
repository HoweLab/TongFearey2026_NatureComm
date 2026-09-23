% runAllFigures.m
% Generates every figure panel (main Figs 1-5, Extended Data Figs 2-9) and
% the statistics / source-data workbook (one tab per figure).
%
% All inputs are set in the "Inputs" section below; nothing else needs to be
% edited. Each figure is a function in figures/ that takes this cfg struct.

clear; close all

%% ===================== Inputs: edit here =====================
dataRoot = 'E:\';   % where you save the processed files

% Output: figures go to <outDir>\<figure name>\, stats to one workbook
cfg.outDir    = fullfile(pwd, 'output');
cfg.statsFile = fullfile(cfg.outDir, 'PaperStats_report.xlsx');

% Linear track, 1 world, and spontaneous locomotion (GUI) - Fig 1-2, Sup 2-5
cfg.track1w = fullfile(dataRoot, 'processedData_2025final', 'processedData_convolv', 'processedDataTrack_all.mat');
cfg.gui     = fullfile(dataRoot, 'processedData_2025final', 'processedData_convolv', 'processedData_gui.mat');
% speed-tuning correlations - Fig 1G, Sup 2
cfg.speedTuning1w  = fullfile(dataRoot, 'Vel_Encoding_SX', 'SpeedTuning_1w.mat');
cfg.speedTuningGui = fullfile(dataRoot, 'Vel_Encoding_SX', 'SpeedTuning_GUI.mat');
% summary ROI folders (cells per cohort) - Sup 3I
cfg.roiTrack = fullfile(dataRoot, 'SummaryROI', 'linearTrackRecording_24to25_1w_deconv');
cfg.roiGui   = fullfile(dataRoot, 'SummaryROI', 'linearTrackRecording_24to25_gui_deconv');
% percent active cells - Sup 4I-P
cfg.prctActive1w  = fullfile(dataRoot, 'processedData_2025final', 'processedData_prctActive', 'processedData_1w.mat');
cfg.prctActiveGui = fullfile(dataRoot, 'processedData_2025final', 'processedData_prctActive', 'processedData_gui.mat');

% Two tracks: tuning subtypes and single-cell data - Fig 3, Sup 6-7
cfg.tuned2w    = fullfile(dataRoot, 'processedData_2025test', 'tuningSubType_2026_2w.mat');
cfg.sc2w       = fullfile(dataRoot, 'dataOrganizedByCell', 'bf2w_2025.mat');
cfg.switchInfo = fullfile(dataRoot, 'notes', 'NovFam_Switch_Info.xlsx');   % also Fig 4, Sup 8

% Novel vs. familiar track (all sessions) - Fig 4, Sup 8
cfg.trackNovFam = fullfile(dataRoot, 'processedData_2025final', 'processedDataTrack_all.mat');

% Infinite track: population, tuning subtypes, single-cell - Fig 5, Sup 9
cfg.infPop   = fullfile(dataRoot, 'processedData_2025final', 'processedData_all.mat');
cfg.infTuned = fullfile(dataRoot, 'processedData_2025test', 'tuningSubType_2025_inf_newtype.mat');
cfg.infSc    = fullfile(dataRoot, 'dataOrganizedByCell', 'bfinf_DistTuned_2025.mat');

% Figures to generate (remove entries to run a subset)
figsToRun = {'Figure1', 'Figure2', 'Figure3', 'Figure4', 'Figure5', ...
             'Sup2', 'Sup3', 'Sup4', 'Sup5', 'Sup6', 'Sup7', 'Sup8', 'Sup9'};

%% ===================== Setup =====================
repoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(repoDir, 'figures'), fullfile(repoDir, 'figTools'), ...
        fullfile(repoDir, 'reportTools'), fullfile(repoDir, 'commonFunction'));

% fail early on a missing input rather than partway through a long run
inputs = rmfield(cfg, {'outDir', 'statsFile'});
names = fieldnames(inputs);
missing = names(~cellfun(@(f) isfile(inputs.(f)) || isfolder(inputs.(f)), names));
if ~isempty(missing)
    error('Input not found: %s', strjoin(strcat(missing, ' = ', ...
        cellfun(@(f) inputs.(f), missing, 'UniformOutput', false)), '; '));
end
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
if isempty(gcp('nocreate')), parpool('Threads'); end

%% ===================== Run =====================
for k = 1:numel(figsToRun)
    fprintf('\n===== %s =====\n', figsToRun{k});
    t = tic;
    close all
    feval(figsToRun{k}, cfg);
    fprintf('%s done (%.1f min)\n', figsToRun{k}, toc(t)/60);
end
