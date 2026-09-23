function Figure1(cfg)
% Main Figure 1: D-G - D and I dynamic & correlation w/ kinematics
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure
% processed data for each session
fpath = cfg.track1w;
corrpath = cfg.speedTuning1w;
% organize data
conditions = {'exp_stage',{'1 world'};
    'objective','20x';
    'sess_type','track';
    'mismatch',0};
track_tbl = mkTblForPlot(fpath,conditions);

% data selection
mice_include = {'tdt','TdTG7', 'a2a'}; % all 3 batches
trial_selection = ~track_tbl.skip_trial & track_tbl.exp_stage == '1 world' ...
                & contains(cellstr(track_tbl.mouse),mice_include);
bout_selection = track_tbl.info_bout_intrack == 1 & contains(cellstr(track_tbl.mouse),mice_include);  
xdata = [1:100; linspace(1,100,100)];
info = struct('startLoc',1,'rewardLoc',100);
% Trig window
on_win = [-31:3*31; linspace(-1,3,125)];
off_win = [3*-31:31; linspace(-3,1,125)];

%% start generate figures for fig1
figDir = fullfile(save2where, 'Fig1');
if ~isfolder(figDir), mkdir(figDir); end

vuBin_tbl = groupsummary(track_tbl(trial_selection,:), ... % track 
    {'mouse','sess','field','trial','vuBin','isD','isI', ...
    'celltype','numFc3'}, ...  
    'mean',{'numFc3','avgFc3','velocity','acceleration','info_trialTimes', 'lick_count'});
binned_velAcc = getPlotDatFromTbl(vuBin_tbl,'vuBin','y','mean_acceleration',  ...
    'vel_y','mean_velocity','xdata',xdata);
binned_data = getPlotDatFromTbl(vuBin_tbl,'vuBin','y','mean_avgFc3', ...
    'vel_y','mean_acceleration','xdata',xdata);
onsets = getPlotDatFromTbl(track_tbl(trial_selection & bout_selection,:), ...
    'onsets','xdata',on_win,'vel_y','velocity');
offsets = getPlotDatFromTbl(track_tbl(trial_selection & bout_selection,:), ...
    'offsets','xdata',off_win,'vel_y','velocity');
trig_dat = struct('onsets',onsets, 'offsets', offsets);

%%% Plotting %%%
% Track velocity
numMouse = numel(unique(vuBin_tbl.mouse));

figure('Position',[900 100 900 700])
plotBinned_new(binned_velAcc,info, ...
    'fromLME',true,'traces',{'velocity'}, 'line_colors',{[0 0 1]})
[~, vel_mouse] = addSingleMouseLinePlot(gca, vuBin_tbl, 'vuBin', 'mean_velocity', 'mean_v', 0, numMouse);
save_img(gcf,figDir,'Fig1_D_bottom_1world_positionBins_velocity'); close(gcf)

% Track acceleration
figure('Position',[900 100 900 700])
plotBinned_new(binned_velAcc,info, ...
    'fromLME',true,'traces',{'dSPN'}, 'line_colors',{[0 0 1]})
[~, acc_mouse] = addSingleMouseLinePlot(gca, vuBin_tbl, 'vuBin', 'mean_acceleration', 'mean_accl', 0, numMouse);
save_img(gcf,figDir,'Fig1_D_top_1world_positionBins_acceleration'); close(gcf)

% Track dspn/ispn seperately
figure('Position',[900 100 900 700])
plotBinned_new(binned_data,info, ...
    'fromLME',true,'traces',{'dSPN'},'plot_vel',false)
save_img(gcf,figDir,'Fig1_E_top_1world_positionBins_dSPN'); close(gcf)
plotBinned_new(binned_data,info, ...
    'fromLME',true,'traces',{'iSPN'}, 'plot_vel',false)
save_img(gcf,figDir,'Fig1_E_bottom_1world_positionBins_iSPN'); close(gcf)

% Triggered D and I
figure('Position',[100 100 900 700])
plotTriggered_new(trig_dat,'traces',{'dSPN'},'fromLME',true, ...
    'plot_vel',false, 'trig_events',fieldnames(trig_dat));
save_img(gcf,figDir,'Fig1_F_top_1world_positionBins_dSPN'); close(gcf)
figure('Position',[100 100 900 700])
plotTriggered_new(trig_dat,'traces',{'iSPN'},'fromLME',true, ...
    'plot_vel',false, 'trig_events',fieldnames(trig_dat));
save_img(gcf,figDir,'Fig1_F_middle_1world_positionBins_iSPN'); close(gcf)

% Triggered velocity
figure('Position',[900 100 900 700])
plotTriggered_new(trig_dat,'traces',{'velocity'},'fromLME',true, ...
    'line_colors',{[0 0 1]},'plot_vel',false, 'trig_events',fieldnames(trig_dat));
save_img(gcf,figDir,'Fig1_F_bottom_1world_triggeredAverageVelocity'); close(gcf)

% correlation figure
track_corr = load(corrpath, 'dspn_cor','ispn_cor');
vals = {'binned_v_eachN_cor', 'binned_v_sess_cor'};
for i = 1:numel(vals)    
    [tho_d_1,p_d_1,tho_i_1,p_i_1]=  grab_celldat(track_corr, vals{i});
    % dspn left; ispn right
    figure('Name','dspn vs. ispn track')
    p1 = boxViolin_colorSig(0.001, tho_d_1, p_d_1, tho_i_1, p_i_1,  get(gcf,'Number'));
    title(sprintf('%s %s; G1 = dspn, G2 = ispn, p = %.4g','track', vals{i}, p1), 'Interpreter','none')
    save_img(gcf, figDir, strcat('Fig1_G_DvsI_','track', '_', vals{i}))
end
close all
%% Count sizes & report stats -> Excel tab 'Fig1'
fig = 'Fig1';
n_trial = countSampleSize(track_tbl(trial_selection,:), 'trial');                 % D, E
n_bout  = countSampleSize(track_tbl(trial_selection & bout_selection,:), 'bout'); % F

R = [statsRow('Figure',fig,'Panel','D-E','Measure','sample size (track traversal, position-binned)', ...
        'Group','dSPN','nStruct',n_trial)
     statsRow('Figure',fig,'Panel','D-E','Measure','sample size (track traversal, position-binned)', ...
        'Group','iSPN','nStruct',n_trial)
     statsRow('Figure',fig,'Panel','F','Measure','sample size (in-track locomotion bouts)', ...
        'Group','dSPN','nStruct',n_bout)
     statsRow('Figure',fig,'Panel','F','Measure','sample size (in-track locomotion bouts)', ...
        'Group','iSPN','nStruct',n_bout)];

% G: box plots (median/IQR), % significant and dSPN vs iSPN rank-sum (the tests shown)
corr_desc = struct('binned_v_eachN_cor','Spearman rho, velocity-binned dF/F vs velocity, per neuron', ...
                   'binned_v_sess_cor', 'Spearman rho, velocity-binned mean dF/F vs velocity, per field');
for i = 1:numel(vals)
    [tho_d_1,p_d_1,tho_i_1,p_i_1] = grab_celldat(track_corr, vals{i});
    n_d = countCorrSampleSize(track_corr.dspn_cor, vals{i}, track_tbl(trial_selection,:));
    n_i = countCorrSampleSize(track_corr.ispn_cor, vals{i}, track_tbl(trial_selection,:));
    R = [R; reportGroupCompare({'dSPN','iSPN'}, {tho_d_1, tho_i_1}, ...
            'pvals',{p_d_1, p_i_1}, 'pThresh',0.001, 'n',{n_d, n_i}, ...
            'Figure',fig, 'Panel','G', 'Measure',corr_desc.(vals{i}))]; %#ok<AGROW>
end

% plotted values (no tests in D-F: LME mean and SEM as drawn)
mouse_bins = unique(vuBin_tbl.vuBin(~isnan(vuBin_tbl.vuBin)));  % bins used by addSingleMouseLinePlot
blocks = {
    'Summary: sample sizes; G box plots, % significant, dSPN vs iSPN rank-sum', R
    'Fig1D top: acceleration (LME mean, SEM) by track position', ...
        sourceDataCurve(binned_velAcc, {'dSPN'}, 'labels',{'acceleration'}, 'xName','position_bin')
    'Fig1D top: acceleration, mean per mouse (gray lines)', ...
        sourceDataPerMouse(acc_mouse, mouse_bins, 'xName','position_bin')
    'Fig1D bottom: velocity (LME mean, SEM) by track position', ...
        sourceDataCurve(binned_velAcc, {'velocity'}, 'xName','position_bin')
    'Fig1D bottom: velocity, mean per mouse (gray lines)', ...
        sourceDataPerMouse(vel_mouse, mouse_bins, 'xName','position_bin')
    'Fig1E: dSPN and iSPN dF/F (LME mean, SEM) by track position', ...
        sourceDataCurve(binned_data, {'dSPN','iSPN'}, 'xName','position_bin')
    'Fig1F: onset-aligned dSPN, iSPN dF/F and velocity (LME mean, SEM)', ...
        sourceDataCurve(trig_dat.onsets, {'dSPN','iSPN','velocity'}, 'xName','time_s')
    'Fig1F: offset-aligned dSPN, iSPN dF/F and velocity (LME mean, SEM)', ...
        sourceDataCurve(trig_dat.offsets, {'dSPN','iSPN','velocity'}, 'xName','time_s')
    'Fig1G top: population correlation per imaging field (red = p<0.001)', ...
        [sourceDataCorrPoints(track_corr.dspn_cor, 'binned_v_sess_cor', 'dSPN', 0.001)
         sourceDataCorrPoints(track_corr.ispn_cor, 'binned_v_sess_cor', 'iSPN', 0.001)]
    'Fig1G bottom: single-neuron correlation (red = p<0.001)', ...
        [sourceDataCorrPoints(track_corr.dspn_cor, 'binned_v_eachN_cor', 'dSPN', 0.001)
         sourceDataCorrPoints(track_corr.ispn_cor, 'binned_v_eachN_cor', 'iSPN', 0.001)]
    'Sample size per mouse: D-E (trials)', n_trial.perMouse
    'Sample size per mouse: F (bouts)', n_bout.perMouse};
writeSourceData(statsFile, fig, blocks);

end

%% Helper
function [tho_d,p_d,tho_i,p_i]=  grab_celldat(corr_struct, varname)
    tho_d = cell2mat(cellfun(@(x) x.tho, {corr_struct.dspn_cor.(varname)}, 'UniformOutput', false)'); 
    p_d = cell2mat(cellfun(@(x) x.p, {corr_struct.dspn_cor.(varname)}, 'UniformOutput', false)'); 
    tho_i = cell2mat(cellfun(@(x) x.tho, {corr_struct.ispn_cor.(varname)}, 'UniformOutput', false)'); 
    p_i = cell2mat(cellfun(@(x) x.p, {corr_struct.ispn_cor.(varname)}, 'UniformOutput', false)'); 
end