function Sup2(cfg)
% Sup2: supplement of main Fig 1;
% Example mouse, licking, and GUI Di/kinematics
% (A-C and L are single-session / single-cell examples: not generated here)
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure
fpath_gui = cfg.gui;
fpath_1w = cfg.track1w;
corrpath_gui = cfg.speedTuningGui;

conditions = {'sess_type','gui'; 'mismatch',0};
gui_tbl = mkTblForPlot(fpath_gui,conditions);
conditions = {'exp_stage',{'1 world'};
    'objective','20x';
    'sess_type','track';
    'mismatch',0};
track_tbl = mkTblForPlot(fpath_1w,conditions);

mice_include = {'tdt','TdTG7', 'a2a'}; % all 3 batches
gui_selection = gui_tbl.info_bout_rewarded ~= 1  &... %
    gui_tbl.info_bout_length > 0 & gui_tbl.info_peak_vel > 15 & ...
    contains(cellstr(gui_tbl.mouse),mice_include); % unrewarded bouts
trial_selection = ~track_tbl.skip_trial & track_tbl.exp_stage == '1 world' ...
                & contains(cellstr(track_tbl.mouse),mice_include);
bout_selection = track_tbl.info_bout_intrack == 1 & contains(cellstr(track_tbl.mouse),mice_include);

xdata = [1:100; linspace(1,100,100)];
fake_info = struct('startLoc',1,'rewardLoc',100);

% Trig window
sr = 31;
on_win = [-31:3*31; linspace(-1,3,125)];
off_win = [3*-31:31; linspace(-3,1,125)];
lick_win = [-2*31:0; linspace(-2,0,63)]; % reward-aligned licking

%% start generate figures for Sup2
figDir = fullfile(save2where, 'Sup2');
if ~isfolder(figDir), mkdir(figDir); end

% (from paperFigScripts/paperFigurePopulation_1w_gui.m)
lick_trig = getPlotDatFromTbl(track_tbl(trial_selection,:), ...
    'rewards','xdata',lick_win,'y','lick_count');
lick_dat = struct('onsets',lick_trig);
pBin_tbl = groupsummary(gui_tbl(gui_selection,:), ...      % gui
    {'mouse','sess','field','bout','progressBin','isD','isI','celltype','numFc3'}, ...
    'mean',{'avgFc3','velocity','acceleration'});
gui_velAcc = getPlotDatFromTbl(pBin_tbl,'progressBin','y','mean_acceleration',  ...
    'vel_y','mean_velocity','xdata',xdata);
gui_binned = getPlotDatFromTbl(pBin_tbl,'progressBin','y','mean_avgFc3', ...
    'vel_y','mean_acceleration','xdata',xdata);
gui_onsets = getPlotDatFromTbl(gui_tbl(gui_selection,:), ...
    'onsets','xdata',on_win,'vel_y','velocity');
gui_offsets = getPlotDatFromTbl(gui_tbl(gui_selection,:), ...
    'offsets','xdata',off_win,'vel_y','velocity');
gui_trig_dat = struct('onsets',gui_onsets, 'offsets', gui_offsets);

%%% Plotting %%%
% D top: reward-aligned licking, model estimate across mice/sessions
figure('Position',[900 100 900 700])
plotTriggered_new(lick_dat,'traces',{'dSPN'},'fromLME',true, ...
    'line_colors',{[0 0 1]},'plot_vel',false,'trig_events',{'onsets'});
save_img(gcf,figDir,'Sup2_D_top_1world_triggered_lick'); close(gcf)

% D bottom: individual mouse averages
lick_tbl = track_tbl(trial_selection & track_tbl.rewards <= 31 & track_tbl.rewards >=-2*31,:);
numMouse = numel(unique(lick_tbl.mouse));
[~, lick_mouse] = addSingleMouseLinePlot([], lick_tbl, 'rewards', 'lick_count','lick_count', 0, numMouse);
save_img(gcf, figDir,'Sup2_D_bottom_1world_lickcount'); close(gcf)

% F: bout peak velocity distribution with the 15 cm/s selection cut-off
[~, bout_all] = bout_descriptive(gui_tbl,31.25,1.5,0);
figure
h_pkvel = histogram(bout_all.pkvel,50);
T_pkvel_hist = sourceDataHist(h_pkvel); % read bars before the figure is closed
hold on; xline(15,'LineWidth',2,'Color','r')
xlabel('cm/s')
title('Bout peak velocity distribution and selection; >15cm/s peak vel')
save_img(gcf,figDir,'Sup2_F_gui_bout_peakvelDistribution'); close(gcf)

% E: selected bout duration distribution
gui_tbl_E = gui_tbl(gui_selection & gui_tbl.isD == 1,:);
[~, bout_sel] = bout_descriptive(gui_tbl_E,31.25,1.5,0);
figure
h_dur = histogram(bout_sel.duration,50);
T_dur_hist = sourceDataHist(h_dur); % read bars before the figure is closed
xlabel('s')
title('GUI selected bout duration distribution')
save_img(gcf,figDir,'Sup2_E_gui_bout_duration_distribution'); close(gcf)

% G, H: GUI velocity / acceleration by bout progress (+ single mice)
numMouse_gui = numel(unique(pBin_tbl.mouse));
% addSingleMouseLinePlot groups by 'trial'; for GUI data each bout is the unit
pBin_mouse = pBin_tbl;
pBin_mouse.trial = pBin_mouse.bout;
figure('Position',[900 100 900 700])
plotBinned_new(gui_velAcc,fake_info, ...
    'fromLME',true,'traces',{'velocity'}, 'line_colors',{[0 0 1]})
[~, vel_mouse] = addSingleMouseLinePlot(gca, pBin_mouse, 'progressBin', 'mean_velocity', 'mean_v', 0, numMouse_gui);
save_img(gcf,figDir,'Sup2_G_gui_boutProgressBin_velocity'); close(gcf)

figure('Position',[900 100 900 700])
plotBinned_new(gui_velAcc,fake_info, ...
    'fromLME',true,'traces',{'dSPN'}, 'line_colors',{[0 0 1]})
[~, acc_mouse] = addSingleMouseLinePlot(gca, pBin_mouse, 'progressBin', 'mean_acceleration', 'mean_accl', 0, numMouse_gui);
save_img(gcf,figDir,'Sup2_H_gui_boutProgressBin_acceleration'); close(gcf)

% I: GUI dspn/ispn seperately by bout progress
figure('Position',[900 100 900 700])
plotBinned_new(gui_binned,fake_info, ...
    'fromLME',true,'traces',{'dSPN'},'plot_vel',false)
save_img(gcf,figDir,'Sup2_I_top_gui_boutProgressBin_dSPN'); close(gcf)
figure('Position',[900 100 900 700])
plotBinned_new(gui_binned,fake_info, ...
    'fromLME',true,'traces',{'iSPN'}, 'plot_vel',false)
save_img(gcf,figDir,'Sup2_I_bottom_gui_boutProgressBin_iSPN'); close(gcf)

% J: GUI triggered D, I and velocity
figure('Position',[100 100 900 700])
plotTriggered_new(gui_trig_dat,'traces',{'dSPN'},'fromLME',true, ...
    'plot_vel',false, 'trig_events',fieldnames(gui_trig_dat));
save_img(gcf,figDir,'Sup2_J_top_gui_trig_dSPN'); close(gcf)
figure('Position',[100 100 900 700])
plotTriggered_new(gui_trig_dat,'traces',{'iSPN'},'fromLME',true, ...
    'plot_vel',false, 'trig_events',fieldnames(gui_trig_dat));
save_img(gcf,figDir,'Sup2_J_middle_gui_trig_iSPN'); close(gcf)
figure('Position',[900 100 900 700])
plotTriggered_new(gui_trig_dat,'traces',{'velocity'},'fromLME',true, ...
    'line_colors',{[0 0 1]},'plot_vel',false, 'trig_events',fieldnames(gui_trig_dat));
save_img(gcf,figDir,'Sup2_J_bottom_gui_triggeredAverageVelocity'); close(gcf)

% K: GUI correlation figure
gui_corr = load(corrpath_gui, 'dspn_cor','ispn_cor');
vals = {'binned_v_eachN_cor', 'binned_v_sess_cor'};
for i = 1:numel(vals)
    [tho_d_1,p_d_1,tho_i_1,p_i_1]=  grab_celldat(gui_corr, vals{i});
    % dspn left; ispn right
    figure('Name','dspn vs. ispn gui')
    p1 = boxViolin_colorSig(0.001, tho_d_1, p_d_1, tho_i_1, p_i_1,  get(gcf,'Number'));
    title(sprintf('%s %s; G1 = dspn, G2 = ispn, p = %.4g','gui', vals{i}, p1), 'Interpreter','none')
    save_img(gcf, figDir, strcat('Sup2_K_DvsI_','gui', '_', vals{i}))
end
close all

%% Count sizes & report stats -> Excel tab 'Sup2'
fig = 'Sup2';
n_trial = countSampleSize(track_tbl(trial_selection,:), 'trial');   % D
n_gui   = countSampleSize(gui_tbl(gui_selection,:), 'bout');        % E, G-J
n_guiAll = countSampleSize(gui_tbl, 'bout');                        % F

R = table();
ss = {'D', 'track traversal trials (licking)', n_trial;
      'E, G-J', 'selected spontaneous locomotion (GUI) bouts', n_gui;
      'F', 'all spontaneous locomotion (GUI) bouts', n_guiAll};
for s = 1:size(ss,1)
    for ct = {'dSPN','iSPN'}
        R = [R; statsRow('Figure',fig,'Panel',ss{s,1},'Measure',['sample size (' ss{s,2} ')'], ...
                'Group',ct{1},'nStruct',ss{s,3})]; %#ok<AGROW>
    end
end
R = [R; statsRow('Figure',fig,'Panel','E','Measure','bouts in histogram', ...
        'N',height(bout_sel),'N_unit',"bouts",'Note','dSPN half of the table only (one row per bout)')
     statsRow('Figure',fig,'Panel','F','Measure','bouts in histogram', ...
        'N',height(bout_all),'N_unit',"bouts", ...
        'Note','all GUI bouts; bout_descriptive on the full table counts each bout in its dSPN and iSPN halves')];

% K: box plots (median/IQR), % significant and dSPN vs iSPN rank-sum (the tests shown)
corr_desc = struct('binned_v_eachN_cor','Spearman rho, velocity-binned dF/F vs velocity, per neuron', ...
                   'binned_v_sess_cor', 'Spearman rho, velocity-binned mean dF/F vs velocity, per field');
for i = 1:numel(vals)
    [tho_d_1,p_d_1,tho_i_1,p_i_1] = grab_celldat(gui_corr, vals{i});
    n_d = countCorrSampleSize(gui_corr.dspn_cor, vals{i}, gui_tbl(gui_selection,:));
    n_i = countCorrSampleSize(gui_corr.ispn_cor, vals{i}, gui_tbl(gui_selection,:));
    R = [R; reportGroupCompare({'dSPN','iSPN'}, {tho_d_1, tho_i_1}, ...
            'pvals',{p_d_1, p_i_1}, 'pThresh',0.001, 'n',{n_d, n_i}, ...
            'Figure',fig, 'Panel','K', 'Measure',corr_desc.(vals{i}))]; %#ok<AGROW>
end

% plotted values (no tests in D-J: LME mean and SEM as drawn)
lick_bins = unique(lick_tbl.rewards(~isnan(lick_tbl.rewards)));      % bins used by addSingleMouseLinePlot
prog_bins = unique(pBin_tbl.progressBin(~isnan(pBin_tbl.progressBin)));
T_lick_mouse = sourceDataPerMouse(lick_mouse, lick_bins, 'xName','reward_sample');
T_lick_mouse = addvars(T_lick_mouse, T_lick_mouse.reward_sample / sr, 'After', 1, 'NewVariableNames', 'time_s');
blocks = {
    'Summary: sample sizes; K box plots, % significant, dSPN vs iSPN rank-sum', R
    'Sup2D top: reward-aligned lick count (LME mean, SEM)', ...
        sourceDataCurve(lick_trig, {'dSPN'}, 'labels',{'lick_count'}, 'xName','time_s')
    'Sup2D bottom: lick count, mean per mouse', T_lick_mouse
    'Sup2E: selected bout duration histogram (s)', T_dur_hist
    'Sup2E: duration of each selected bout', sourceDataBouts(gui_tbl_E, bout_sel, {'duration'})
    'Sup2F: bout peak velocity histogram (cm/s; cut-off 15)', T_pkvel_hist
    'Sup2F: peak velocity of each bout', sourceDataBouts(gui_tbl, bout_all, {'pkvel','rewarded'})
    'Sup2G: velocity (LME mean, SEM) by bout progress', ...
        sourceDataCurve(gui_velAcc, {'velocity'}, 'xName','progress_bin')
    'Sup2G: velocity, mean per mouse (gray lines)', ...
        sourceDataPerMouse(vel_mouse, prog_bins, 'xName','progress_bin')
    'Sup2H: acceleration (LME mean, SEM) by bout progress', ...
        sourceDataCurve(gui_velAcc, {'dSPN'}, 'labels',{'acceleration'}, 'xName','progress_bin')
    'Sup2H: acceleration, mean per mouse (gray lines)', ...
        sourceDataPerMouse(acc_mouse, prog_bins, 'xName','progress_bin')
    'Sup2I: dSPN and iSPN dF/F (LME mean, SEM) by bout progress', ...
        sourceDataCurve(gui_binned, {'dSPN','iSPN'}, 'xName','progress_bin')
    'Sup2J: onset-aligned dSPN, iSPN dF/F and velocity (LME mean, SEM)', ...
        sourceDataCurve(gui_trig_dat.onsets, {'dSPN','iSPN','velocity'}, 'xName','time_s')
    'Sup2J: offset-aligned dSPN, iSPN dF/F and velocity (LME mean, SEM)', ...
        sourceDataCurve(gui_trig_dat.offsets, {'dSPN','iSPN','velocity'}, 'xName','time_s')
    'Sup2K top: population correlation per imaging field (red = p<0.001)', ...
        [sourceDataCorrPoints(gui_corr.dspn_cor, 'binned_v_sess_cor', 'dSPN', 0.001)
         sourceDataCorrPoints(gui_corr.ispn_cor, 'binned_v_sess_cor', 'iSPN', 0.001)]
    'Sup2K bottom: single-neuron correlation (red = p<0.001)', ...
        [sourceDataCorrPoints(gui_corr.dspn_cor, 'binned_v_eachN_cor', 'dSPN', 0.001)
         sourceDataCorrPoints(gui_corr.ispn_cor, 'binned_v_eachN_cor', 'iSPN', 0.001)]
    'Sample size per mouse: D (track trials)', n_trial.perMouse
    'Sample size per mouse: E, G-J (selected GUI bouts)', n_gui.perMouse};
writeSourceData(statsFile, fig, blocks);

end

%% Helper
function [tho_d,p_d,tho_i,p_i]=  grab_celldat(corr_struct, varname)
    tho_d = cell2mat(cellfun(@(x) x.tho, {corr_struct.dspn_cor.(varname)}, 'UniformOutput', false)');
    p_d = cell2mat(cellfun(@(x) x.p, {corr_struct.dspn_cor.(varname)}, 'UniformOutput', false)');
    tho_i = cell2mat(cellfun(@(x) x.tho, {corr_struct.ispn_cor.(varname)}, 'UniformOutput', false)');
    p_i = cell2mat(cellfun(@(x) x.p, {corr_struct.ispn_cor.(varname)}, 'UniformOutput', false)');
end
