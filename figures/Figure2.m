function Figure2(cfg)
% Main Figure2 - DI imbalance for GUI and track & kinematic matching
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure
% processed data for each session
fpath = cfg.track1w;
fpath_gui = cfg.gui;
% organize data
conditions = {'exp_stage',{'1 world'};
    'objective','20x';
    'sess_type','track';
    'mismatch',0};
track_tbl = mkTblForPlot(fpath,conditions);
conditions = {'sess_type','gui';
    'mismatch',0};
gui_tbl = mkTblForPlot(fpath_gui,conditions);

% data selection
mice_include = {'tdt','TdTG7', 'a2a'}; % all 3 batches
trial_selection = ~track_tbl.skip_trial & track_tbl.exp_stage == '1 world' ...
                & contains(cellstr(track_tbl.mouse),mice_include);
bout_selection = track_tbl.info_bout_intrack == 1 & contains(cellstr(track_tbl.mouse),mice_include);
gui_selection = gui_tbl.info_bout_rewarded ~= 1  &... %
    gui_tbl.info_bout_length > 0 & gui_tbl.info_peak_vel > 15 & ...
    contains(cellstr(gui_tbl.mouse),mice_include); % unrewarded bouts

xdata = [1:100; linspace(1,100,100)];
fake_info = struct('startLoc',1,'rewardLoc',100);
% Trig window
sr = 31;
on_win = [-31:3*31; linspace(-1,3,125)];
off_win = [3*-31:31; linspace(-3,1,125)];

%% start generate figures for fig2
figDir = fullfile(save2where, 'Fig2');
if ~isfolder(figDir), mkdir(figDir); end

vuBin_tbl = groupsummary(track_tbl(trial_selection,:), ... % track
    {'mouse','sess','field','trial','vuBin','isD','isI', ...
    'celltype','numFc3'}, ...
    'mean',{'numFc3','avgFc3','velocity','acceleration','info_trialTimes', 'lick_count'});
pBin_tbl = groupsummary(gui_tbl(gui_selection,:), ...      % gui
    {'mouse','sess','field','bout','progressBin','isD','isI','celltype','numFc3'}, ...
    'mean',{'avgFc3','velocity','acceleration'}); % calculate mean for avgFc3 and velocity
vu_binned = getPlotDatFromTbl(vuBin_tbl,'vuBin','y','mean_avgFc3', ...
    'vel_y','mean_acceleration','xdata',xdata);
gui_binned = getPlotDatFromTbl(pBin_tbl,'progressBin','y','mean_avgFc3', ...
    'vel_y','mean_acceleration','xdata',xdata);
vu_onsets = getPlotDatFromTbl(track_tbl(trial_selection & bout_selection,:), ...
    'onsets','xdata',on_win,'vel_y','acceleration');
vu_offsets = getPlotDatFromTbl(track_tbl(trial_selection & bout_selection,:), ...
    'offsets','xdata',off_win,'vel_y','acceleration');
vu_trig_dat = struct('onsets',vu_onsets, 'offsets', vu_offsets);
gui_onsets = getPlotDatFromTbl(gui_tbl(gui_selection,:), ...
    'onsets','xdata',on_win,'vel_y','acceleration');
gui_offsets = getPlotDatFromTbl(gui_tbl(gui_selection,:), ...
    'offsets','xdata',off_win,'vel_y','acceleration');
gui_trig_dat = struct('onsets',gui_onsets, 'offsets', gui_offsets);

%%% Plotting %%%
% A: 1w DI
figure('Position',[900 100 900 700])
plotBinned_new(vu_binned,fake_info, ...
    'fromLME',true,'traces',{'diff_di'}, ...
    'plot_vel',true,'plot_sig',true)
save_img(gcf,figDir,'Fig2_A_1world_positionBins_dff_accl_aligned')

% C: GUI DI
figure('Position',[900 100 900 700])
plotBinned_new(gui_binned,fake_info, ...
    'fromLME',true,'traces',{'diff_di'}, ...
    'plot_vel',true,'plot_sig',true)
save_img(gcf,figDir,'Fig2_C_gui_positionBins_dff_accl_aligned')

% B: 1w trig (significance bars as in the paper; p reported in source data)
figure('Position',[100 100 900 700])
plotTriggered_new(vu_trig_dat,'traces',{'diff_di'},'fromLME',true, ...
    'plot_vel',true,'plot_sig',true, 'trig_events',fieldnames(vu_trig_dat));
save_img(gcf,figDir,'Fig2_B_1w_trig_dff_accl_aligned')

% D: Gui trig
figure('Position',[100 100 900 700])
plotTriggered_new(gui_trig_dat,'traces',{'diff_di'},'fromLME',true, ...
    'plot_vel',true,'plot_sig',true, 'trig_events',fieldnames(vu_trig_dat));
save_img(gcf,figDir,'Fig2_D_gui_trig_dff_accl_aligned')
close all

% E: Cross correlation of LME mean DI diff vs. acceleration/velocity (track only)
% (from NewAnalyses_SX/paperFigure_GUI_1w_AccVel_Corr.m)
% top: position bins
[sessout_1w, avgout_1w] = crosscor_DI_behav(track_tbl(trial_selection,:), ...
    'onlyAvg', 1, 'whichFc3', 'avgFc3');
plot_crosscor(avgout_1w, sessout_1w)
save_img(gcf, figDir, 'Fig2_E_top_1w_LME_vel_acc_DI_crosscor'); close(gcf)
% middle/bottom: onset/offset triggered
[sessout_1w_onset, avgout_1w_onset] = crosscor_DI_behav(track_tbl(trial_selection & bout_selection,:), ...
    'x_data', on_win, 'bin_unit', 'onsets', 'BinOrTrig', 'trig', 'onlyAvg', 1);
[sessout_1w_offset, avgout_1w_offset] = crosscor_DI_behav(track_tbl(trial_selection & bout_selection,:), ...
    'x_data', off_win, 'bin_unit', 'offsets', 'BinOrTrig', 'trig', 'onlyAvg', 1);
plot_crosscor(avgout_1w_onset, sessout_1w_onset)
save_img(gcf, figDir, 'Fig2_E_middle_1w_onset_LME_vel_acc_DI_crosscor'); close(gcf)
plot_crosscor(avgout_1w_offset, sessout_1w_offset)
save_img(gcf, figDir, 'Fig2_E_bottom_1w_offset_LME_vel_acc_DI_crosscor'); close(gcf)

% F-K: DI at high vs. low onset acceleration / offset deceleration (tertiles)
% (from NewAnalyses_SX/paperFigure_GUI_1w_AccVel_Corr.m; selection as in
%  paperFigScripts/Count_exp_fig_SampleSize.m)
win_events = {'onsets','offsets'};
stats_vars = {'acc_peak','dec_peak'};      % onset: peak accel; offset: peak decel
exp_name = {'Fig2_1w','Fig2_gui'};         % used by plot_AccRange_DI for file names
accRange = select_VelAccRange_tbl({track_tbl(trial_selection & bout_selection,:), ...
    gui_tbl(gui_selection,:)}, 'numRange', 3, ...
    'BinOrTrig', 'Trig', 'win_events', win_events, ...
    'wind', {-1*sr:2*sr, -2*sr:sr});
% H, K: acceleration per tertile (+ DI per tertile); saved inside
trig_dat_nosig = plot_AccRange_DI(accRange, figDir, 'y', 'avgFc3', 'vel_y', 'acceleration', ...
    'bin_unit', win_events, 'wind_name', win_events, 'binName', 'Trig', ...
    'stats_vars', stats_vars, 'exp_name', exp_name);
% significance for high (3) vs. low (1) within each experiment
trig_dat_sig = plot_AccRange_DI(accRange, figDir, 'y', 'avgFc3', 'vel_y', 'acceleration', ...
    'bin_unit', win_events, 'wind_name', win_events, 'binName', 'Trig', ...
    'stats_vars', stats_vars, 'ifgetsig', 1, 'levlplotted', [1,3], 'ComparisonTp', 'bt_level');

% F (1w onset), G (gui onset), I (1w offset), J (gui offset)
fgij = {'F','I'; 'G','J'};                 % rows: 1w, gui; cols: onset, offset
field2plot = {'iSPN', 'dSPN','diff_di'};
sigfield = {'diff_I_exp01', 'diff_D_exp01','diff_DI_exp01'};
for i = 1:numel(exp_name)
    for j = 1:numel(stats_vars)
        plot_dat = struct();
        plot_dat.low = trig_dat_nosig{i,j,1}{1,1};
        plot_dat.high = trig_dat_nosig{i,j,3}{1,1};
        sig = trig_dat_sig{i,j};
        titleN = [exp_name{i}(1:5), fgij{i,j}, exp_name{i}(5:end), '_', stats_vars{j}, '_high_lowaccl_'];
        plot_SigOnTrace(plot_dat, sig, field2plot, sigfield, sig.bins, titleN, figDir, sr)
    end
end
close all

%% Count sizes & report stats -> Excel tab 'Fig2'
fig = 'Fig2';
n_trial = countSampleSize(track_tbl(trial_selection,:), 'trial');                 % A, E top
n_bout  = countSampleSize(track_tbl(trial_selection & bout_selection,:), 'bout'); % B, E middle/bottom
n_gui   = countSampleSize(gui_tbl(gui_selection,:), 'bout');                      % C, D

R = table();
ss = {'A, E top', 'track traversal, position-binned', n_trial;
      'B, E middle-bottom', 'in-track locomotion bouts', n_bout;
      'C-D', 'spontaneous locomotion (GUI) bouts', n_gui};
for s = 1:size(ss,1)
    for ct = {'dSPN','iSPN'}
        R = [R; statsRow('Figure',fig,'Panel',ss{s,1},'Measure',['sample size (' ss{s,2} ')'], ...
                'Group',ct{1},'nStruct',ss{s,3})]; %#ok<AGROW>
    end
end

% E: cross-correlation peak (marked in the plots; no test shown)
R = [R; reportCrossCorr(avgout_1w, 'Figure',fig, 'Panel','E top', ...
        'Measure','position-binned', 'n',n_trial, 'lagUnit','bin')];
R = [R; reportCrossCorr(avgout_1w_onset, 'Figure',fig, 'Panel','E middle', ...
        'Measure','onset-aligned', 'n',n_bout, 'lagScale',1/sr, 'lagUnit','s')];
R = [R; reportCrossCorr(avgout_1w_offset, 'Figure',fig, 'Panel','E bottom', ...
        'Measure','offset-aligned', 'n',n_bout, 'lagScale',1/sr, 'lagUnit','s')];

% F-K: tertile cut-offs and sample size of each tertile
hk = {'H','K'};                            % acceleration profiles: onset, offset
exp_label = {'track','GUI'};
lev = struct('name',{'low','mid','high'}, 'id',{1,2,3});
for j = 1:numel(stats_vars)
    st = accRange.trial_stats{j};
    cut = st.rangeVal(:, strcmp(st.statsName, stats_vars{j}));
    R = [R; statsRow('Figure',fig, 'Panel',[fgij{1,j} ',' fgij{2,j} ',' hk{j}], ...
            'Measure',sprintf('%s tertile cut-offs (from track, applied to both)', stats_vars{j}), ...
            'Note',sprintf(['low <= %.4g < mid <= %.4g < high <= %.4g cm/s^2; ' ...
            'GUI %s above the track maximum are excluded (level 4)'], cut, win_events{j}))]; %#ok<AGROW>
    for i = 1:numel(exp_name)
        tbl_w = accRange.bin_tbl{i,j};
        for L = 1:numel(lev)
            nL = countSampleSize(tbl_w(tbl_w.(stats_vars{j}) == lev(L).id,:), 'bout');
            R = [R; statsRow('Figure',fig, 'Panel',[fgij{i,j} ',' hk{j}], ...
                    'Measure',sprintf('sample size (%s %s, %s tertile)', exp_label{i}, win_events{j}, lev(L).name), ...
                    'Group',sprintf('%s %s', lev(L).name, stats_vars{j}), 'nStruct',nL, ...
                    'Note',sprintf('dSPN = %d, iSPN = %d', nL.nDSPN, nL.nISPN))]; %#ok<AGROW>
        end
    end
end

%%% plotted values %%%
% A-D: dSPN - iSPN (LME mean, SEM, Holm p for the significance bars) and acceleration (gray)
DI = {'diff_di'}; DI_lbl = {'dSPN_minus_iSPN'};
curve_DI_acc = @(d, xName) [sourceDataCurve(d, DI, 'labels',DI_lbl, 'pval',true, 'xName',xName), ...
    removevars(sourceDataCurve(d, {'velocity'}, 'labels',{'acceleration'}, 'xName',xName), xName)];
blocks = {
    'Summary: sample sizes; E cross-correlation peaks; F-K tertile cut-offs', R
    'Fig2A: track, by track position', curve_DI_acc(vu_binned, 'position_bin')
    'Fig2B: track, onset-aligned', curve_DI_acc(vu_trig_dat.onsets, 'time_s')
    'Fig2B: track, offset-aligned', curve_DI_acc(vu_trig_dat.offsets, 'time_s')
    'Fig2C: GUI, by bout progress', curve_DI_acc(gui_binned, 'progress_bin')
    'Fig2D: GUI, onset-aligned', curve_DI_acc(gui_trig_dat.onsets, 'time_s')
    'Fig2D: GUI, offset-aligned', curve_DI_acc(gui_trig_dat.offsets, 'time_s')
    'Fig2E top: cross-correlation, position lag', sourceDataXcorr(avgout_1w, 1, 'lag_bin')
    'Fig2E middle: cross-correlation, onset-aligned lag', sourceDataXcorr(avgout_1w_onset, 1/sr, 'lag_s')
    'Fig2E bottom: cross-correlation, offset-aligned lag', sourceDataXcorr(avgout_1w_offset, 1/sr, 'lag_s')};

% F, G, I, J: low vs high tertile (as plot_SigOnTrace), p = high vs low (gray bars)
fld_lbl = {'iSPN','dSPN','dSPN_minus_iSPN'};
for j = 1:numel(stats_vars)
    for i = 1:numel(exp_name)
        sig = trig_dat_sig{i,j};
        x = sig.bins / sr;
        Tl = sourceDataCurve(trig_dat_nosig{i,j,1}{1,1}, field2plot, ...
            'labels',strcat('low_', fld_lbl), 'x',x, 'xName','time_s');
        Th = sourceDataCurve(trig_dat_nosig{i,j,3}{1,1}, field2plot, ...
            'labels',strcat('high_', fld_lbl), 'x',x, 'xName','time_s');
        Ts = sourceDataCurve(sig, sigfield, 'labels',strcat('high_vs_low_', fld_lbl), ...
            'pval',true, 'x',x, 'xName','time_s');
        Ts = Ts(:, endsWith(Ts.Properties.VariableNames, '_p_Holm'));
        blocks(end+1,:) = {sprintf('Fig2%s: %s %s, low vs high %s tertile', ...
            fgij{i,j}, exp_label{i}, win_events{j}, stats_vars{j}), [Tl, Th(:,2:end), Ts]};
    end
end

% H, K: acceleration of each tertile (as plot_AccRange_DI)
for j = 1:numel(stats_vars)
    for i = 1:numel(exp_name)
        x = accRange.bin_idx{j} / sr;
        T = table(x(:), 'VariableNames', {'time_s'});
        for L = 1:numel(lev)
            Ta = sourceDataCurve(trig_dat_nosig{i,j,lev(L).id}{1,1}, {'velocity'}, ...
                'labels',{[lev(L).name '_acceleration']}, 'x',x, 'xName','time_s');
            T = [T, Ta(:,2:end)]; %#ok<AGROW>
        end
        blocks(end+1,:) = {sprintf('Fig2%s: %s %s-aligned acceleration per %s tertile', ...
            hk{j}, exp_label{i}, win_events{j}, stats_vars{j}), T};
    end
end

blocks = [blocks
    {'Sample size per mouse: A (track trials)', n_trial.perMouse
     'Sample size per mouse: B, E middle-bottom (track bouts)', n_bout.perMouse
     'Sample size per mouse: C-D (GUI bouts)', n_gui.perMouse}];
writeSourceData(statsFile, fig, blocks);

end
