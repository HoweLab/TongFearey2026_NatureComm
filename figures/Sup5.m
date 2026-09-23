function Sup5(cfg)
% Sup 5: Supplement of main Fig 2
% A: track - spontaneous difference of the DI imbalance
% B-C: dSPN / iSPN for high vs. low acceleration groups (companion of Fig 2F-J)
% D-H: anticipatory (pre-reward) licking vs. no pre-licking trials
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure
fpath_gui = cfg.gui;
fpath_1w = cfg.track1w;

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
lick_win = [-62,31]; % reward-aligned licking window (frames)

%% start generate figures for Sup5
figDir = fullfile(save2where, 'Sup5');
if ~isfolder(figDir), mkdir(figDir); end

% A: track vs. GUI difference of DI diff (from NewAnalyses_SX/ReviewerComm.m)
vuBin_tbl = groupsummary(track_tbl(trial_selection,:), ...
    {'mouse','sess','field','trial','vuBin','isD','isI', ...
    'celltype'}, ... % grouping variables
    'mean',{'avgFc3','velocity','acceleration','info_trialTimes'});
pBin_tbl = groupsummary(gui_tbl(gui_selection,:), ...
    {'mouse','sess','field','bout','progressBin','isD','isI','celltype'}, ... % grouping variables
    'mean',{'avgFc3','velocity','acceleration'});
[A_binned, figs] = get_DIdiff_glm_SX(vuBin_tbl, pBin_tbl);
save_img(figs.figExp01,figDir,'Sup5_A_top_1w_gui_DIdiff_Sigtest'); close all
[A_onsets, figs] = get_DIdiff_glm_SX(track_tbl(bout_selection & trial_selection,:), gui_tbl(gui_selection,:), ...
    'fake_info', struct('startLoc',-1,'rewardLoc',3),...
    'xdata',on_win,'y','avgFc3','field','onsets','vel_y',[]);
save_img(figs.figExp01,figDir,'Sup5_A_bottom_1w_gui_onsets_DIdiff_Sigtest'); close all
[A_offsets, figs] = get_DIdiff_glm_SX(track_tbl(bout_selection & trial_selection,:), gui_tbl(gui_selection,:), ...
    'fake_info', struct('startLoc',-3,'rewardLoc',1),...
    'xdata',off_win,'y','avgFc3','field','offsets','vel_y',[]);
save_img(figs.figExp01,figDir,'Sup5_A_bottom_1w_gui_offsets_DIdiff_Sigtest'); close all

% B, C: dSPN / iSPN per acceleration tertile (same pipeline as Figure 2F-K;
% from NewAnalyses_SX/paperFigure_GUI_1w_AccVel_Corr.m); the per-tertile
% dSPN/iSPN figures are saved inside plot_AccRange_DI
win_events = {'onsets','offsets'};
stats_vars = {'acc_peak','dec_peak'};      % onset: peak accel; offset: peak decel
exp_name = {'Sup5_BC_1w','Sup5_BC_gui'};
accRange = select_VelAccRange_tbl({track_tbl(trial_selection & bout_selection,:), ...
    gui_tbl(gui_selection,:)}, 'numRange', 3, ...
    'BinOrTrig', 'Trig', 'win_events', win_events, ...
    'wind', {-1*sr:2*sr, -2*sr:sr});
trig_dat_nosig = plot_AccRange_DI(accRange, figDir, 'y', 'avgFc3', 'vel_y', 'acceleration', ...
    'bin_unit', win_events, 'wind_name', win_events, 'binName', 'Trig', ...
    'stats_vars', stats_vars, 'exp_name', exp_name);
close all

% D-H: pre-lick vs. no pre-lick trials (from NewAnalyses_SX/ReviewerComm.m)
[prelick_bin, no_prelick_bin, lickfigs] = find_lick_trialTypes(track_tbl, trial_selection,'avgFc3',lick_win, sr);
save_img(lickfigs.fig_lick_trig, figDir, 'Sup5_D_top_1w_lick_counts_NoPrelickVsPrelick_rewardTrig')
save_img(lickfigs.fig_lick_sp, figDir, 'Sup5_D_bottom_1w_lick_counts_NoPrelickVsPrelick_vuBin')
close all
[figs_lick, ~, lick_binned] = plot_subtable({no_prelick_bin, prelick_bin},{'mean_avgFc3'}, ...
    {'1w noprelick';'1w prelick'});
match_yaxis_across_figures(figs_lick)
save_img(figs_lick{1},figDir,'Sup5_EF_1w_trackbin_NoPrelickTrials')
save_img(figs_lick{2},figDir,'Sup5_GH_1w_trackbin_PrelickTrials')
close all

%% Count sizes & report stats -> Excel tab 'Sup5'
fig = 'Sup5';
n_trial = countSampleSize(track_tbl(trial_selection,:), 'trial');                 % A top
n_bout  = countSampleSize(track_tbl(trial_selection & bout_selection,:), 'bout'); % A bottom
n_gui   = countSampleSize(gui_tbl(gui_selection,:), 'bout');                      % A
n_nopre = countSampleSize(no_prelick_bin, 'trial');                               % D-F
n_pre   = countSampleSize(prelick_bin, 'trial');                                  % D, G-H

R = table();
ss = {'A top', 'track traversal, position-binned', n_trial;
      'A bottom', 'in-track locomotion bouts', n_bout;
      'A', 'spontaneous locomotion (GUI) bouts', n_gui;
      'D-F', 'trials without pre-reward licking', n_nopre;
      'D, G-H', 'trials with pre-reward licking', n_pre};
for s = 1:size(ss,1)
    for ct = {'dSPN','iSPN'}
        R = [R; statsRow('Figure',fig,'Panel',ss{s,1},'Measure',['sample size (' ss{s,2} ')'], ...
                'Group',ct{1},'nStruct',ss{s,3})]; %#ok<AGROW>
    end
end
% B, C: tertile cut-offs and sample size of each plotted tertile
bc = {'B','C'}; exp_label = {'track','GUI'};
lev = struct('name',{'low','high'}, 'id',{1,3});
for j = 1:numel(stats_vars)
    st = accRange.trial_stats{j};
    cut = st.rangeVal(:, strcmp(st.statsName, stats_vars{j}));
    R = [R; statsRow('Figure',fig, 'Panel',bc{j}, ...
            'Measure',sprintf('%s tertile cut-offs (from track, applied to both)', stats_vars{j}), ...
            'Note',sprintf(['low <= %.4g < mid <= %.4g < high <= %.4g cm/s^2; ' ...
            'GUI %s above the track maximum are excluded (level 4)'], cut, win_events{j}))]; %#ok<AGROW>
    for i = 1:numel(exp_name)
        tbl_w = accRange.bin_tbl{i,j};
        for L = 1:numel(lev)
            nL = countSampleSize(tbl_w(tbl_w.(stats_vars{j}) == lev(L).id,:), 'bout');
            R = [R; statsRow('Figure',fig, 'Panel',bc{j}, ...
                    'Measure',sprintf('sample size (%s %s, %s tertile)', exp_label{i}, win_events{j}, lev(L).name), ...
                    'Group',sprintf('%s %s', lev(L).name, stats_vars{j}), 'nStruct',nL, ...
                    'Note',sprintf('dSPN = %d, iSPN = %d', nL.nDSPN, nL.nISPN))]; %#ok<AGROW>
        end
    end
end

%%% plotted values %%%
TvsS = {'diff_DI_exp01'}; TvsS_lbl = {'track_minus_spont_DIdiff'};
blocks = {
    'Summary: sample sizes; B-C tertile cut-offs', R
    'Sup5A top: (dSPN - iSPN) track - spontaneous, by position/progress bin (LME, SEM, Holm p)', ...
        sourceDataCurve(A_binned, TvsS, 'labels',TvsS_lbl, 'pval',true, 'xName','bin')
    'Sup5A bottom: (dSPN - iSPN) track - spontaneous, onset-aligned (LME, SEM, Holm p)', ...
        sourceDataCurve(A_onsets, TvsS, 'labels',TvsS_lbl, 'pval',true, 'xName','time_s')
    'Sup5A bottom: (dSPN - iSPN) track - spontaneous, offset-aligned (LME, SEM, Holm p)', ...
        sourceDataCurve(A_offsets, TvsS, 'labels',TvsS_lbl, 'pval',true, 'xName','time_s')};
for j = 1:numel(stats_vars)
    for i = 1:numel(exp_name)
        x = accRange.bin_idx{j} / sr;
        T = table(x(:), 'VariableNames', {'time_s'});
        for L = [2 1] % high row, then low row
            Tl = sourceDataCurve(trig_dat_nosig{i,j,lev(L).id}{1,1}, {'dSPN','iSPN'}, ...
                'labels',strcat(lev(L).name, '_', {'dSPN','iSPN'}), 'x',x, 'xName','time_s');
            T = [T, Tl(:,2:end)]; %#ok<AGROW>
        end
        blocks(end+1,:) = {sprintf('Sup5%s: %s %s-aligned dSPN, iSPN per %s tertile (LME mean, SEM)', ...
            bc{j}, exp_label{i}, win_events{j}, stats_vars{j}), T};
    end
end
blocks = [blocks
    {'Sup5D top: reward-aligned lick count (mean, SEM across trials)', lickfigs.src_trig
     'Sup5D bottom: lick count by track position (mean, SEM across trials)', lickfigs.src_sp
     'Sup5E: no pre-lick trials, dSPN and iSPN dF/F (LME mean, SEM)', ...
        sourceDataCurve(lick_binned{1,1}, {'dSPN','iSPN'}, 'xName','position_bin')
     'Sup5F: no pre-lick trials, dSPN - iSPN (LME mean, SEM, Holm p)', ...
        sourceDataCurve(lick_binned{1,1}, {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, 'pval',true, 'xName','position_bin')
     'Sup5G: pre-lick trials, dSPN and iSPN dF/F (LME mean, SEM)', ...
        sourceDataCurve(lick_binned{2,1}, {'dSPN','iSPN'}, 'xName','position_bin')
     'Sup5H: pre-lick trials, dSPN - iSPN (LME mean, SEM, Holm p)', ...
        sourceDataCurve(lick_binned{2,1}, {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, 'pval',true, 'xName','position_bin')
     'Sample size per mouse: A top (track trials)', n_trial.perMouse
     'Sample size per mouse: A (GUI bouts)', n_gui.perMouse
     'Sample size per mouse: D-F (no pre-lick trials)', n_nopre.perMouse
     'Sample size per mouse: D, G-H (pre-lick trials)', n_pre.perMouse}];
writeSourceData(statsFile, fig, blocks);

end

%% Helper
function [prelick_bin, no_prelick_bin, figs] = find_lick_trialTypes(track_tbl, trial_selection, var_fc, lick_win, fr)
    % copied from NewAnalyses_SX/ReviewerComm.m; changes: numFc3 added to the
    % grouping variables (constant per field, for sample-size counts) and the
    % plotted lick values returned in figs.src_sp / figs.src_trig
    % var_fc  = 'avgFc3' for 1w and gui for instance
    if ischar(var_fc) | isstring(var_fc)
        var_fc = {var_fc};
    end
    resp_vars = [var_fc(:).' {'velocity','acceleration','lick_count'}];
    [allidx, alltype] = findgroups(track_tbl(:,{'mouse','fname','field','trial','sess'}));
    pre_lick = track_tbl(trial_selection & track_tbl.phase ==2 & track_tbl.reward==0 ...
            & track_tbl.y_pos <92 & track_tbl.lick &~isnan(track_tbl.vuBin),:);
    [~,presess] = findgroups(pre_lick(:,{'mouse','fname','field','trial','sess'})); % find sessions that have pre-lick
    prelick_idx = ismember(allidx, find(ismember(alltype, presess, 'rows')));
    lick_tbl =  track_tbl(prelick_idx, : );
    prelick_bin = groupsummary(lick_tbl, {'mouse','sess','field','trial','vuBin','isD','isI','celltype','numFc3'}, ...
                        'mean', resp_vars);

    % no prelicking
    no_pre = track_tbl(trial_selection & ~prelick_idx & track_tbl.phase == 3 ...
                        & track_tbl.lick & track_tbl.rewards >0,:);
    [~,no_presess] = findgroups(no_pre(:,{'mouse','fname','field','trial','sess'})); % find sessions that have pre-lick
    idx_match = ismember(allidx, find(ismember(alltype, no_presess, 'rows')));
    no_pre_tbl = track_tbl(idx_match,:);
    no_prelick_bin = groupsummary(no_pre_tbl, {'mouse','sess','field','trial','vuBin','isD','isI','celltype','numFc3'}, ...
                        'mean',resp_vars);
    % no intersection test
    [~,t1] = findgroups(lick_tbl(:,{'mouse','fname','field','trial','sess'})); % find sessions that have pre-lick
    [~,t2] = findgroups(no_pre_tbl(:,{'mouse','fname','field','trial','sess'})); % find sessions that have pre-lick
    if ~isempty(find(ismember(t1, t2, 'rows'), 1))
        error('Overlap between trials selected for pre-lick and no-prelick trials')
    end
    fprintf ('Prelick: n trials = %d; n sess = %d \n', height(t1), numel(unique(t1.fname)))
    fprintf ('No prelick: n trials = %d; n sess = %d \n', height(t2), numel(unique(t2.fname)))

    %% plot real prelick counts binned by position
    figs.fig_lick_sp = figure('Position',[900 100 900 700]);
    pre_count = groupsummary(prelick_bin, {'celltype','vuBin','isD','isI'},{'mean','std'},{'mean_lick_count'});
    nopre_count = groupsummary(no_prelick_bin, {'celltype','vuBin','isD','isI'},{'mean','std'},{'mean_lick_count'});
    isD = pre_count.isD == 1 & ~isnan(pre_count.vuBin);
    plot_error(gca,1:100, pre_count.mean_mean_lick_count(isD),...
                 pre_count.std_mean_lick_count(isD)./sqrt(pre_count.GroupCount(isD)))
    figs.src_sp = table((1:100)', pre_count.mean_mean_lick_count(isD), ...
        pre_count.std_mean_lick_count(isD)./sqrt(pre_count.GroupCount(isD)), ...
        'VariableNames', {'position_bin','prelick_mean','prelick_SEM'});
    isD = nopre_count.isD == 1 & ~isnan(nopre_count.vuBin);
    plot_error(gca,1:100, nopre_count.mean_mean_lick_count(isD),...
                 nopre_count.std_mean_lick_count(isD)./sqrt(nopre_count.GroupCount(isD)))
    figs.src_sp.noprelick_mean = nopre_count.mean_mean_lick_count(isD);
    figs.src_sp.noprelick_SEM = nopre_count.std_mean_lick_count(isD)./sqrt(nopre_count.GroupCount(isD));
    ylabel('Lick counts'); xlabel('virmen bins')
    title('Lick counts for prelick vs. noPrelick trials binned by trak pos')

    % plot triggered licking window
    figs.fig_lick_trig = figure('Position',[900 100 900 700]);
    pre_count = groupsummary(lick_tbl(lick_tbl.isD == 1 & lick_tbl.rewards <= lick_win(2) & lick_tbl.rewards >=lick_win(1),:), ...
        {'celltype','rewards','isD','isI'},{'mean','std'},{'lick_count'});
    nopre_count = groupsummary(no_pre_tbl(no_pre_tbl.isD == 1 & no_pre_tbl.rewards <= lick_win(2) & no_pre_tbl.rewards >=lick_win(1),:), ...
        {'celltype','rewards','isD','isI'},{'mean','std'},{'lick_count'});
    plot_error(gca,(lick_win(1):lick_win(2))./fr, pre_count.mean_lick_count,...
                 pre_count.std_lick_count./sqrt(pre_count.GroupCount))
    plot_error(gca,(lick_win(1):lick_win(2))./fr , nopre_count.mean_lick_count,...
                 nopre_count.std_lick_count./sqrt(nopre_count.GroupCount))
    figs.src_trig = table((lick_win(1):lick_win(2))'./fr, ...
        pre_count.mean_lick_count, pre_count.std_lick_count./sqrt(pre_count.GroupCount), ...
        nopre_count.mean_lick_count, nopre_count.std_lick_count./sqrt(nopre_count.GroupCount), ...
        'VariableNames', {'time_s','prelick_mean','prelick_SEM','noprelick_mean','noprelick_SEM'});
    ylabel('Lick counts'); xlabel('time (s)')
    title('Reward triggered average')
end
