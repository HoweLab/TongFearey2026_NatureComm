function Sup3(cfg)
% Sup 3: Supplement of main Fig 2
% Validation for DI imbalance using A2A mice (A-H) and fraction of iSPNs
% per cohort (I)
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure
fpath_gui = cfg.gui;
fpath_1w = cfg.track1w;
roipath_track = cfg.roiTrack;
roipath_gui   = cfg.roiGui;

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

% A-H: A2A mice only
a2a = {'a2a'};
trial_sel_a2a = trial_selection & contains(cellstr(track_tbl.mouse),a2a);
bout_sel_a2a  = bout_selection  & contains(cellstr(track_tbl.mouse),a2a);
gui_sel_a2a   = gui_selection   & contains(cellstr(gui_tbl.mouse),a2a);
% I: cohorts
cohorts = {'D1-tdTomato', {'tdt','TdTG7'};
           'A2A-cre',     {'a2a'}};

xdata = [1:100; linspace(1,100,100)];
fake_info = struct('startLoc',1,'rewardLoc',100);

% Trig window
sr = 31;
on_win = [-31:3*31; linspace(-1,3,125)];
off_win = [3*-31:31; linspace(-3,1,125)];

%% start generate figures for Sup3
figDir = fullfile(save2where, 'Sup3');
if ~isfolder(figDir), mkdir(figDir); end

% (from paperFigScripts/paperFigurePopulation_1w_gui.m, mice_include = {'a2a'})
vuBin_tbl = groupsummary(track_tbl(trial_sel_a2a,:), ... % track
    {'mouse','sess','field','trial','vuBin','isD','isI', ...
    'celltype','numFc3'}, ...
    'mean',{'numFc3','avgFc3','velocity','acceleration'});
pBin_tbl = groupsummary(gui_tbl(gui_sel_a2a,:), ...      % gui
    {'mouse','sess','field','bout','progressBin','isD','isI','celltype','numFc3'}, ...
    'mean',{'avgFc3','velocity','acceleration'});
exps = struct();
exps.track.binned_velAcc = getPlotDatFromTbl(vuBin_tbl,'vuBin','y','mean_acceleration', ...
    'vel_y','mean_velocity','xdata',xdata);
exps.track.binned_data = getPlotDatFromTbl(vuBin_tbl,'vuBin','y','mean_avgFc3', ...
    'vel_y','mean_acceleration','xdata',xdata);
exps.gui.binned_velAcc = getPlotDatFromTbl(pBin_tbl,'progressBin','y','mean_acceleration', ...
    'vel_y','mean_velocity','xdata',xdata);
exps.gui.binned_data = getPlotDatFromTbl(pBin_tbl,'progressBin','y','mean_avgFc3', ...
    'vel_y','mean_acceleration','xdata',xdata);
trig_tbls = struct('track', track_tbl(trial_sel_a2a & bout_sel_a2a,:), ...
                   'gui',   gui_tbl(gui_sel_a2a,:));
for e = fieldnames(trig_tbls)'
    t = trig_tbls.(e{1});
    exps.(e{1}).trig_dat = struct( ...          % dSPN, iSPN, diff_di, velocity
        'onsets',  getPlotDatFromTbl(t,'onsets','xdata',on_win,'vel_y','velocity'), ...
        'offsets', getPlotDatFromTbl(t,'offsets','xdata',off_win,'vel_y','velocity'));
    exps.(e{1}).trig_acc = struct( ...          % acceleration only (stored as 'velocity')
        'onsets',  getPlotDatFromTbl(t,'onsets','xdata',on_win,'y',[],'vel_y','acceleration'), ...
        'offsets', getPlotDatFromTbl(t,'offsets','xdata',off_win,'y',[],'vel_y','acceleration'));
end
clear trig_tbls t

%%% Plotting %%%
% A-D: track traversal; E-H: spontaneous locomotion (GUI)
pan = struct('track',{{'A','B','C','D'}}, 'gui',{{'E','F','G','H'}});
xlbl = struct('track','positionBins', 'gui','boutProgressBin');
for e = {'track','gui'}
    ex = exps.(e{1}); p = pan.(e{1});
    % A / E: dSPN, iSPN, velocity, acceleration binned
    for tr = {'dSPN','iSPN'}
        figure('Position',[900 100 900 700])
        plotBinned_new(ex.binned_data,fake_info, ...
            'fromLME',true,'traces',tr,'plot_vel',false)
        save_img(gcf,figDir,sprintf('Sup3_%s_%s_%s_%s',p{1},e{1},xlbl.(e{1}),tr{1})); close(gcf)
    end
    figure('Position',[900 100 900 700])
    plotBinned_new(ex.binned_velAcc,fake_info, ...
        'fromLME',true,'traces',{'velocity'}, 'line_colors',{[0 0 1]})
    save_img(gcf,figDir,sprintf('Sup3_%s_%s_%s_velocity',p{1},e{1},xlbl.(e{1}))); close(gcf)
    figure('Position',[900 100 900 700])
    plotBinned_new(ex.binned_velAcc,fake_info, ...
        'fromLME',true,'traces',{'dSPN'}, 'line_colors',{[0 0 1]})
    save_img(gcf,figDir,sprintf('Sup3_%s_%s_%s_acceleration',p{1},e{1},xlbl.(e{1}))); close(gcf)

    % B / F: DI diff binned, with significance
    figure('Position',[900 100 900 700])
    plotBinned_new(ex.binned_data,fake_info, ...
        'fromLME',true,'traces',{'diff_di'}, 'plot_sig',true, ...
        'line_colors',{[0 0 0]},'plot_vel',false);
    save_img(gcf,figDir,sprintf('Sup3_%s_%s_%s_DIdiff',p{2},e{1},xlbl.(e{1}))); close(gcf)

    % C / G: triggered dSPN, iSPN, velocity, acceleration
    for tr = {'dSPN','iSPN'}
        figure('Position',[100 100 900 700])
        plotTriggered_new(ex.trig_dat,'traces',tr,'fromLME',true, ...
            'plot_vel',false, 'trig_events',fieldnames(ex.trig_dat));
        save_img(gcf,figDir,sprintf('Sup3_%s_%s_trig_%s',p{3},e{1},tr{1})); close(gcf)
    end
    figure('Position',[900 100 900 700])
    plotTriggered_new(ex.trig_dat,'traces',{'velocity'},'fromLME',true, ...
        'line_colors',{[0 0 1]},'plot_vel',false, 'trig_events',fieldnames(ex.trig_dat));
    save_img(gcf,figDir,sprintf('Sup3_%s_%s_trig_velocity',p{3},e{1})); close(gcf)
    figure('Position',[900 100 900 700])
    plotTriggered_new(ex.trig_acc,'traces',{'velocity'},'fromLME',true, ...
        'line_colors',{[0 0 1]},'plot_vel',false, 'trig_events',fieldnames(ex.trig_acc));
    save_img(gcf,figDir,sprintf('Sup3_%s_%s_trig_acceleration',p{3},e{1})); close(gcf)

    % D / H: triggered DI diff, with significance
    figure('Position',[100 100 900 700])
    plotTriggered_new(ex.trig_dat,'traces',{'diff_di'},'fromLME',true, ...
        'plot_sig',true,'line_colors',{[0 0 0]}, ...
        'plot_vel',false, 'trig_events',fieldnames(ex.trig_dat));
    save_img(gcf,figDir,sprintf('Sup3_%s_%s_trig_DIdiff',p{4},e{1})); close(gcf)
end

% I: fraction of iSPNs per imaging field, track vs. spontaneous, per cohort
skiplogical = 1;
[cat_track, grp_idx_track, fname_track] = find_cells_summaryROIs(track_tbl, roipath_track,'Fc3',skiplogical); % this is slow
[cat_gui, grp_idx_gui, fname_gui] = find_cells_summaryROIs(gui_tbl, roipath_gui,'Fc3',skiplogical); % this is slow
cellFrac = cell(size(cohorts,1), 2); % {cohort, track/gui} per-field tables
pI = nan(size(cohorts,1),1);
for c = 1:size(cohorts,1)
    idx_1w = contains(cellstr(track_tbl.mouse),cohorts{c,2});
    idx_gui = gui_selection & contains(cellstr(gui_tbl.mouse),cohorts{c,2});
    cellFrac{c,1} = count_cells_field(cat_track, grp_idx_track, fname_track, track_tbl.celltype, idx_1w);
    cellFrac{c,2} = count_cells_field(cat_gui, grp_idx_gui, fname_gui, gui_tbl.celltype, idx_gui);

    figure('Name',cohorts{c,1}); fnum = get(gcf,'Number');
    boxViolin(cellFrac{c,1}.fraction_iSPN, cellFrac{c,2}.fraction_iSPN, fnum);
    [pI(c),~] = ranksum(cellFrac{c,1}.fraction_iSPN', cellFrac{c,2}.fraction_iSPN');
    title(['prct of ispn ' cohorts{c,1} ': track (left); gui (right), p =',num2str(pI(c))])
    save_img(gcf, figDir, ['Sup3_I_prct_ispn_track_vs_GUI_' cohorts{c,1}])
    close gcf
end

%% Count sizes & report stats -> Excel tab 'Sup3'
fig = 'Sup3';
n_trial = countSampleSize(track_tbl(trial_sel_a2a,:), 'trial');              % A, B
n_bout  = countSampleSize(track_tbl(trial_sel_a2a & bout_sel_a2a,:), 'bout'); % C, D
n_gui   = countSampleSize(gui_tbl(gui_sel_a2a,:), 'bout');                    % E-H

R = table();
ss = {'A-B', 'A2A track traversal, position-binned', n_trial;
      'C-D', 'A2A in-track locomotion bouts', n_bout;
      'E-H', 'A2A spontaneous locomotion (GUI) bouts', n_gui};
for s = 1:size(ss,1)
    for ct = {'dSPN','iSPN'}
        R = [R; statsRow('Figure',fig,'Panel',ss{s,1},'Measure',['sample size (' ss{s,2} ')'], ...
                'Group',ct{1},'nStruct',ss{s,3})]; %#ok<AGROW>
    end
end
% I: box plots (median/IQR) and track vs. spontaneous rank-sum per cohort (the test shown)
for c = 1:size(cohorts,1)
    R = [R; reportGroupCompare({'Track','Spont'}, ...
            {cellFrac{c,1}.fraction_iSPN, cellFrac{c,2}.fraction_iSPN}, ...
            'Figure',fig, 'Panel','I', ...
            'Measure',sprintf('fraction of iSPNs per imaging field, %s', cohorts{c,1}))]; %#ok<AGROW>
end

% plotted values (A, C, E, G: no tests; B, D, F, H: Holm p for the significance bars)
blocks = {'Summary: sample sizes; I box plots and track vs spont rank-sum', R};
lbl = struct('track',{{'track position','position_bin'}}, 'gui',{{'bout progress','progress_bin'}});
for e = {'track','gui'}
    ex = exps.(e{1}); p = pan.(e{1}); xn = lbl.(e{1}){2};
    blocks = [blocks
        {sprintf('Sup3%s: A2A %s dSPN, iSPN dF/F (LME mean, SEM) by %s', p{1}, e{1}, lbl.(e{1}){1}), ...
            sourceDataCurve(ex.binned_data, {'dSPN','iSPN'}, 'xName',xn)
         sprintf('Sup3%s: A2A %s velocity and acceleration (LME mean, SEM) by %s', p{1}, e{1}, lbl.(e{1}){1}), ...
            [sourceDataCurve(ex.binned_velAcc, {'velocity'}, 'xName',xn), ...
             removevars(sourceDataCurve(ex.binned_velAcc, {'dSPN'}, 'labels',{'acceleration'}, 'xName',xn), xn)]
         sprintf('Sup3%s: A2A %s dSPN - iSPN dF/F (LME mean, SEM, Holm p) by %s', p{2}, e{1}, lbl.(e{1}){1}), ...
            sourceDataCurve(ex.binned_data, {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, 'pval',true, 'xName',xn)}]; %#ok<AGROW>
    for ev = {'onsets','offsets'}
        blocks = [blocks
            {sprintf('Sup3%s: A2A %s %s-aligned dSPN, iSPN, velocity, acceleration (LME mean, SEM)', p{3}, e{1}, ev{1}), ...
                [sourceDataCurve(ex.trig_dat.(ev{1}), {'dSPN','iSPN','velocity'}, 'xName','time_s'), ...
                 removevars(sourceDataCurve(ex.trig_acc.(ev{1}), {'velocity'}, 'labels',{'acceleration'}, 'xName','time_s'), 'time_s')]
             sprintf('Sup3%s: A2A %s %s-aligned dSPN - iSPN (LME mean, SEM, Holm p)', p{4}, e{1}, ev{1}), ...
                sourceDataCurve(ex.trig_dat.(ev{1}), {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, 'pval',true, 'xName','time_s')}]; %#ok<AGROW>
    end
end
for c = 1:size(cohorts,1)
    blocks = [blocks
        {sprintf('Sup3I: %s, cells per imaging field (track, then spont)', cohorts{c,1}), ...
            [addvars(cellFrac{c,1}, repmat("Track", height(cellFrac{c,1}), 1), 'Before', 1, 'NewVariableNames', 'condition')
             addvars(cellFrac{c,2}, repmat("Spont", height(cellFrac{c,2}), 1), 'Before', 1, 'NewVariableNames', 'condition')]}]; %#ok<AGROW>
end
blocks = [blocks
    {'Sample size per mouse: A-B (A2A track trials)', n_trial.perMouse
     'Sample size per mouse: C-D (A2A track bouts)', n_bout.perMouse
     'Sample size per mouse: E-H (A2A GUI bouts)', n_gui.perMouse}];
writeSourceData(statsFile, fig, blocks);

end

%% Helper
function T = count_cells_field(cat_dat, grp_idx, fnames, celltype, idx)
% number of dSPNs / iSPNs per imaging field (cells kept after skipLogical)
% and the fraction of iSPNs, for the rows selected by idx
    numCell = cellfun(@(x) length(x), cat_dat(idx));
    sessIdx = grp_idx(idx);
    ct = celltype(idx);
    n = max(grp_idx);
    nD = accumarray(sessIdx(ct == 'dSPN'), numCell(ct == 'dSPN'), [n 1], @mean);
    nI = accumarray(sessIdx(ct == 'iSPN'), numCell(ct == 'iSPN'), [n 1], @mean);
    used = unique(sessIdx);
    T = table(string(fnames.fname(used)), nD(used), nI(used), ...
        'VariableNames', {'field','n_dSPN','n_iSPN'});
    T.fraction_iSPN = T.n_iSPN ./ (T.n_dSPN + T.n_iSPN);
end
