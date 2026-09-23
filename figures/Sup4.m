function Sup4(cfg)
% Sup 4: Supplement of main Fig 2
% DI imbalance reproduced with deconvolved activity (A-H) and
% percent of active cells (I-P), track traversal and spontaneous locomotion
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure
% one entry per activity measure: data files, y variable, panels (track; gui)
datasets = struct( ...
    'name',   {'deconv', 'prctActive'}, ...
    'label',  {'deconvolved dF/F', '% active'}, ...
    'fpath_1w',  {cfg.track1w, cfg.prctActive1w}, ...
    'fpath_gui', {cfg.gui, cfg.prctActiveGui}, ...
    'y',      {'decFc3', 'Fc3_activeZ'}, ...
    'scale',  {30, 1}, ...  % deconv: amp/frame -> amp/second (fr = 30, as paperFigurePopulation_1w_gui.m)
    'panels', {{'A','B','C','D'; 'E','F','G','H'}, {'I','J','K','L'; 'M','N','O','P'}});

mice_include = {'tdt','TdTG7', 'a2a'}; % all 3 batches
xdata = [1:100; linspace(1,100,100)];
fake_info = struct('startLoc',1,'rewardLoc',100);

% Trig window
sr = 31;
on_win = [-31:3*31; linspace(-1,3,125)];
off_win = [3*-31:31; linspace(-3,1,125)];

%% start generate figures for Sup4
% deconvolved: paperFigScripts/paperFigurePopulation_1w_gui.m ("ALL Deconvolved plots")
% % active:    NewAnalyses_SX/ReviewerComm.m ("1w vs. GUI", Fc3_prct)
figDir = fullfile(save2where, 'Sup4');
if ~isfolder(figDir), mkdir(figDir); end

res = struct();
for d = 1:numel(datasets)
    ds = datasets(d);
    conditions = {'exp_stage',{'1 world'};
        'objective','20x';
        'sess_type','track';
        'mismatch',0};
    track_tbl = mkTblForPlot(ds.fpath_1w,conditions);
    conditions = {'sess_type','gui'; 'mismatch',0};
    gui_tbl = mkTblForPlot(ds.fpath_gui,conditions);
    track_tbl.(ds.y) = track_tbl.(ds.y) * ds.scale;
    gui_tbl.(ds.y) = gui_tbl.(ds.y) * ds.scale;

    gui_selection = gui_tbl.info_bout_rewarded ~= 1  &... %
        gui_tbl.info_bout_length > 0 & gui_tbl.info_peak_vel > 15 & ...
        contains(cellstr(gui_tbl.mouse),mice_include); % unrewarded bouts
    trial_selection = ~track_tbl.skip_trial & track_tbl.exp_stage == '1 world' ...
                    & contains(cellstr(track_tbl.mouse),mice_include);
    bout_selection = track_tbl.info_bout_intrack == 1 & contains(cellstr(track_tbl.mouse),mice_include);

    if d == 1 % sample sizes (I-P use the same mice, sessions, fields and neurons)
        n_trial = countSampleSize(track_tbl(trial_selection,:), 'trial');
        n_bout  = countSampleSize(track_tbl(trial_selection & bout_selection,:), 'bout');
        n_gui   = countSampleSize(gui_tbl(gui_selection,:), 'bout');
    end

    vuBin_tbl = groupsummary(track_tbl(trial_selection,:), ...
        {'mouse','sess','field','trial','vuBin','isD','isI','celltype'}, ...
        'mean',{ds.y});
    pBin_tbl = groupsummary(gui_tbl(gui_selection,:), ...
        {'mouse','sess','field','bout','progressBin','isD','isI','celltype'}, ...
        'mean',{ds.y});
    res.(ds.name).track.binned = getPlotDatFromTbl(vuBin_tbl,'vuBin','y',['mean_' ds.y], ...
        'vel_y',[],'xdata',xdata);
    res.(ds.name).gui.binned = getPlotDatFromTbl(pBin_tbl,'progressBin','y',['mean_' ds.y], ...
        'vel_y',[],'xdata',xdata);
    trig_tbls = struct('track', track_tbl(trial_selection & bout_selection,:), ...
                       'gui',   gui_tbl(gui_selection,:));
    for e = {'track','gui'}
        res.(ds.name).(e{1}).trig_dat = struct( ...
            'onsets',  getPlotDatFromTbl(trig_tbls.(e{1}),'onsets','xdata',on_win,'vel_y',[],'y',ds.y), ...
            'offsets', getPlotDatFromTbl(trig_tbls.(e{1}),'offsets','xdata',off_win,'vel_y',[],'y',ds.y));
    end
    clear track_tbl gui_tbl trig_tbls vuBin_tbl pBin_tbl

    %%% Plotting %%%
    xlbl = struct('track','positionBins', 'gui','boutProgressBin');
    for k = 1:2
        e = {'track','gui'}; e = e{k};
        p = ds.panels(k,:);
        r = res.(ds.name).(e);
        % A/E/I/M: dSPN, iSPN binned
        for tr = {'dSPN','iSPN'}
            figure('Position',[900 100 900 700])
            plotBinned_new(r.binned,fake_info, ...
                'fromLME',true,'traces',tr,'plot_vel',false)
            save_img(gcf,figDir,sprintf('Sup4_%s_%s_%s_%s_%s',p{1},ds.name,e,xlbl.(e),tr{1})); close(gcf)
        end
        % B/F/J/N: DI diff binned, with significance
        figure('Position',[900 100 900 700])
        plotBinned_new(r.binned,fake_info, ...
            'fromLME',true,'traces',{'diff_di'}, 'plot_sig',true, ...
            'line_colors',{[0 0 0]},'plot_vel',false);
        save_img(gcf,figDir,sprintf('Sup4_%s_%s_%s_%s_DIdiff',p{2},ds.name,e,xlbl.(e))); close(gcf)
        % C/G/K/O: triggered dSPN, iSPN
        for tr = {'dSPN','iSPN'}
            figure('Position',[100 100 900 700])
            plotTriggered_new(r.trig_dat,'traces',tr,'fromLME',true, ...
                'plot_vel',false, 'trig_events',{'onsets','offsets'});
            save_img(gcf,figDir,sprintf('Sup4_%s_%s_%s_trig_%s',p{3},ds.name,e,tr{1})); close(gcf)
        end
        % D/H/L/P: triggered DI diff, with significance
        figure('Position',[100 100 900 700])
        plotTriggered_new(r.trig_dat,'traces',{'diff_di'},'fromLME',true, ...
            'plot_sig',true,'line_colors',{[0 0 0]}, ...
            'plot_vel',false, 'trig_events',{'onsets','offsets'});
        save_img(gcf,figDir,sprintf('Sup4_%s_%s_%s_trig_DIdiff',p{4},ds.name,e)); close(gcf)
    end
end
close all

%% Count sizes & report stats -> Excel tab 'Sup4'
fig = 'Sup4';
R = table();
ss = {'A-B, I-J', 'track traversal, position-binned', n_trial;
      'C-D, K-L', 'in-track locomotion bouts', n_bout;
      'E-H, M-P', 'spontaneous locomotion (GUI) bouts', n_gui};
for s = 1:size(ss,1)
    for ct = {'dSPN','iSPN'}
        R = [R; statsRow('Figure',fig,'Panel',ss{s,1},'Measure',['sample size (' ss{s,2} ')'], ...
                'Group',ct{1},'nStruct',ss{s,3}, ...
                'Note','counted on the convolved-data tables; I-P use the same mice, sessions, fields and neurons')]; %#ok<AGROW>
    end
end

% plotted values (A/C-type panels: no tests; B/D-type panels: Holm p for the significance bars)
blocks = {'Summary: sample sizes', R};
lbl = struct('track',{{'track position','position_bin'}}, 'gui',{{'bout progress','progress_bin'}});
for d = 1:numel(datasets)
    ds = datasets(d);
    for k = 1:2
        e = {'track','gui'}; e = e{k};
        p = ds.panels(k,:); r = res.(ds.name).(e); xn = lbl.(e){2};
        blocks = [blocks
            {sprintf('Sup4%s: %s %s dSPN, iSPN (LME mean, SEM) by %s', p{1}, e, ds.label, lbl.(e){1}), ...
                sourceDataCurve(r.binned, {'dSPN','iSPN'}, 'xName',xn)
             sprintf('Sup4%s: %s %s dSPN - iSPN (LME mean, SEM, Holm p) by %s', p{2}, e, ds.label, lbl.(e){1}), ...
                sourceDataCurve(r.binned, {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, 'pval',true, 'xName',xn)}]; %#ok<AGROW>
        for ev = {'onsets','offsets'}
            blocks = [blocks
                {sprintf('Sup4%s: %s %s %s-aligned dSPN, iSPN (LME mean, SEM)', p{3}, e, ds.label, ev{1}), ...
                    sourceDataCurve(r.trig_dat.(ev{1}), {'dSPN','iSPN'}, 'xName','time_s')
                 sprintf('Sup4%s: %s %s %s-aligned dSPN - iSPN (LME mean, SEM, Holm p)', p{4}, e, ds.label, ev{1}), ...
                    sourceDataCurve(r.trig_dat.(ev{1}), {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, 'pval',true, 'xName','time_s')}]; %#ok<AGROW>
        end
    end
end
blocks = [blocks
    {'Sample size per mouse: track trials', n_trial.perMouse
     'Sample size per mouse: track bouts', n_bout.perMouse
     'Sample size per mouse: GUI bouts', n_gui.perMouse}];
writeSourceData(statsFile, fig, blocks);

end
