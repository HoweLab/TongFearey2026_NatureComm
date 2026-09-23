function Figure4(cfg)
% Main Figure4 - novel vs. familiar track: imbalance and locomotion
% A-B: dSPN - iSPN binned by track position (familiar / novel)
% C-E: high-acceleration onsets; F-H: high-deceleration offsets
% I: velocity binned by track position
% (J-M behavioural metrics are not in this repository)
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure

% novel (first exposure) and familiar (expert) trials, one session per mouse
fpath = cfg.trackNovFam;
switchInfoPath = cfg.switchInfo;
dat = selectNovelFamiliar('fpath', fpath, 'switchInfoPath', switchInfoPath);
switch_tbl = dat.tbl;
swBinned_tbl = dat.binned;

xdata = [1:100; 1:100];
fake_info = struct('startLoc',1,'rewardLoc',100);
sr = 31;
conds = {'familiar', 0; 'novel', 1};          % name, info_trial_novel
cond_color = {[0 0 0], [0.5 0.2 0.7]};        % black / purple, as in the paper

%% Plotting
% (from NewAnalyses_SX/paperFigurePopulation_2w_Novel_Familar_SX.m)
figDir = fullfile(save2where, 'Fig4');
if ~isfolder(figDir), mkdir(figDir); end

% A, B: dSPN - iSPN binned by track position, per condition
binned_DI = struct();
for c = 1:2
    binned_DI.(conds{c,1}) = getPlotDatFromTbl( ...
        swBinned_tbl(swBinned_tbl.info_trial_novel == conds{c,2},:), ...
        'vuBin','compare','celltype','xdata',xdata,'y','mean_avgFc3','vel_y',[]);
    figure('Position',[900 100 900 700])
    plotBinned_new(binned_DI.(conds{c,1}),fake_info,'fromLME',true,'traces',{'diff_di'}, ...
        'plot_vel',false,'line_colors',cond_color(c),'plot_sig',true)
    title(sprintf('Fig4%s %s track', char('A'+c-1), conds{c,1}),'Interpreter','none')
    save_img(gcf,figDir,sprintf('Fig4_%s_2world_positionBins_DIdiff_%s', char('A'+c-1), conds{c,1}))
    close(gcf)
end

% C-H: onset / offset triggered activity for high-acceleration bouts
% (upper half of onsets/offsets ranked by peak acceleration / deceleration)
win_events = {'onsets','offsets'};
stats_vars = {'acc_peak','dec_peak'};
grp_vars = [{'mouse','field','sess','bout'};{'mouse','field','sess','bout'}];
exp_name = {'Fig4_fam','Fig4_nov'};
numRange = 2;                 % halves; level 2 = high acceleration / deceleration
bout_sel = switch_tbl.info_bout_intrack == 1;
accRange = select_VelAccRange_tbl( ...
    {switch_tbl(bout_sel & dat.fam,:), switch_tbl(bout_sel & dat.nov,:)}, ...
    'numRange', numRange, 'BinOrTrig', 'Trig', 'win_events', win_events, ...
    'wind', {-1*sr:2*sr, -2*sr:sr}, 'grp_vars', grp_vars);
trig_dat = plot_AccRange_DI(accRange, figDir, 'y','avgFc3', 'vel_y','acceleration', ...
    'bin_unit', win_events, 'wind_name', win_events, 'binName', 'Trig', ...
    'stats_vars', stats_vars, 'exp_name', exp_name);
close all

hi = numRange;                        % high-acceleration half
pan = {'C','D','E'; 'F','G','H'};     % rows: onset, offset
for w = 1:numel(win_events)
    x = accRange.bin_idx{w} / sr;
    fi = struct('startLoc',min(x),'rewardLoc',max(x));
    % C, F: acceleration for the two conditions
    figure('Position',[900 100 900 700]); hold on
    for c = 1:2
        v = trig_dat{c,w,hi}{1,1}.velocity;
        plot_error(gca, x, v.mu, v.sem, 'Color', cond_color{c});
    end
    title(sprintf('Fig4%s %s-aligned acceleration (black familiar, purple novel)', pan{w,1}, win_events{w}))
    save_img(gcf,figDir,sprintf('Fig4_%s_2world_%s_acceleration_famVsNov', pan{w,1}, win_events{w}))
    close(gcf)
    for c = 1:2
        d = trig_dat{c,w,hi}{1,1};
        d.xdata = x; % plot in seconds
        % D, G: dSPN - iSPN
        figure('Position',[900 100 900 700])
        plotBinned_new(d,fi,'fromLME',true,'traces',{'diff_di'}, ...
            'plot_vel',false,'line_colors',cond_color(c),'plot_sig',true)
        title(sprintf('Fig4%s %s %s-aligned dSPN - iSPN', pan{w,2}, conds{c,1}, win_events{w}),'Interpreter','none')
        save_img(gcf,figDir,sprintf('Fig4_%s_2world_%s_DIdiff_%s', pan{w,2}, win_events{w}, conds{c,1}))
        close(gcf)
        % E, H: dSPN and iSPN
        figure('Position',[900 100 900 700])
        plotBinned_new(d,fi,'fromLME',true,'traces',{'dSPN','iSPN'},'plot_vel',false)
        title(sprintf('Fig4%s %s %s-aligned dSPN and iSPN', pan{w,3}, conds{c,1}, win_events{w}),'Interpreter','none')
        save_img(gcf,figDir,sprintf('Fig4_%s_2world_%s_D&I_%s', pan{w,3}, win_events{w}, conds{c,1}))
        close(gcf)
    end
end

% I: velocity binned by track position, one gray line per session
binned_vel = getPlotDatFromTbl(swBinned_tbl(swBinned_tbl.celltype == 'dSPN',:),'vuBin', ...
    'compare','novelty','y','mean_velocity','xdata',xdata);
vel_sess = struct();
for c = 1:2
    figure('Position',[900 100 900 700])
    plotBinned_new(binned_vel,fake_info,'fromLME',true,'traces',conds(c,1), ...
        'plot_vel',false,'line_colors',cond_color(c), 'plot_title',[conds{c,1} ' track'])
    [~, vel_sess.(conds{c,1})] = addSessionLinePlot(gca, ...
        switch_tbl(switch_tbl.info_trial_novel == conds{c,2},:), 'vuBin', 'velocity');
    save_img(gcf,figDir,sprintf('Fig4_I_2world_positionBins_velocity_%s', conds{c,1}))
    close(gcf)
end
close all

%% Count sizes & report stats -> Excel tab 'Fig4'
fig = 'Fig4';
n_fam = countSampleSize(switch_tbl(dat.fam,:), 'trial');
n_nov = countSampleSize(switch_tbl(dat.nov,:), 'trial');
n_fam_bout = countSampleSize(switch_tbl(bout_sel & dat.fam,:), 'bout');
n_nov_bout = countSampleSize(switch_tbl(bout_sel & dat.nov,:), 'bout');

R = [statsRow('Figure',fig,'Panel','A, I','Measure','sample size (familiar track trials)','nStruct',n_fam)
     statsRow('Figure',fig,'Panel','B, I','Measure','sample size (novel track trials)','nStruct',n_nov)
     statsRow('Figure',fig,'Panel','C-H','Measure','sample size (familiar in-track bouts)','nStruct',n_fam_bout)
     statsRow('Figure',fig,'Panel','C-H','Measure','sample size (novel in-track bouts)','nStruct',n_nov_bout)];
% C-H: cut-off and sample size of the high-acceleration half actually plotted
for w = 1:numel(win_events)
    st = accRange.trial_stats{w};
    cut = st.rangeVal(:, strcmp(st.statsName, stats_vars{w}));
    R = [R; statsRow('Figure',fig,'Panel',strjoin(pan(w,:),','), ...
            'Measure',sprintf('%s median split cut-off (from familiar, applied to both)', stats_vars{w}), ...
            'Note',sprintf('low <= %.4g < high <= %.4g cm/s^2', cut))]; %#ok<AGROW>
    for c = 1:2
        tbl_w = accRange.bin_tbl{c,w};
        nL = countSampleSize(tbl_w(tbl_w.(stats_vars{w}) == hi,:), 'bout');
        R = [R; statsRow('Figure',fig,'Panel',strjoin(pan(w,:),','), ...
                'Measure',sprintf('sample size (%s, high %s half)', conds{c,1}, stats_vars{w}), ...
                'Group',conds{c,1},'nStruct',nL)]; %#ok<AGROW>
    end
end

% plotted values (no tests in C, E, H, I; A-B, D, G: Holm p for the significance bars)
pos_bins = unique(swBinned_tbl.vuBin(~isnan(swBinned_tbl.vuBin)));
blocks = {'Summary: sample sizes; C-H acceleration cut-offs', R
    'Fig4I: velocity (LME mean, SEM) by track position, familiar and novel', ...
        sourceDataCurve(binned_vel, {'familiar','novel'}, 'xName','position_bin')};
for c = 1:2
    blocks = [blocks
        {sprintf('Fig4%s: %s track, dSPN - iSPN by track position (LME mean, SEM, Holm p)', char('A'+c-1), conds{c,1}), ...
            sourceDataCurve(binned_DI.(conds{c,1}), {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, ...
                'pval',true, 'xName','position_bin')
         sprintf('Fig4I: %s track velocity, mean per session (gray lines)', conds{c,1}), ...
            sourceDataPerMouse(vel_sess.(conds{c,1}), pos_bins, 'xName','position_bin')}]; %#ok<AGROW>
end
for w = 1:numel(win_events)
    x = accRange.bin_idx{w} / sr;
    for c = 1:2
        d = trig_dat{c,w,hi}{1,1};
        blocks = [blocks
            {sprintf('Fig4%s: %s %s-aligned acceleration (LME mean, SEM)', pan{w,1}, conds{c,1}, win_events{w}), ...
                sourceDataCurve(d, {'velocity'}, 'labels',{'acceleration'}, 'x',x, 'xName','time_s')
             sprintf('Fig4%s: %s %s-aligned dSPN - iSPN (LME mean, SEM, Holm p)', pan{w,2}, conds{c,1}, win_events{w}), ...
                sourceDataCurve(d, {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, 'pval',true, 'x',x, 'xName','time_s')
             sprintf('Fig4%s: %s %s-aligned dSPN and iSPN (LME mean, SEM)', pan{w,3}, conds{c,1}, win_events{w}), ...
                sourceDataCurve(d, {'dSPN','iSPN'}, 'x',x, 'xName','time_s')}]; %#ok<AGROW>
    end
end
blocks = [blocks
    {'Sessions and trials picked per mouse (novel and familiar)', dat.perMouse
     'Sample size per mouse: familiar trials', n_fam.perMouse
     'Sample size per mouse: novel trials', n_nov.perMouse}];
writeSourceData(statsFile, fig, blocks);

end
