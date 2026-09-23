function Sup8(cfg)
% Sup 8: Supplement for main Fig4
% A-B: dSPN and iSPN binned by track position (familiar / novel)
% C-E: low-amplitude onsets; F-H: low-amplitude offsets
% I: velocity profiles of each mouse, familiar vs. novel
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
% (from NewAnalyses_SX/paperFigurePopulation_2w_Novel_Familar_SX.m; same
%  pipeline as Figure4, here for the low-amplitude half of the transitions)
figDir = fullfile(save2where, 'Sup8');
if ~isfolder(figDir), mkdir(figDir); end

% A, B: dSPN and iSPN binned by track position, per condition
binned_DI = struct();
for c = 1:2
    binned_DI.(conds{c,1}) = getPlotDatFromTbl( ...
        swBinned_tbl(swBinned_tbl.info_trial_novel == conds{c,2},:), ...
        'vuBin','compare','celltype','xdata',xdata,'y','mean_avgFc3','vel_y',[]);
    figure('Position',[900 100 900 700])
    plotBinned_new(binned_DI.(conds{c,1}),fake_info,'fromLME',true, ...
        'traces',{'dSPN','iSPN'},'plot_vel',false)
    title(sprintf('Sup8%s %s track', char('A'+c-1), conds{c,1}),'Interpreter','none')
    save_img(gcf,figDir,sprintf('Sup8_%s_2world_positionBins_D&I_%s', char('A'+c-1), conds{c,1}))
    close(gcf)
end

% C-H: onset / offset triggered activity for low-amplitude bouts
% (lower half of onsets/offsets ranked by peak acceleration / deceleration)
win_events = {'onsets','offsets'};
stats_vars = {'acc_peak','dec_peak'};
grp_vars = [{'mouse','field','sess','bout'};{'mouse','field','sess','bout'}];
exp_name = {'Sup8_fam','Sup8_nov'};
numRange = 2;                 % halves; level 1 = low acceleration / deceleration
bout_sel = switch_tbl.info_bout_intrack == 1;
accRange = select_VelAccRange_tbl( ...
    {switch_tbl(bout_sel & dat.fam,:), switch_tbl(bout_sel & dat.nov,:)}, ...
    'numRange', numRange, 'BinOrTrig', 'Trig', 'win_events', win_events, ...
    'wind', {-1*sr:2*sr, -2*sr:sr}, 'grp_vars', grp_vars);
trig_dat = plot_AccRange_DI(accRange, figDir, 'y','avgFc3', 'vel_y','acceleration', ...
    'bin_unit', win_events, 'wind_name', win_events, 'binName', 'Trig', ...
    'stats_vars', stats_vars, 'exp_name', exp_name);
close all

lo = 1;                               % low-amplitude half
pan = {'C','D','E'; 'F','G','H'};     % rows: onset, offset
for w = 1:numel(win_events)
    x = accRange.bin_idx{w} / sr;
    fi = struct('startLoc',min(x),'rewardLoc',max(x));
    % C, F: acceleration for the two conditions
    figure('Position',[900 100 900 700]); hold on
    for c = 1:2
        v = trig_dat{c,w,lo}{1,1}.velocity;
        plot_error(gca, x, v.mu, v.sem, 'Color', cond_color{c});
    end
    title(sprintf('Sup8%s %s-aligned acceleration (black familiar, purple novel)', pan{w,1}, win_events{w}))
    save_img(gcf,figDir,sprintf('Sup8_%s_2world_%s_acceleration_famVsNov', pan{w,1}, win_events{w}))
    close(gcf)
    for c = 1:2
        d = trig_dat{c,w,lo}{1,1};
        d.xdata = x; % plot in seconds
        % D, G: dSPN - iSPN
        figure('Position',[900 100 900 700])
        plotBinned_new(d,fi,'fromLME',true,'traces',{'diff_di'}, ...
            'plot_vel',false,'line_colors',cond_color(c),'plot_sig',true)
        title(sprintf('Sup8%s %s %s-aligned dSPN - iSPN', pan{w,2}, conds{c,1}, win_events{w}),'Interpreter','none')
        save_img(gcf,figDir,sprintf('Sup8_%s_2world_%s_DIdiff_%s', pan{w,2}, win_events{w}, conds{c,1}))
        close(gcf)
        % E, H: dSPN and iSPN
        figure('Position',[900 100 900 700])
        plotBinned_new(d,fi,'fromLME',true,'traces',{'dSPN','iSPN'},'plot_vel',false)
        title(sprintf('Sup8%s %s %s-aligned dSPN and iSPN', pan{w,3}, conds{c,1}, win_events{w}),'Interpreter','none')
        save_img(gcf,figDir,sprintf('Sup8_%s_2world_%s_D&I_%s', pan{w,3}, win_events{w}, conds{c,1}))
        close(gcf)
    end
end

% I: velocity profile of each mouse, familiar vs. novel
allmice = unique(swBinned_tbl.mouse);
ncol = ceil(sqrt(numel(allmice))); nrow = ceil(numel(allmice)/ncol);
vel_each = struct();
figure('Position',[50 50 1600 900])
for m = 1:numel(allmice)
    cur = swBinned_tbl(swBinned_tbl.mouse == allmice(m) & swBinned_tbl.celltype == 'dSPN',:);
    vel_each.(matlab.lang.makeValidName(string(allmice(m)))) = getPlotDatFromTbl(cur,'vuBin', ...
        'compare','novelty','y','mean_velocity','xdata',xdata, ...
        'randomEffects','(1|sess) + (1|sess:field)');
    v = vel_each.(matlab.lang.makeValidName(string(allmice(m))));
    subplot(nrow,ncol,m); hold on
    for c = 1:2
        plot_error(gca, v.xdata, v.(conds{c,1}).mu, v.(conds{c,1}).sem, 'Color', cond_color{c});
    end
    title(string(allmice(m)),'Interpreter','none'); xlabel('Track position'); ylabel('Velocity (cm/s)')
end
sgtitle('Sup8I velocity per mouse (black familiar, purple novel)')
save_img(gcf,figDir,'Sup8_I_2world_positionBins_velocity_eachMouse')
close all

%% Count sizes & report stats -> Excel tab 'Sup8'
fig = 'Sup8';
n_fam = countSampleSize(switch_tbl(dat.fam,:), 'trial');
n_nov = countSampleSize(switch_tbl(dat.nov,:), 'trial');
n_fam_bout = countSampleSize(switch_tbl(bout_sel & dat.fam,:), 'bout');
n_nov_bout = countSampleSize(switch_tbl(bout_sel & dat.nov,:), 'bout');

R = [statsRow('Figure',fig,'Panel','A, I','Measure','sample size (familiar track trials)','nStruct',n_fam)
     statsRow('Figure',fig,'Panel','B, I','Measure','sample size (novel track trials)','nStruct',n_nov)
     statsRow('Figure',fig,'Panel','C-H','Measure','sample size (familiar in-track bouts)','nStruct',n_fam_bout)
     statsRow('Figure',fig,'Panel','C-H','Measure','sample size (novel in-track bouts)','nStruct',n_nov_bout)];
% C-H: cut-off and sample size of the low-amplitude half actually plotted
for w = 1:numel(win_events)
    st = accRange.trial_stats{w};
    cut = st.rangeVal(:, strcmp(st.statsName, stats_vars{w}));
    R = [R; statsRow('Figure',fig,'Panel',strjoin(pan(w,:),','), ...
            'Measure',sprintf('%s median split cut-off (from familiar, applied to both)', stats_vars{w}), ...
            'Note',sprintf('low <= %.4g < high <= %.4g cm/s^2', cut))]; %#ok<AGROW>
    for c = 1:2
        tbl_w = accRange.bin_tbl{c,w};
        nL = countSampleSize(tbl_w(tbl_w.(stats_vars{w}) == lo,:), 'bout');
        R = [R; statsRow('Figure',fig,'Panel',strjoin(pan(w,:),','), ...
                'Measure',sprintf('sample size (%s, low %s half)', conds{c,1}, stats_vars{w}), ...
                'Group',conds{c,1},'nStruct',nL)]; %#ok<AGROW>
    end
end

% plotted values (no tests in A-C, E, F, H, I; D and G: Holm p for the bars)
blocks = {'Summary: sample sizes; C-H acceleration cut-offs', R};
for c = 1:2
    blocks = [blocks
        {sprintf('Sup8%s: %s track, dSPN and iSPN by track position (LME mean, SEM)', char('A'+c-1), conds{c,1}), ...
            sourceDataCurve(binned_DI.(conds{c,1}), {'dSPN','iSPN'}, 'xName','position_bin')}]; %#ok<AGROW>
end
for w = 1:numel(win_events)
    x = accRange.bin_idx{w} / sr;
    for c = 1:2
        d = trig_dat{c,w,lo}{1,1};
        blocks = [blocks
            {sprintf('Sup8%s: %s %s-aligned acceleration (LME mean, SEM)', pan{w,1}, conds{c,1}, win_events{w}), ...
                sourceDataCurve(d, {'velocity'}, 'labels',{'acceleration'}, 'x',x, 'xName','time_s')
             sprintf('Sup8%s: %s %s-aligned dSPN - iSPN (LME mean, SEM, Holm p)', pan{w,2}, conds{c,1}, win_events{w}), ...
                sourceDataCurve(d, {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, 'pval',true, 'x',x, 'xName','time_s')
             sprintf('Sup8%s: %s %s-aligned dSPN and iSPN (LME mean, SEM)', pan{w,3}, conds{c,1}, win_events{w}), ...
                sourceDataCurve(d, {'dSPN','iSPN'}, 'x',x, 'xName','time_s')}]; %#ok<AGROW>
    end
end
for m = 1:numel(allmice)
    mouse = matlab.lang.makeValidName(string(allmice(m)));
    blocks(end+1,:) = {sprintf('Sup8I: %s velocity, familiar and novel (LME mean, SEM)', mouse), ...
        sourceDataCurve(vel_each.(mouse), {'familiar','novel'}, 'xName','position_bin')};
end
blocks = [blocks
    {'Sessions and trials picked per mouse (novel and familiar)', dat.perMouse
     'Sample size per mouse: familiar trials', n_fam.perMouse
     'Sample size per mouse: novel trials', n_nov.perMouse}];
writeSourceData(statsFile, fig, blocks);

end
