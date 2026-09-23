function Sup9(cfg)
% Sup 9: Supplement for main Fig5
% B: distance/time field-centre distribution (dSPN vs iSPN)
% C: field-centre location per neuron and % tuned cells per session
% D: tuning-type pie charts (plotted, not reported)
% E: onset/offset triggered dSPN, iSPN and velocity; F: their difference
% G: acceleration by visual position and by bout distance
% (A examples: not generated here)
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure

% population data (E-G: all recorded neurons), tuning counts and single-cell
% data (B-D); the tuning-subtype table holds only tuned neurons, so it is
% not used for the population panels
popPath = cfg.infPop;
fpath   = cfg.infTuned;
scpath  = cfg.infSc;
dat = loadInfDataset('popPath', popPath, 'fpath', fpath, 'scpath', scpath);
pop_tbl = dat.pop.tbl;
pop_trial = dat.pop.trial_selection;
pop_bout  = dat.pop.trial_selection & dat.pop.bout_selection;   % binned panels (G)
pop_trig  = dat.pop.trial_selection & dat.pop.trig_selection;   % triggered panels (E, F)
binType = dat.pop.binType;
sc = dat.sc;

x_bout  = [1:100; linspace(1,150,100)];     % bout distance (cm)
x_track = [1:100; linspace(-120,40,100)];   % visual position (arb. unit)
plt_info_bout  = struct('startLoc',1,'rewardLoc',150);
plt_info_track = struct('startLoc',-120,'rewardLoc',40);
sr = 31;
on_win = [-31:3*31; linspace(-1,3,125)];
off_win = [3*-31:31; linspace(-3,1,125)];
ct_name = {'dSPN','iSPN'};

%% Plotting
% (from paperFigScripts/paperFigureSingleCell_Inf.m and
%  paperFigScripts/paperFigureSubpopulation_Inf.m)
figDir = fullfile(save2where, 'Sup9');
if ~isfolder(figDir), mkdir(figDir); end

% peak (field centre) of each tuned neuron, in bins and in cm
[~,dpeak] = nanmax(sc.drm_w1,[],2);
[~,ipeak] = nanmax(sc.irm_w1,[],2);
dist_peak = {dpeak(sc.dstb(:,1) & ~sc.dboth), ipeak(sc.istb(:,1) & ~sc.iboth)}; % distance tuned ONLY

% B: distribution of distance/time field centres, dSPN vs iSPN (distance tuned ONLY)
figure
h1 = histogram(dist_peak{1},10,'binLimits',[1 100],'Normalization','Probability', ...
    'DisplayStyle','Stairs','EdgeColor',[0.2 0.4 0.8],'LineWidth',2);
hold on;
h2 = histogram(dist_peak{2},10,'binLimits',[1 100],'Normalization','Probability',...
    'DisplayStyle','Stairs','EdgeColor',[0.8 0.6 0.2],'LineWidth',2);
T_fieldcentre = [renamevars(sourceDataHist(h1),'count','probability_dSPN'), ...
    renamevars(sourceDataHist(h2),{'bin_left','bin_right','count'}, ...
    {'bin_left_iSPN','bin_right_iSPN','probability_iSPN'})]; % read bars before closing
[~, pB] = kstest2(dist_peak{1}, dist_peak{2});
legend('dSPN','iSPN'); xlabel('Field center distance (bins)')
title(['Ks-test p = ' num2str(pB)])
save_img(gcf,figDir,'Sup9_B_sc_inf_histogram_PeakLoc_OnlyDistTune_DvsI'); close(gcf)

% C top: field width of visual-position tuned neurons, dSPN vs iSPN
% (sc_inf_OnlyDistTune_boxplotFieldWidthDvsI in paperFigureSingleCell_Inf.m)
dFieldWidth = computeFieldWidth(sc.dspn);
iFieldWidth = computeFieldWidth(sc.ispn);
C_top = {dFieldWidth(sc.dstb(:,1) & ~sc.dboth, 1), iFieldWidth(sc.istb(:,1) & ~sc.iboth,1)};
figure; boxViolin(C_top{1}, C_top{2}, get(gcf,'Number'));
pC_top = ranksum(C_top{1}, C_top{2});
title(['Only distance tuned field width: left dSPN, right iSPN; ranksum p = ' num2str(pC_top)])
ylabel('Field width')
save_img(gcf,figDir,'Sup9_C_top_sc_inf_OnlyDistTune_boxplotFieldWidth'); close(gcf)

% C bottom: percentage of distance/time tuned cells per session
possrecs = unique([unique(sc.drecnum), unique(sc.irecnum)]);
edges = min(possrecs)-0.5:1:max(possrecs)+0.5;
dnum = histcounts(sc.drecnum, edges);
inum = histcounts(sc.irecnum, edges);
dnum_dist = histcounts(sc.drecnum(sc.dstb(:,1)' & ~sc.dboth'), edges);
inum_dist = histcounts(sc.irecnum(sc.istb(:,1)' & ~sc.iboth'), edges);
C_bot = {100 * (dnum_dist./dnum)', 100 * (inum_dist./inum)'};
figure; boxScatterLine(C_bot{1}, C_bot{2}, get(gcf,'Number'));
ok = ~isnan(C_bot{1}) & ~isnan(C_bot{2});
pC_bot = signrank(C_bot{1}(ok), C_bot{2}(ok));
title(['% distance/time cells per session; signrank p = ' num2str(pC_bot)])
save_img(gcf,figDir,'Sup9_C_bottom_sc_inf_prctSess_DistTune_DvsI_signrank'); close(gcf)

% D: tuning-type pie charts (not reported)
counts = struct();
counts.dSPN = [sum(dat.info.nta_d), sum(dat.info.only_bDist_d), sum(dat.info.only_track_d), sum(sc.dboth)];
counts.iSPN = [sum(dat.info.nta_i), sum(dat.info.only_bDist_i), sum(dat.info.only_track_i), sum(sc.iboth)];
ntna = [size(sc.dstb,1) - sum(counts.dSPN), size(sc.istb,1) - sum(counts.iSPN)];
for c = 1:2
    figure
    piechart([counts.(ct_name{c}), ntna(c)], ...
        {'Not tuned/active', 'dist only', 'track only', 'both', 'Not tuned/not active'})
    title([ct_name{c} ': inf track'])
    save_img(gcf, figDir, sprintf('Sup9_D_Inf_%s_celltype_piechart', ct_name{c})); close(gcf)
end

% E, F: onset / offset triggered activity of all recorded neurons
trig_dat = struct( ...
    'onsets',  getPlotDatFromTbl(pop_tbl(pop_trig,:), ...
        'onsets','xdata',on_win,'y','avgFc3'), ...
    'offsets', getPlotDatFromTbl(pop_tbl(pop_trig,:), ...
        'offsets','xdata',off_win,'y','avgFc3'));
figure('Position',[100 100 900 700])
plotTriggered_new(trig_dat,'traces',{'dSPN','iSPN'},'fromLME',true, ...
    'plot_vel',false,'trig_events',{'onsets','offsets'});
save_img(gcf,figDir,'Sup9_E_inf_TriggeredAvg_D&I'); close(gcf)
figure('Position',[900 100 900 700])
plotTriggered_new(trig_dat,'traces',{'velocity'},'fromLME',true, ...
    'line_colors',{[0 0 1]},'plot_vel',false,'trig_events',{'onsets','offsets'});
save_img(gcf,figDir,'Sup9_E_inf_TriggeredAvg_velocity'); close(gcf)
figure('Position',[100 100 900 700])
plotTriggered_new(trig_dat,'traces',{'diff_di'},'fromLME',true, ...
    'plot_sig',true,'line_colors',{[0 0 0]},'plot_vel',false,'trig_events',{'onsets','offsets'});
save_img(gcf,figDir,'Sup9_F_inf_TriggeredAvg_DIdiff'); close(gcf)

% G: acceleration binned by visual position (top) and bout distance (bottom)
vuBin_tbl = groupsummary(pop_tbl(pop_trial,:), ...
    {'mouse','sess','field','repeat','vuBin','isD','isI','celltype'}, ...
    'mean',{'velocity','acceleration'});
pBin_tbl = groupsummary(pop_tbl(pop_bout,:), ...
    {'mouse','sess','field','bout',binType,'isD','isI','celltype'}, ...
    'mean',{'velocity','acceleration'});
velAcc = struct();
velAcc.track = getPlotDatFromTbl(vuBin_tbl,'vuBin','y','mean_acceleration', ...
    'vel_y','mean_velocity','xdata',x_track);
velAcc.bout = getPlotDatFromTbl(pBin_tbl,binType,'y','mean_acceleration', ...
    'vel_y','mean_velocity','xdata',x_bout);
figure('Position',[900 100 900 700])
plotBinned_new(velAcc.track,plt_info_track,'fromLME',true,'traces',{'dSPN'},'line_colors',{[0 0 1]})
xlabel('Inf Track unit (arb.)'); ylabel('Acceleration (cm/s^2)')
save_img(gcf,figDir,'Sup9_G_top_inf_TrackPosition_accel'); close(gcf)
figure('Position',[900 100 900 700])
plotBinned_new(velAcc.bout,plt_info_bout,'fromLME',true,'traces',{'dSPN'},'line_colors',{[0 0 1]})
xlabel('Bout Distance Bin (cm)'); ylabel('Acceleration (cm/s^2)')
save_img(gcf,figDir,'Sup9_G_bottom_inf_BoutDistBin_accel'); close(gcf)
close all

%% Count sizes & report stats -> Excel tab 'Sup9'
% D (pie charts) is not reported
fig = 'Sup9';
n_trial = countSampleSize(pop_tbl(pop_trial,:), 'repeat');
n_bout  = countSampleSize(pop_tbl(pop_trig,:), 'bout');
R = [statsRow('Figure',fig,'Panel','G top','Measure','sample size (visual position, trials)','nStruct',n_trial)
     statsRow('Figure',fig,'Panel','E-G','Measure','sample size (bouts > 150 cm)','nStruct',n_bout)
     statsRow('Figure',fig,'Panel','B','Measure','distance/time field centre distribution', ...
        'Group','dSPN vs iSPN','Test',"Kolmogorov-Smirnov",'pValue',pB, ...
        'N',numel(dist_peak{1}) + numel(dist_peak{2}),'N_unit',"neurons (both types)")];
R = [R; reportGroupCompare(ct_name, C_top, 'Figure',fig,'Panel','C top', ...
        'Measure','field width of visual-position tuned neurons, one point per neuron')
     reportGroupCompare(ct_name, C_bot, 'paired',true, 'Figure',fig,'Panel','C bottom', ...
        'Measure','% distance/time tuned cells per session')];

% plotted values (no tests in E and G; F: Holm p for the bars)
blocks = {'Summary: sample sizes; B, C tests', R
    'Sup9B: distance/time field centre histogram (probability per bin)', T_fieldcentre
    'Sup9C top: field width (visual-position tuned), one point per neuron', ...
        sourceDataPoints(ct_name, C_top, 'valueName','field_width')
    'Sup9C bottom: % distance/time tuned cells, one point per session (paired)', ...
        sourceDataPoints(ct_name, C_bot, 'paired',true, 'valueName','pct_distance_tuned')
    'Sup9G top: acceleration by visual position (LME mean, SEM)', ...
        sourceDataCurve(velAcc.track, {'dSPN'}, 'labels',{'acceleration'}, 'xName','track_position')
    'Sup9G bottom: acceleration by bout distance (LME mean, SEM)', ...
        sourceDataCurve(velAcc.bout, {'dSPN'}, 'labels',{'acceleration'}, 'xName','bout_distance_cm')};
for ev = {'onsets','offsets'}
    blocks = [blocks
        {sprintf('Sup9E: %s-aligned dSPN, iSPN and velocity (LME mean, SEM)', ev{1}), ...
            sourceDataCurve(trig_dat.(ev{1}), {'dSPN','iSPN','velocity'}, 'xName','time_s')
         sprintf('Sup9F: %s-aligned dSPN - iSPN (LME mean, SEM, Holm p)', ev{1}), ...
            sourceDataCurve(trig_dat.(ev{1}), {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, ...
                'pval',true, 'xName','time_s')}]; %#ok<AGROW>
end
blocks = [blocks
    {'Sample size per mouse: population repeats (visual position)', n_trial.perMouse
     'Sample size per mouse: population bouts (> 150 cm)', n_bout.perMouse}];
writeSourceData(statsFile, fig, blocks);

end
