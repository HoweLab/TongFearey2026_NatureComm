function Figure3(cfg)
% Main Figure3 - 2w single cell raster and 2w track imbalance
% B: velocity per track; D-E: population rasters; F-H: dSPN/iSPN per
% subpopulation; I-K: their dSPN - iSPN difference
% (A and C are schematics / single-cell examples: not generated here)
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure
dspn_cmap = parula; ispn_cmap = parula; % color map

% Load 2w single cell data
fpath  = cfg.tuned2w;
scpath = cfg.sc2w;
Tpath  = cfg.switchInfo;
dat = load2wDataset('fpath', fpath, 'scpath', scpath, 'Tpath', Tpath);
track_2w_tbl = dat.tbl;
trial_selection = dat.trial_selection;

xdata = [1:100; linspace(1,100,100)];
fake_info = struct('startLoc',1,'rewardLoc',100);

%% Plotting
figDir = fullfile(save2where, 'Fig3');
if ~isfolder(figDir), mkdir(figDir); end

% (from paperFigScripts/paperFigureSubPopulation.m, "COMPARE TWO WORLDS")
vuBin_tbl = groupsummary(track_2w_tbl(trial_selection,:), ...
    {'mouse','sess','field','trial','vuBin','isD','isI','celltype','info_trial_novel'}, ...
    'mean',{'avgFc3_rmp_overlap','avgFc3_tnrmp','avgFc3_nta','velocity','acceleration'});
fam_idx = vuBin_tbl.info_trial_novel == 0; % w1 - familiar world = Track 1
nov_idx = vuBin_tbl.info_trial_novel == 1; % w2 - novel world   = Track 2

binned_vel = getPlotDatFromTbl(vuBin_tbl(vuBin_tbl.celltype == 'dSPN',:),'vuBin', ...
    'compare','novelty','y','mean_velocity','xdata',xdata);

% B: velocity from two tracks
figure('Position',[900 100 900 700])
plotBinned_new(binned_vel,fake_info,'fromLME',true,'traces',{'familiar'}, ...
    'plot_vel',false,'line_colors',{[0 0 1],[0.2 0.8 0.8]}, ...
    'plot_title', 'Track 1 (familiar world)')
numMice = numel(unique(track_2w_tbl(trial_selection,:).mouse));
[~, vel_mouse_w1] = addSingleMouseLinePlot(gca, track_2w_tbl(trial_selection & track_2w_tbl.info_trial_novel == 0,:), ...
    'vuBin', 'velocity','fam', 0, numMice);
save_img(gcf,figDir,'Fig3_B_2w_positionBins_vel_track1_Fam'); close(gcf)

figure('Position',[900 100 900 700])
plotBinned_new(binned_vel,fake_info,'fromLME',true,'traces',{'novel'}, ...
    'plot_vel',false,'line_colors',{[0 0 1],[0.2 0.8 0.8]}, ...
    'plot_title', 'Track 2 (novel world)')
[~, vel_mouse_w2] = addSingleMouseLinePlot(gca, track_2w_tbl(trial_selection & track_2w_tbl.info_trial_novel == 1,:), ...
    'vuBin', 'velocity','novel', 0, numMice);
save_img(gcf,figDir,'Fig3_B_2w_positionBins_vel_track2_Novel'); close(gcf)

% D, E: Plot raster (single-cell data already restricted to valid neurons)
drm_w1 = dat.sc.drm_w1;  irm_w1 = dat.sc.irm_w1;
drm_w2 = dat.sc.drm_w2;  irm_w2 = dat.sc.irm_w2;
dstb   = dat.sc.dstb;    istb   = dat.sc.istb;
dremap = dat.sc.dremap;  iremap = dat.sc.iremap;
iPlot  = 0;

% D: track-sensitive (remapping) neurons
dw1 = sortMatrixCorMat(smoothdata(drm_w1(dstb(:,1) & dremap,:),2), ...
    smoothdata(drm_w2(dstb(:,1) & dremap,:),2), 0); % neurons tuned to world 1 and track sensitive
dw2 = sortMatrixCorMat(smoothdata(drm_w2(dstb(:,2) & dremap,:),2), ...
    smoothdata(drm_w1(dstb(:,2) & dremap,:),2), 0); % neurons tuned to world 2 and track sensitive
figure(1);
subplot(221);colormap(dspn_cmap); imagesc(dw1.maxn(dw1.sortind,:));
title(['number of cell : ' num2str(numel(dw1.sortind))])
subplot(222);colormap(dspn_cmap); imagesc(dw1.maxn2(dw1.sortind,:));
subplot(223);colormap(dspn_cmap); imagesc(dw2.maxn2(dw2.sortind,:));
title(['number of cell : ' num2str(numel(dw2.sortind))])
subplot(224);colormap(dspn_cmap); imagesc(dw2.maxn(dw2.sortind,:));

% ispn  remap
iw1 = sortMatrixCorMat(smoothdata(irm_w1(istb(:,1) & iremap,:),2), ...
    smoothdata(irm_w2(istb(:,1) & iremap,:),2),iPlot);
iw2 = sortMatrixCorMat(smoothdata(irm_w2(istb(:,2) & iremap,:),2), ...
    smoothdata(irm_w1(istb(:,2) & iremap,:),2),iPlot);
figure(2);
subplot(221);colormap(ispn_cmap); imagesc(iw1.maxn(iw1.sortind,:));
title(['number of cell : ' num2str(numel(iw1.sortind))])
subplot(222);colormap(ispn_cmap); imagesc(iw1.maxn2(iw1.sortind,:));
subplot(223);colormap(ispn_cmap); imagesc(iw2.maxn2(iw2.sortind,:));
title(['number of cell : ' num2str(numel(iw2.sortind))])
subplot(224);colormap(ispn_cmap); imagesc(iw2.maxn(iw2.sortind,:));

% E: track-insensitive (no remapping) neurons
% dspn no remap
dstbonly = any(dstb,2) & ~dremap; % stable cells (either world) and NO remap
dw1_nr = sortMatrixCorMat(smoothdata(drm_w1(dstbonly,:),2), ...
    smoothdata(drm_w2(dstbonly,:),2),iPlot);
dw2_nr = sortMatrixCorMat(smoothdata(drm_w2(dstbonly,:),2), ...
    smoothdata(drm_w1(dstbonly,:),2),iPlot);
figure(3)
subplot(221);colormap(dspn_cmap); imagesc(dw1_nr.maxn(dw1_nr.sortind,:));
title(['number of cell : ' num2str(numel(dw1_nr.sortind))])
subplot(222);colormap(dspn_cmap); imagesc(dw1_nr.maxn2(dw1_nr.sortind,:));
subplot(223);colormap(dspn_cmap); imagesc(dw2_nr.maxn2(dw1_nr.sortind,:));
title(['number of cell : ' num2str(numel(dw2_nr.sortind))])
subplot(224);colormap(dspn_cmap); imagesc(dw2_nr.maxn(dw1_nr.sortind,:));

% ispn no remap
istbonly = any(istb,2) & ~iremap;
iw1_nr = sortMatrixCorMat(smoothdata(irm_w1(istbonly,:),2), ...
    smoothdata(irm_w2(istbonly,:),2),iPlot);
iw2_nr = sortMatrixCorMat(smoothdata(irm_w2(istbonly,:),2), ...
    smoothdata(irm_w1(istbonly,:),2),iPlot);
figure(4)
subplot(221);colormap(ispn_cmap); imagesc(iw1_nr.maxn(iw1_nr.sortind,:));
title(['number of cell : ' num2str(numel(iw1_nr.sortind))])
subplot(222);colormap(ispn_cmap); imagesc(iw1_nr.maxn2(iw1_nr.sortind,:));
subplot(223);colormap(ispn_cmap); imagesc(iw2_nr.maxn2(iw2_nr.sortind,:));
title(['number of cell : ' num2str(numel(iw2_nr.sortind))])
subplot(224);colormap(ispn_cmap);imagesc(iw2_nr.maxn(iw2_nr.sortind,:));

save_img(figure(1),figDir,'Fig3_D_sc_2w_raster_dSPN_remap')
save_img(figure(2),figDir,'Fig3_D_sc_2w_raster_iSPN_remap')
save_img(figure(3),figDir,'Fig3_E_sc_2w_raster_dSPN_NOremap')
save_img(figure(4),figDir,'Fig3_E_sc_2w_raster_iSPN_NOremap')
close all

% F-K: subpopulation activity on Track 1 (familiar) trials
% F, I: position tuned and track sensitive (world-1 tuned)
% G, J: position tuned but track insensitive
% H, K: not position tuned but active
pops = struct( ...
    'name',  {'rmp_w1', 'tnrmp', 'nta'}, ...
    'y',     {'mean_avgFc3_rmp_overlap', 'mean_avgFc3_tnrmp', 'mean_avgFc3_nta'}, ...
    'panel', {{'F','I'}, {'G','J'}, {'H','K'}}, ...
    'label', {'tuned, track-sensitive', 'tuned, track-insensitive', 'non-position tuned (active)'}, ...
    'count', {{'rmp_w1_d','rmp_w1_i'}, {'tnrmp_d','tnrmp_i'}, {'nta_d','nta_i'}});
binned_pop = struct();
for p = 1:numel(pops)
    pp = pops(p);
    binned_pop.(pp.name) = getPlotDatFromTbl(vuBin_tbl(fam_idx,:),'vuBin','y',pp.y, ...
        'vel_y',[],'xdata',xdata);
    % F, G, H: dSPN and iSPN
    figure('Position',[900 100 900 700]);
    plotBinned_new(binned_pop.(pp.name),fake_info,'fromLME',true,'traces',{'dSPN'},'plot_vel',false)
    plotBinned_new(binned_pop.(pp.name),fake_info,'fromLME',true,'traces',{'iSPN'}, ...
        'line_colors',{[0 1 0]},'plot_vel',false);
    title(sprintf('Fig3%s %s', pp.panel{1}, pp.label),'Interpreter','none')
    save_img(gcf,figDir,sprintf('Fig3_%s_2w_positionBin_D&I_%s', pp.panel{1}, pp.name)); close(gcf)
    % I, J, K: dSPN - iSPN, with significance
    figure('Position',[900 100 900 700]);
    plotBinned_new(binned_pop.(pp.name),fake_info,'fromLME',true,'traces',{'diff_di'}, ...
        'plot_sig',true,'line_colors',{[0 0 0]},'plot_vel',false);
    title(sprintf('Fig3%s %s', pp.panel{2}, pp.label),'Interpreter','none')
    save_img(gcf,figDir,sprintf('Fig3_%s_2w_positionBin_Diff_%s', pp.panel{2}, pp.name)); close(gcf)
end
close all

%% Count sizes & report stats -> Excel tab 'Fig3'
% raster panels (D, E): only the number of cells; line plots (B, F-K): plotted values
fig = 'Fig3';
n_w1 = countSampleSize(track_2w_tbl(trial_selection & track_2w_tbl.info_trial_novel == 0,:), 'trial');
n_w2 = countSampleSize(track_2w_tbl(trial_selection & track_2w_tbl.info_trial_novel == 1,:), 'trial');
n_all = countSampleSize(track_2w_tbl(trial_selection,:), 'trial');

R = [statsRow('Figure',fig,'Panel','B, F-K','Measure','sample size (2 world, all trials)','nStruct',n_all)
     statsRow('Figure',fig,'Panel','B, F-K','Measure','sample size (Track 1 / familiar trials)','nStruct',n_w1)
     statsRow('Figure',fig,'Panel','B','Measure','sample size (Track 2 / novel trials)','nStruct',n_w2)];

% D, E: number of cells in each raster
rasters = {'D','dSPN tuned in Track 1, track-sensitive', numel(dw1.sortind)
           'D','dSPN tuned in Track 2, track-sensitive', numel(dw2.sortind)
           'D','iSPN tuned in Track 1, track-sensitive', numel(iw1.sortind)
           'D','iSPN tuned in Track 2, track-sensitive', numel(iw2.sortind)
           'E','dSPN position tuned, track-insensitive', numel(dw1_nr.sortind)
           'E','iSPN position tuned, track-insensitive', numel(iw1_nr.sortind)};
for r = 1:size(rasters,1)
    R = [R; statsRow('Figure',fig,'Panel',rasters{r,1},'Measure','number of neurons in raster', ...
            'Group',rasters{r,2},'N',rasters{r,3},'N_unit',"neurons",'N_neurons',rasters{r,3})]; %#ok<AGROW>
end

% F-K: number of neurons per subpopulation (per-session counts of the
% sessions used, from load2wDataset)
for p = 1:numel(pops)
    pp = pops(p);
    for c = 1:2
        ct = {'dSPN','iSPN'};
        nCell = sum([dat.info.(pp.count{c})]);
        R = [R; statsRow('Figure',fig,'Panel',[pp.panel{1} ',' pp.panel{2}], ...
                'Measure',sprintf('number of neurons (%s)', pp.label), 'Group',ct{c}, ...
                'N',nCell, 'N_unit',"neurons", 'N_neurons',nCell, 'nStruct',n_w1)]; %#ok<AGROW>
    end
end

% plotted values (B, F-H: no tests; I-K: Holm p for the significance bars)
% bins used by addSingleMouseLinePlot, per track
binsOf = @(nov) unique(track_2w_tbl.vuBin(trial_selection & track_2w_tbl.info_trial_novel == nov ...
    & ~isnan(track_2w_tbl.vuBin)));
blocks = {
    'Summary: sample sizes; D-E raster cell counts; F-K subpopulation cell counts', R
    'Fig3B: velocity (LME mean, SEM) by track position, Track 1 and Track 2', ...
        sourceDataCurve(binned_vel, {'familiar','novel'}, ...
            'labels',{'track1_familiar','track2_novel'}, 'xName','position_bin')
    'Fig3B: Track 1 velocity, mean per mouse (colored lines)', ...
        sourceDataPerMouse(vel_mouse_w1, binsOf(0), 'xName','position_bin')
    'Fig3B: Track 2 velocity, mean per mouse (colored lines)', ...
        sourceDataPerMouse(vel_mouse_w2, binsOf(1), 'xName','position_bin')};
for p = 1:numel(pops)
    pp = pops(p); b = binned_pop.(pp.name);
    blocks = [blocks
        {sprintf('Fig3%s: %s, dSPN and iSPN (LME mean, SEM)', pp.panel{1}, pp.label), ...
            sourceDataCurve(b, {'dSPN','iSPN'}, 'xName','position_bin')
         sprintf('Fig3%s: %s, dSPN - iSPN (LME mean, SEM, Holm p)', pp.panel{2}, pp.label), ...
            sourceDataCurve(b, {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, 'pval',true, 'xName','position_bin')}]; %#ok<AGROW>
end
blocks = [blocks
    {'Sample size per mouse: Track 1 (familiar) trials', n_w1.perMouse
     'Sample size per mouse: Track 2 (novel) trials', n_w2.perMouse}];
writeSourceData(statsFile, fig, blocks);

end
