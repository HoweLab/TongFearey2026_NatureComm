function Figure5(cfg)
% Main Figure5 - infinite track: imbalance organized by distance/time progress
% A bottom: distributions of reward and bout distances
% B: velocity by visual position / bout distance
% D-E: population rasters (number of neurons reported only)
% F-I: all neurons binned by bout distance (F, G) and visual position (H, I)
% J-M: distance/time tuned (J, K) and non-distance/time tuned (L, M)
% N: difference between the two populations' imbalance
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure
dspn_cmap = parula; ispn_cmap = parula;

% two datasets: population (all neurons) for A, B, F-I; tuning subtypes for J-N
popPath = cfg.infPop;
fpath   = cfg.infTuned;
scpath  = cfg.infSc;
dat = loadInfDataset('popPath', popPath, 'fpath', fpath, 'scpath', scpath);
pop_tbl = dat.pop.tbl;            % all recorded neurons (population panels)
tun_tbl = dat.tuned.tbl;          % avgFc3_<tuning type> columns (subpopulation panels)
binType = dat.pop.binType;        % bout-distance bin column of the population table

x_bout  = [1:100; linspace(1,150,100)];     % bout distance (cm)
x_track = [1:100; linspace(-120,40,100)];   % visual position (arb. unit)
plt_info_bout  = struct('startLoc',1,'rewardLoc',150);
plt_info_bout50 = struct('startLoc',50,'rewardLoc',150); % after the first 50 cm
plt_info_track = struct('startLoc',-120,'rewardLoc',40);

%% Plotting
% (from paperFigScripts/paperFigureSubpopulation_Inf.m and
%  paperFigScripts/paperFigureSingleCell_Inf.m)
figDir = fullfile(save2where, 'Fig5');
if ~isfolder(figDir), mkdir(figDir); end

% A bottom: bout distance distribution with the 150 cm cut-off
figure
h_boutdist = histogram(dat.pop.bout.totldist, 50);
T_boutdist = sourceDataHist(h_boutdist);  % read bars before the figure is closed
hold on; xline(dat.min_boutdist,'LineWidth',2,'Color','r')
xlabel('bout distance (cm)'); ylabel('running bout counts')
save_img(gcf,figDir,'Fig5_A_inf_bout_distance_distribution'); close(gcf)

% A bottom: reward distance per trial (from paperFigurePopulation_Inf_BoutDist.m)
T_rewdist = table();
if ismember('rewdist', pop_tbl.Properties.VariableNames)
    allmice = unique(pop_tbl.mouse);
    allrew = cell(1, numel(allmice));
    for i = 1:numel(allmice)
        curtbl = pop_tbl(pop_tbl.mouse == allmice(i),:);
        [sessIdx,sess] = findgroups(curtbl(:,{'sess'}));
        rewdist = [];
        tracklgth = unique(curtbl.trackLgth);
        for s = 1:height(sess)
            t = curtbl(sessIdx == s,:);
            [grpsN,tidN] = findgroups(t(:,{'sess','trial'}));
            [~, min_idx] = min(str2double(string(tidN.trial)));
            alltrialid = unique(grpsN(~isnan(grpsN)));
            alltrialid = alltrialid(~ismember(alltrialid, min_idx));
            for j = 1:length(alltrialid)
                % virmen unit -> meter
                rewdist = [rewdist, max(t.rewdist(grpsN==alltrialid(j))) / (160/tracklgth)]; %#ok<AGROW>
            end
        end
        allrew{i} = rewdist;
    end
    figure
    h_rewdist = histogram(cell2mat(allrew) * 100, 50);
    T_rewdist = sourceDataHist(h_rewdist);
    xlabel('Real distance traveled to reward for each trial in cm'); ylabel('trial count')
    save_img(gcf,figDir,'Fig5_A_inf_reward_distance_distribution'); close(gcf)
else
    warning('No rewdist column in this table: reward-distance histogram skipped')
end

% binned tables
% population: by visual position (repeats) and by bout distance (bouts)
pop_trial = dat.pop.trial_selection;
pop_bout  = dat.pop.trial_selection & dat.pop.bout_selection;
vuBin_tbl = groupsummary(pop_tbl(pop_trial,:), ...
    {'mouse','sess','field','repeat','vuBin','isD','isI','celltype'}, ...
    'mean',{'avgFc3','velocity','acceleration'});
pBin_tbl = groupsummary(pop_tbl(pop_bout,:), ...
    {'mouse','sess','field','bout',binType,'isD','isI','celltype'}, ...
    'mean',{'avgFc3','velocity','acceleration','lick_count'});
tun_bout = dat.tuned.bout_selection;
tunBin_tbl = groupsummary(tun_tbl(tun_bout,:), ...
    {'mouse','sess','field','bout','boutDistBin','isD','isI','celltype'}, ...
    'mean',[strcat('avgFc3_',dat.ext), {'velocity','acceleration','lick_count'}]);

% B: velocity binned by visual position (top) and bout distance (bottom)
velAcc = struct();
velAcc.track = getPlotDatFromTbl(vuBin_tbl,'vuBin','y','mean_acceleration', ...
    'vel_y','mean_velocity','xdata',x_track);
velAcc.bout = getPlotDatFromTbl(pBin_tbl,binType,'y','mean_acceleration', ...
    'vel_y','mean_velocity','xdata',x_bout);
figure('Position',[900 100 900 700])
plotBinned_new(velAcc.track,plt_info_track,'fromLME',true,'traces',{'velocity'},'line_colors',{[0 0 1]})
xlabel('Inf Track unit (arb.)')
save_img(gcf,figDir,'Fig5_B_top_inf_TrackPosition_velocity'); close(gcf)
figure('Position',[900 100 900 700])
plotBinned_new(velAcc.bout,plt_info_bout,'fromLME',true,'traces',{'velocity'},'line_colors',{[0 0 1]})
xlabel('Bout Distance Bin (cm)')
save_img(gcf,figDir,'Fig5_B_bottom_inf_BoutDistBin_velocity'); close(gcf)

% D, E: rasters of distance/time tuned (D) and visual-position tuned (E) neurons
sc = dat.sc;
iPlot = false;
dw1 = sortMatrixCorMat(smoothdata(sc.drm_w1(sc.dstb(:,1),:),2), ...
    smoothdata(sc.drm_w2(sc.dstb(:,1),:),2),iPlot); % distance tuned
dw2 = sortMatrixCorMat(smoothdata(sc.drm_w2(sc.dstb(:,2),:),2), ...
    smoothdata(sc.drm_w1(sc.dstb(:,2),:),2),iPlot); % position tuned
iw1 = sortMatrixCorMat(smoothdata(sc.irm_w1(sc.istb(:,1),:),2), ...
    smoothdata(sc.irm_w2(sc.istb(:,1),:),2),iPlot);
iw2 = sortMatrixCorMat(smoothdata(sc.irm_w2(sc.istb(:,2),:),2), ...
    smoothdata(sc.irm_w1(sc.istb(:,2),:),2),iPlot);
figure(1);
subplot(221);colormap(dspn_cmap); imagesc(dw1.maxn(dw1.sortind,:));
title(['number of cell : ' num2str(numel(dw1.sortind))])
subplot(222);colormap(dspn_cmap); imagesc(dw1.maxn2(dw1.sortind,:));
subplot(223);colormap(dspn_cmap); imagesc(dw2.maxn(dw2.sortind,:));
title(['number of cell : ' num2str(numel(dw2.sortind))])
subplot(224);colormap(dspn_cmap); imagesc(dw2.maxn2(dw2.sortind,:));
sgtitle('dspn: top distance/time tuned (D), bottom visual position tuned (E)')
figure(2);
subplot(221);colormap(ispn_cmap); imagesc(iw1.maxn(iw1.sortind,:));
title(['number of cell : ' num2str(numel(iw1.sortind))])
subplot(222);colormap(ispn_cmap); imagesc(iw1.maxn2(iw1.sortind,:));
subplot(223);colormap(ispn_cmap); imagesc(iw2.maxn(iw2.sortind,:));
title(['number of cell : ' num2str(numel(iw2.sortind))])
subplot(224);colormap(ispn_cmap); imagesc(iw2.maxn2(iw2.sortind,:));
sgtitle('ispn: top distance/time tuned (D), bottom visual position tuned (E)')
save_img(figure(1),figDir,'Fig5_DE_sc_inf_raster_dSPN_onlyDistOrPos')
save_img(figure(2),figDir,'Fig5_DE_sc_inf_raster_iSPN_onlyDistOrPos')
close all

% F-M: population and subpopulation activity
% F, G: all neurons by bout distance; H, I: all neurons by visual position
% J, K: distance/time tuned only; L, M: non-distance/time tuned
plot_sets = struct( ...
    'name',  {'all_bout','all_track','distTuned','nonTuned'}, ...
    'y',     {'mean_avgFc3','mean_avgFc3','mean_avgFc3_only_bDist','mean_avgFc3_nta'}, ...
    'src',   {'pop_bout','pop_track','tuned','tuned'}, ...
    'panel', {{'F','G'},{'H','I'},{'J','K'},{'L','M'}}, ...
    'label', {'all neurons, bout distance','all neurons, visual position', ...
              'distance/time tuned','non-distance/time tuned'});
binned = struct();
for s = 1:numel(plot_sets)
    ss = plot_sets(s);
    switch ss.src
        case 'pop_track' % all neurons, visual position
            tbl = vuBin_tbl; bin = 'vuBin'; xd = x_track;
            fi = plt_info_track; xlab = 'Inf Track unit (arb.)';
        case 'pop_bout'  % all neurons, bout distance
            tbl = pBin_tbl; bin = binType; xd = x_bout;
            fi = plt_info_bout; xlab = 'Bout Distance Bin (cm)';
        case 'tuned'     % tuning subtypes, bout distance, after the first 50 cm
            tbl = tunBin_tbl; bin = 'boutDistBin'; xd = x_bout;
            fi = plt_info_bout50; xlab = 'Bout Distance Bin (cm)';
    end
    binned.(ss.name) = getPlotDatFromTbl(tbl,bin,'y',ss.y,'vel_y','mean_velocity','xdata',xd);
    figure('Position',[900 100 900 700])
    plotBinned_new(binned.(ss.name),fi,'fromLME',true,'traces',{'dSPN','iSPN'},'plot_vel',false)
    title(sprintf('Fig5%s %s', ss.panel{1}, ss.label),'Interpreter','none'); xlabel(xlab)
    save_img(gcf,figDir,sprintf('Fig5_%s_inf_D&I_%s', ss.panel{1}, ss.name)); close(gcf)
    figure('Position',[900 100 900 700])
    plotBinned_new(binned.(ss.name),fi,'fromLME',true,'traces',{'diff_di'}, ...
        'plot_sig',true,'line_colors',{[0 0 0]},'plot_vel',false)
    title(sprintf('Fig5%s %s', ss.panel{2}, ss.label),'Interpreter','none'); xlabel(xlab)
    save_img(gcf,figDir,sprintf('Fig5_%s_inf_DIdiff_%s', ss.panel{2}, ss.name)); close(gcf)
end

% N: distance/time tuned vs. non-distance/time tuned imbalance
pBin_dist = tunBin_tbl; pBin_dist.mean_avgFc3 = tunBin_tbl.mean_avgFc3_only_bDist;
pBin_nta  = tunBin_tbl; pBin_nta.mean_avgFc3  = tunBin_tbl.mean_avgFc3_nta;
[N_dat, figs] = get_DIdiff_glm_SX(pBin_dist, pBin_nta, 'field','boutDistBin', ...
    'xdata', x_bout, 'fake_info', plt_info_bout50, 'vel_y', []);
save_img(figs.figExp01,figDir,'Fig5_N_inf_BoutDistBin_sigDiff_distTune_vs_nonTuned')
close all

%% Count sizes & report stats -> Excel tab 'Fig5'
% D-E rasters: only the number of neurons
fig = 'Fig5';
n_trial = countSampleSize(pop_tbl(pop_trial,:), 'repeat');
n_bout  = countSampleSize(pop_tbl(pop_bout,:), 'bout');
n_tuned = countSampleSize(tun_tbl(tun_bout,:), 'bout');
R = [statsRow('Figure',fig,'Panel','B top, H, I','Measure','sample size (population, visual position)','nStruct',n_trial)
     statsRow('Figure',fig,'Panel','B bottom, F, G','Measure','sample size (population, bouts > 150 cm)','nStruct',n_bout)
     statsRow('Figure',fig,'Panel','J-N','Measure','sample size (tuning subtypes, bouts > 150 cm)','nStruct',n_tuned)
     statsRow('Figure',fig,'Panel','A','Measure','bouts in distance histogram', ...
        'N',height(dat.pop.bout),'N_unit',"bouts",'Note','all locomotion bouts, before the 150 cm cut-off')];
% D, E: number of neurons in each raster
rasters = {'D','dSPN distance/time tuned', numel(dw1.sortind)
           'D','iSPN distance/time tuned', numel(iw1.sortind)
           'E','dSPN visual position tuned', numel(dw2.sortind)
           'E','iSPN visual position tuned', numel(iw2.sortind)};
for r = 1:size(rasters,1)
    R = [R; statsRow('Figure',fig,'Panel',rasters{r,1},'Measure','number of neurons in raster', ...
            'Group',rasters{r,2},'N',rasters{r,3},'N_unit',"neurons",'N_neurons',rasters{r,3})]; %#ok<AGROW>
end
% F-N: number of neurons per population (per-session counts of the sessions used)
% (all-neuron counts for F-I come from the sample-size rows above, via numFc3)
cellcounts = {'J,K,N','distance/time tuned only','only_bDist_d','only_bDist_i'
              'L,M,N','non-distance/time tuned','nta_d','nta_i'
              'D,E','distance/time tuned (any)','bDist_d','bDist_i'
              'D,E','visual position tuned (any)','track_d','track_i'};
for r = 1:size(cellcounts,1)
    for c = 1:2
        ct = {'dSPN','iSPN'};
        nCell = sum(dat.info.(cellcounts{r,2+c}));
        R = [R; statsRow('Figure',fig,'Panel',cellcounts{r,1}, ...
                'Measure',sprintf('number of neurons (%s)', cellcounts{r,2}), 'Group',ct{c}, ...
                'N',nCell,'N_unit',"neurons",'N_neurons',nCell)]; %#ok<AGROW>
    end
end

% plotted values (no tests in B, F, H, J, L; G, I, K, M, N: Holm p for the bars)
blocks = {'Summary: sample sizes; D-E raster and population neuron counts', R
    'Fig5A: bout distance histogram (cut-off 150 cm)', T_boutdist
    'Fig5B top: velocity by visual position (LME mean, SEM)', ...
        sourceDataCurve(velAcc.track, {'velocity'}, 'xName','track_position')
    'Fig5B bottom: velocity by bout distance (LME mean, SEM)', ...
        sourceDataCurve(velAcc.bout, {'velocity'}, 'xName','bout_distance_cm')};
if ~isempty(T_rewdist)
    blocks(end+1,:) = {'Fig5A: reward distance histogram (cm)', T_rewdist};
end
for s = 1:numel(plot_sets)
    ss = plot_sets(s);
    xn = 'bout_distance_cm'; if strcmp(ss.src,'pop_track'), xn = 'track_position'; end
    blocks = [blocks
        {sprintf('Fig5%s: %s, dSPN and iSPN (LME mean, SEM)', ss.panel{1}, ss.label), ...
            sourceDataCurve(binned.(ss.name), {'dSPN','iSPN'}, 'xName',xn)
         sprintf('Fig5%s: %s, dSPN - iSPN (LME mean, SEM, Holm p)', ss.panel{2}, ss.label), ...
            sourceDataCurve(binned.(ss.name), {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, ...
                'pval',true, 'xName',xn)}]; %#ok<AGROW>
end
blocks = [blocks
    {'Fig5N: (dSPN - iSPN) distance/time tuned - non-tuned (LME mean, SEM, Holm p)', ...
        sourceDataCurve(N_dat, {'diff_DI_exp01'}, 'labels',{'tuned_minus_nonTuned_DIdiff'}, ...
            'pval',true, 'xName','bout_distance_cm')
     'Sample size per mouse: population repeats (visual position)', n_trial.perMouse
     'Sample size per mouse: population bouts (> 150 cm)', n_bout.perMouse
     'Sample size per mouse: tuning-subtype bouts (> 150 cm)', n_tuned.perMouse}];
writeSourceData(statsFile, fig, blocks);

end
