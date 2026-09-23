function Sup6(cfg)
% Sup 6: Supplement for main Fig3
% Characterization of track-sensitive and track-insensitive position tuning
% C: cross-track correlation matrices (number of neurons reported only)
% D: correlation vs. position; E: per-neuron cross-track correlation
% F: spatial error per session; G: % track-sensitive per session
% H: field centers; I: mean dF/F per neuron; J: field width
% (A-B are single-cell examples: not generated here)
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure
dspn_cmap = parula; ispn_cmap = parula; % color map

% Load 2w single cell data (same sessions/mice excluded as Fig3)
fpath  = cfg.tuned2w;
scpath = cfg.sc2w;
Tpath  = cfg.switchInfo;
dat = load2wDataset('fpath', fpath, 'scpath', scpath, 'Tpath', Tpath);

drm_w1 = dat.sc.drm_w1;  irm_w1 = dat.sc.irm_w1;
drm_w2 = dat.sc.drm_w2;  irm_w2 = dat.sc.irm_w2;
dstb   = dat.sc.dstb;    istb   = dat.sc.istb;
dremap = dat.sc.dremap;  iremap = dat.sc.iremap;
drecnum = dat.sc.drecnum; irecnum = dat.sc.irecnum;
dFieldWidth = dat.sc.dFieldWidth; iFieldWidth = dat.sc.iFieldWidth;
dspn = dat.sc.dspn; ispn = dat.sc.ispn;
% track-insensitive: position tuned in either track, no remapping
dstbonly = any(dstb,2) & ~dremap;
istbonly = any(istb,2) & ~iremap;

%% Plotting
% (from paperFigScripts/paperFigureSingleCell_1w_2w.m, "2 worlds" section)
figDir = fullfile(save2where, 'Sup6');
if ~isfolder(figDir), mkdir(figDir); end

% C: correlation/error matrices for track-sensitive and track-insensitive
iPlot = true;
sortMatrixCorMat(smoothdata(drm_w1(dremap,:),2), smoothdata(drm_w2(dremap,:),2),iPlot);
save_img(gcf,figDir,'Sup6_C_sc_2w_dSPN_reconerrmat_remap'); close(gcf)
sortMatrixCorMat(smoothdata(irm_w1(iremap,:),2), smoothdata(irm_w2(iremap,:),2),iPlot);
save_img(gcf,figDir,'Sup6_C_sc_2w_iSPN_reconerrmat_remap'); close(gcf)
sortMatrixCorMat(smoothdata(drm_w1(dstbonly,:),2), smoothdata(drm_w2(dstbonly,:),2),iPlot);
save_img(gcf,figDir,'Sup6_C_sc_2w_dSPN_reconerrmat_noremap'); close(gcf)
sortMatrixCorMat(smoothdata(irm_w1(istbonly,:),2), smoothdata(irm_w2(istbonly,:),2),iPlot);
save_img(gcf,figDir,'Sup6_C_sc_2w_iSPN_reconerrmat_noremap'); close(gcf)

% per-session correlation matrices (fields with > 10 neurons of that class)
iPlot = false;
possrecs = unique([unique(drecnum);unique(irecnum)]);
emptyOut = struct('maxn',[],'sortind',[],'cm',[],'err',[],'maxn2',[]);
[dout_rmp, dout_stbonly, iout_rmp, iout_stbonly] = deal(repmat(emptyOut, 1, numel(possrecs)));
[dnumremap, dnumstb, inumremap, inumstb] = deal(nan(1, numel(possrecs)));
for nr = 1 : length(possrecs)
    currec = possrecs(nr);
    dt = find(drecnum == currec & dremap);
    dnumremap(nr) = numel(dt);
    if length(dt)>10, dout_rmp(nr) = sortMatrixCorMat(drm_w1(dt,:),drm_w2(dt,:),iPlot); end
    dt = find(drecnum==currec & dstbonly);
    dnumstb(nr) = numel(dt) + dnumremap(nr);
    if length(dt)>10, dout_stbonly(nr) = sortMatrixCorMat(drm_w1(dt,:),drm_w2(dt,:),iPlot); end

    dt = find(irecnum==currec & iremap);
    inumremap(nr) = numel(dt);
    if length(dt)>10, iout_rmp(nr) = sortMatrixCorMat(irm_w1(dt,:),irm_w2(dt,:),iPlot); end
    dt = find(irecnum==currec & istbonly);
    inumstb(nr) = numel(dt) + inumremap(nr);
    if length(dt)>10, iout_stbonly(nr) = sortMatrixCorMat(irm_w1(dt,:),irm_w2(dt,:),iPlot); end
end
dworldcor = cell2mat(arrayfun(@(x) x.rms.worldcor, dspn,'uniformoutput', false)');
iworldcor = cell2mat(arrayfun(@(x) x.rms.worldcor, ispn,'uniformoutput', false)');

% D: correlation (rho) at corresponding positions, per tuning class
getErr = @(o) cell2mat(arrayfun(@(x) nanmean(abs(x.err)), o, 'uniformoutput',0)');
getR2  = @(o) cell2mat(arrayfun(@(x) abs(diag(x.cm)), o, 'UniformOutput',0));
D = struct();
D.sensitive.d = getR2(dout_rmp);      D.sensitive.i = getR2(iout_rmp);
D.insensitive.d = getR2(dout_stbonly); D.insensitive.i = getR2(iout_stbonly);
for c = {'sensitive','insensitive'}
    figure('Position',[100 100 900 700]);hold on;
    shadedErrorBar(1:101,nanmean(D.(c{1}).d'),stderror(D.(c{1}).d'),'lineProps','r')
    shadedErrorBar(1:101,nanmean(D.(c{1}).i'),stderror(D.(c{1}).i'),'lineProps','g')
    title(sprintf('Track-%s: correlation vs position (red dSPN, green iSPN)', c{1}))
    save_img(gcf,figDir,sprintf('Sup6_D_sc_2w_RsquareVSposition_%s', c{1})); close(gcf)
end

% E: cross-track correlation per neuron, track-sensitive vs track-insensitive
E = struct('d',{{dworldcor(dremap), dworldcor(dstbonly)}}, ...
           'i',{{iworldcor(iremap), iworldcor(istbonly)}});
ct_name = {'dSPN','iSPN'};
for c = 1:2
    v = E.(lower(ct_name{c}(1)));
    figure; boxViolin(v{1}, v{2}, get(gcf,'Number'));
    p = ranksum(v{1}, v{2});
    title(sprintf('%s correlation: left = track-sensitive, right = track-insensitive; ranksum p = %s', ...
        ct_name{c}, num2str(p)))
    save_img(gcf,figDir,sprintf('Sup6_E_sc_2w_wcorrRemapVsTune_%s', ct_name{c})); close(gcf)
end

% F: spatial error per session, track-sensitive vs track-insensitive (paired)
F = struct('d',{{getErr(dout_rmp), getErr(dout_stbonly)}}, ...
           'i',{{getErr(iout_rmp), getErr(iout_stbonly)}});
for c = 1:2
    v = F.(lower(ct_name{c}(1)));
    figure; boxScatterLine(v{1}, v{2}, get(gcf,'Number'));
    ok = ~isnan(v{1}) & ~isnan(v{2}); % fields with both classes
    p = signrank(v{1}(ok), v{2}(ok));
    title(sprintf('%s error: left = track-sensitive, right = track-insensitive; signrank p = %s', ...
        ct_name{c}, num2str(p)))
    save_img(gcf,figDir,sprintf('Sup6_F_sc_2w_meanErrRemapVSnoremap_%s_signrank', ct_name{c})); close(gcf)
end

% G: percentage of position-tuned neurons that are track-sensitive, per session
dremapprct = (dnumremap./dnumstb)';
iremapprct = (inumremap./inumstb)';
figure; boxScatterLine(dremapprct(2:end), iremapprct(2:end), get(gcf,'Number'));
pG = signrank(dremapprct(2:end), iremapprct(2:end));
title(['Sign rank test p = ' num2str(pG)]);
save_img(gcf,figDir,'Sup6_G_sc_2w_prctSessRemapDandI_signrank'); close(gcf)

% H: distribution of field centres
[~,dpeak(:,1)] = nanmax(drm_w1,[],2); % peak position in each track, per cell
[~,dpeak(:,2)] = nanmax(drm_w2,[],2);
[~,ipeak(:,1)] = nanmax(irm_w1,[],2);
[~,ipeak(:,2)] = nanmax(irm_w2,[],2);
H = struct();
H.dSPN.insensitive = [dpeak(dstb(:,1)==1 & ~dremap,1); dpeak(dstb(:,2)==1 & ~dremap,2)];
H.dSPN.sensitive   = [dpeak(dstb(:,1)==1 &  dremap,1); dpeak(dstb(:,2)==1 &  dremap,2)];
H.iSPN.insensitive = [ipeak(istb(:,1)==1 & ~iremap,1); ipeak(istb(:,2)==1 & ~iremap,2)];
H.iSPN.sensitive   = [ipeak(istb(:,1)==1 &  iremap,1); ipeak(istb(:,2)==1 &  iremap,2)];
H_hist = struct(); pH = nan(1,2);
for c = 1:2
    v = H.(ct_name{c});
    figure;
    h1 = histogram(v.insensitive,10,'binLimits',[1 100],'Normalization','Probability', ...
        'DisplayStyle','Stairs','EdgeColor',[0.2 0.4 0.8],'LineWidth',2);
    hold on;
    h2 = histogram(v.sensitive,10,'binLimits',[1 100],'Normalization','Probability',...
        'DisplayStyle','Stairs','EdgeColor',[0.8 0.6 0.2],'LineWidth',2);
    H_hist.(ct_name{c}) = [renamevars(sourceDataHist(h1),'count','probability_insensitive'), ...
        renamevars(sourceDataHist(h2),{'bin_left','bin_right','count'}, ...
        {'bin_left_sens','bin_right_sens','probability_sensitive'})]; % read bars before closing
    legend('track-insensitive', 'track-sensitive')
    [~, pH(c)] = kstest2(v.insensitive, v.sensitive);
    title([ct_name{c} ' Ks-test p = ' num2str(pH(c))])
    save_img(gcf,figDir,sprintf('Sup6_H_sc_2w_histogramPeakLoc%s', ct_name{c})); close(gcf)
end

% I: mean activity per neuron, track-sensitive vs track-insensitive
[avg_dspn, ~] = get_singleCell_sess_dff(dspn, 'allrm', 1:100);
[avg_ispn, ~] = get_singleCell_sess_dff(ispn, 'allrm', 1:100);
dtnrmp = any(dstb,2) & ~dremap;
itnrmp = any(istb,2) & ~iremap;
I = struct('d',{{avg_dspn(dremap)', avg_dspn(dtnrmp)'}}, ...
           'i',{{avg_ispn(iremap)', avg_ispn(itnrmp)'}});
for c = 1:2
    v = I.(lower(ct_name{c}(1)));
    figure; boxViolin(v{1}, v{2}, get(gcf,'Number'));
    p = ranksum(v{1}, v{2});
    title(sprintf('%s activity: left = track-sensitive, right = track-insensitive; ranksum p = %s', ...
        ct_name{c}, num2str(p)))
    save_img(gcf,figDir,sprintf('Sup6_I_2w_%s_scActivity_comparison', ct_name{c})); close(gcf)
end

% J: field width, dSPN vs iSPN within each tuning class
J = struct();
J.sensitive.d = [dFieldWidth(dstb(:,1)==1 &  dremap,1); dFieldWidth(dstb(:,2)==1 &  dremap,2)];
J.sensitive.i = [iFieldWidth(istb(:,1)==1 &  iremap,1); iFieldWidth(istb(:,2)==1 &  iremap,2)];
J.insensitive.d = [dFieldWidth(dstb(:,1)==1 & ~dremap,1); dFieldWidth(dstb(:,2)==1 & ~dremap,2)];
J.insensitive.i = [iFieldWidth(istb(:,1)==1 & ~iremap,1); iFieldWidth(istb(:,2)==1 & ~iremap,2)];
for c = {'sensitive','insensitive'}
    figure; boxViolin(J.(c{1}).d, J.(c{1}).i, get(gcf,'Number'));
    p = ranksum(J.(c{1}).d, J.(c{1}).i);
    title(sprintf('Track-%s: left = dSPN, right = iSPN; ranksum p = %s', c{1}, num2str(p)))
    save_img(gcf,figDir,sprintf('Sup6_J_sc_2w_%s_DvsI_boxplotFieldWidth', c{1})); close(gcf)
end
close all

%% Count sizes & report stats -> Excel tab 'Sup6'
% C: only the number of neurons; E-J: every plotted point and its test
fig = 'Sup6';
n_trial = countSampleSize(dat.tbl(dat.trial_selection,:), 'trial');
R = statsRow('Figure',fig,'Panel','C-J','Measure','sample size (2 world, all trials)','nStruct',n_trial);

cellcounts = {'C','dSPN track-sensitive (remapping)', sum(dremap)
              'C','iSPN track-sensitive (remapping)', sum(iremap)
              'C','dSPN track-insensitive (tuned, no remapping)', sum(dstbonly)
              'C','iSPN track-insensitive (tuned, no remapping)', sum(istbonly)};
for r = 1:size(cellcounts,1)
    R = [R; statsRow('Figure',fig,'Panel',cellcounts{r,1},'Measure','number of neurons in correlation matrix', ...
            'Group',cellcounts{r,2},'N',cellcounts{r,3},'N_unit',"neurons",'N_neurons',cellcounts{r,3})]; %#ok<AGROW>
end
% D: number of fields contributing to each averaged curve
for c = {'sensitive','insensitive'}
    for ct = 1:2
        v = D.(c{1}).(lower(ct_name{ct}(1)));
        R = [R; statsRow('Figure',fig,'Panel','D', ...
                'Measure',sprintf('fields averaged, track-%s', c{1}), 'Group',ct_name{ct}, ...
                'N',size(v,2),'N_unit',"fields")]; %#ok<AGROW>
    end
end
% E-J: the tests drawn on the panels
for c = 1:2
    R = [R; reportGroupCompare({'track-sensitive','track-insensitive'}, E.(lower(ct_name{c}(1))), ...
            'Figure',fig,'Panel','E','Measure',sprintf('%s cross-track correlation per neuron', ct_name{c}))
         reportGroupCompare({'track-sensitive','track-insensitive'}, F.(lower(ct_name{c}(1))), ...
            'paired',true, 'Figure',fig,'Panel','F', ...
            'Measure',sprintf('%s spatial error per session (%% of track)', ct_name{c}))
         reportGroupCompare({'track-sensitive','track-insensitive'}, I.(lower(ct_name{c}(1))), ...
            'Figure',fig,'Panel','I','Measure',sprintf('%s mean dF/F per neuron', ct_name{c}))
         statsRow('Figure',fig,'Panel','H','Measure',sprintf('%s field centre distribution', ct_name{c}), ...
            'Group','track-sensitive vs track-insensitive','Test',"Kolmogorov-Smirnov", ...
            'pValue',pH(c),'N',numel(H.(ct_name{c}).sensitive) + numel(H.(ct_name{c}).insensitive), ...
            'N_unit',"neurons (both classes)")]; %#ok<AGROW>
end
R = [R; reportGroupCompare({'dSPN','iSPN'}, {dremapprct(2:end), iremapprct(2:end)}, 'paired',true, ...
        'Figure',fig,'Panel','G','Measure','fraction of position-tuned neurons that are track-sensitive, per session')];
for c = {'sensitive','insensitive'}
    R = [R; reportGroupCompare({'dSPN','iSPN'}, {J.(c{1}).d, J.(c{1}).i}, ...
            'Figure',fig,'Panel','J','Measure',sprintf('field width, track-%s (%% of track)', c{1}))]; %#ok<AGROW>
end

% plotted values (no data points for the C matrices)
blocks = {'Summary: sample sizes; C neuron counts; E-J tests', R
    'Sup6D: correlation vs position, track-sensitive (mean, SEM across fields)', ...
        table((1:101)', nanmean(D.sensitive.d')', stderror(D.sensitive.d')', ...
              nanmean(D.sensitive.i')', stderror(D.sensitive.i')', ...
              'VariableNames',{'position_bin','dSPN_mean','dSPN_SEM','iSPN_mean','iSPN_SEM'})
    'Sup6D: correlation vs position, track-insensitive (mean, SEM across fields)', ...
        table((1:101)', nanmean(D.insensitive.d')', stderror(D.insensitive.d')', ...
              nanmean(D.insensitive.i')', stderror(D.insensitive.i')', ...
              'VariableNames',{'position_bin','dSPN_mean','dSPN_SEM','iSPN_mean','iSPN_SEM'})};
for c = 1:2
    blocks = [blocks
        {sprintf('Sup6E: %s cross-track correlation, one point per neuron', ct_name{c}), ...
            sourceDataPoints({'track-sensitive','track-insensitive'}, E.(lower(ct_name{c}(1))), 'valueName','rho')
         sprintf('Sup6F: %s spatial error, one point per field (paired)', ct_name{c}), ...
            sourceDataPoints({'track-sensitive','track-insensitive'}, F.(lower(ct_name{c}(1))), ...
                'paired',true, 'valueName','error_pct_of_track')
         sprintf('Sup6H: %s field centre histogram (probability per bin)', ct_name{c}), H_hist.(ct_name{c})
         sprintf('Sup6I: %s mean dF/F, one point per neuron', ct_name{c}), ...
            sourceDataPoints({'track-sensitive','track-insensitive'}, I.(lower(ct_name{c}(1))), 'valueName','mean_dFF')}]; %#ok<AGROW>
end
blocks = [blocks
    {'Sup6G: fraction track-sensitive, one point per session (paired)', ...
        sourceDataPoints({'dSPN','iSPN'}, {dremapprct(2:end), iremapprct(2:end)}, ...
            'paired',true, 'valueName','fraction_track_sensitive')
     'Sup6J: field width, track-sensitive, one point per neuron', ...
        sourceDataPoints({'dSPN','iSPN'}, {J.sensitive.d, J.sensitive.i}, 'valueName','field_width_pct')
     'Sup6J: field width, track-insensitive, one point per neuron', ...
        sourceDataPoints({'dSPN','iSPN'}, {J.insensitive.d, J.insensitive.i}, 'valueName','field_width_pct')
     'Sample size per mouse: 2 world trials', n_trial.perMouse}];
writeSourceData(statsFile, fig, blocks);

end
