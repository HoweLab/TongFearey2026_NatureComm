function Sup7(cfg)
% Sup 7: Supplement for main Fig3
% Single cell from 2 track experiments with stats
% A-C: onset/offset triggered dSPN and iSPN for each subpopulation
% D-F: their dSPN - iSPN difference
% K-M: difference of the imbalance between subpopulations (TS-TI, TS-NPT, TI-NPT)
% (G-J are single-mouse examples: not generated here)
% cfg: input/output paths, set in runAllFigures.m

save2where = cfg.outDir;
statsFile = cfg.statsFile; % one workbook, one tab per figure

% Load 2w data (same sessions/mice excluded as Fig3)
fpath  = cfg.tuned2w;
scpath = cfg.sc2w;
Tpath  = cfg.switchInfo;
dat = load2wDataset('fpath', fpath, 'scpath', scpath, 'Tpath', Tpath);
track_2w_tbl = dat.tbl;
trial_selection = dat.trial_selection;
bout_selection = track_2w_tbl.info_bout_intrack == 1; % rewarded and unrewarded bouts
fam_idx = track_2w_tbl.info_trial_novel == 0;         % Track 1, familiar world

xdata = [1:100; linspace(1,100,100)];
fake_info = struct('startLoc',1,'rewardLoc',100);
% Trig window
sr = 31;
on_win = [-31:3*31; linspace(-1,3,125)];
off_win = [3*-31:31; linspace(-3,1,125)];

% subpopulations: track-sensitive (TS), track-insensitive (TI), non-position tuned (NPT)
pops = struct( ...
    'name',  {'TS','TI','NPT'}, ...
    'y',     {'avgFc3_rmp_overlap', 'avgFc3_tnrmp', 'avgFc3_nta'}, ...
    'panel', {{'A','D'}, {'B','E'}, {'C','F'}}, ...
    'label', {'track-sensitive', 'track-insensitive', 'non-position tuned (active)'}, ...
    'count', {{'rmp_w1_d','rmp_w1_i'}, {'tnrmp_d','tnrmp_i'}, {'nta_d','nta_i'}});

%% Plotting
% (from paperFigScripts/paperFigureSubPopulation.m: "triggered average for
%  movement bouts" and "plot significance between DI difference for three types")
figDir = fullfile(save2where, 'Sup7');
if ~isfolder(figDir), mkdir(figDir); end

% A-F: onset / offset triggered activity per subpopulation, Track 1 trials
trig_sel = trial_selection & bout_selection & fam_idx;
trig_pop = struct();
for p = 1:numel(pops)
    pp = pops(p);
    trig_pop.(pp.name) = struct( ...
        'onsets',  getPlotDatFromTbl(track_2w_tbl(trig_sel,:),'onsets','xdata',on_win, 'y',pp.y), ...
        'offsets', getPlotDatFromTbl(track_2w_tbl(trig_sel,:),'offsets','xdata',off_win, 'y',pp.y));
    if p == 1 % velocity is plotted once, under panel A
        figure('Position',[900 100 900 700])
        plotTriggered_new(trig_pop.(pp.name),'traces',{'velocity'},'fromLME',true, ...
            'line_colors',{[0 0 1]},'plot_vel',false,'trig_events',{'onsets','offsets'});
        save_img(gcf,figDir,'Sup7_A_2world_triggeredAverageVelocity'); close(gcf)
    end
    % A, B, C: dSPN and iSPN
    figure('Position',[100 100 900 700])
    plotTriggered_new(trig_pop.(pp.name),'traces',{'dSPN','iSPN'},'fromLME',true, ...
        'plot_vel',false,'trig_events',{'onsets','offsets'});
    sgtitle(sprintf('Sup7%s %s', pp.panel{1}, pp.label),'Interpreter','none')
    save_img(gcf,figDir,sprintf('Sup7_%s_2world_triggeredAverageD&I_%s', pp.panel{1}, pp.name)); close(gcf)
    % D, E, F: dSPN - iSPN, with significance
    figure('Position',[100 100 900 700])
    plotTriggered_new(trig_pop.(pp.name),'traces',{'diff_di'},'fromLME',true, ...
        'plot_sig',true,'line_colors',{[0 0 0]}, ...
        'plot_vel',false,'trig_events',{'onsets','offsets'});
    sgtitle(sprintf('Sup7%s %s', pp.panel{2}, pp.label),'Interpreter','none')
    save_img(gcf,figDir,sprintf('Sup7_%s_2world_triggeredAverageDiff_%s', pp.panel{2}, pp.name)); close(gcf)
end
close all

% K-M: difference of the dSPN-iSPN imbalance between subpopulations
vuBin_tbl = groupsummary(track_2w_tbl(trial_selection,:), ...
    {'mouse','sess','field','trial','vuBin','isD','isI','celltype','info_trial_novel'}, ...
    'mean',{'avgFc3_rmp_overlap','avgFc3_tnrmp','avgFc3_nta','velocity','acceleration'});
fam_bin = vuBin_tbl.info_trial_novel == 0; % Track 1 trials
pop_tbl = struct();
for p = 1:numel(pops)
    t = vuBin_tbl(fam_bin,:);
    t.mean_avgFc3 = t.(['mean_' pops(p).y]);
    pop_tbl.(pops(p).name) = t;
end
cmp = {'K','TS','TI'; 'L','TS','NPT'; 'M','TI','NPT'};
cmp_dat = struct();
for k = 1:size(cmp,1)
    [cmp_dat.(cmp{k,1}), figs] = get_DIdiff_glm_SX(pop_tbl.(cmp{k,2}), pop_tbl.(cmp{k,3}), ...
        'y','mean_avgFc3','field','vuBin', 'plotornot', 1, 'vel_y', []);
    sgtitle(figs.figExp01, sprintf('Sup7%s: %s - %s', cmp{k,1}, cmp{k,2}, cmp{k,3}))
    save_img(figs.figExp01, figDir, sprintf('Sup7_%s_2w_positionBin_DIdiff_%s_vs_%s_sig', cmp{k,:}));
    close all
end

%% Count sizes & report stats -> Excel tab 'Sup7'
fig = 'Sup7';
n_trial = countSampleSize(track_2w_tbl(trial_selection & fam_idx,:), 'trial');   % K-M
n_bout  = countSampleSize(track_2w_tbl(trig_sel,:), 'bout');                     % A-F

R = [statsRow('Figure',fig,'Panel','A-F','Measure','sample size (Track 1 in-track locomotion bouts)','nStruct',n_bout)
     statsRow('Figure',fig,'Panel','K-M','Measure','sample size (Track 1 trials, position-binned)','nStruct',n_trial)];
% number of neurons per subpopulation (sessions used, from load2wDataset)
for p = 1:numel(pops)
    pp = pops(p);
    for c = 1:2
        ct = {'dSPN','iSPN'};
        nCell = sum([dat.info.(pp.count{c})]);
        R = [R; statsRow('Figure',fig,'Panel',[pp.panel{1} ',' pp.panel{2} ',K-M'], ...
                'Measure',sprintf('number of neurons (%s)', pp.label), 'Group',ct{c}, ...
                'N',nCell,'N_unit',"neurons",'N_neurons',nCell,'nStruct',n_bout)]; %#ok<AGROW>
    end
end

% plotted values (A-C: no tests; D-F and K-M: Holm p for the significance bars)
blocks = {'Summary: sample sizes and neurons per subpopulation', R};
for p = 1:numel(pops)
    pp = pops(p);
    for ev = {'onsets','offsets'}
        blocks = [blocks
            {sprintf('Sup7%s: %s, %s-aligned dSPN and iSPN (LME mean, SEM)', pp.panel{1}, pp.label, ev{1}), ...
                sourceDataCurve(trig_pop.(pp.name).(ev{1}), {'dSPN','iSPN'}, 'xName','time_s')
             sprintf('Sup7%s: %s, %s-aligned dSPN - iSPN (LME mean, SEM, Holm p)', pp.panel{2}, pp.label, ev{1}), ...
                sourceDataCurve(trig_pop.(pp.name).(ev{1}), {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, ...
                    'pval',true, 'xName','time_s')}]; %#ok<AGROW>
    end
end
blocks = [blocks
    {'Sup7A: onset-aligned velocity (LME mean, SEM)', ...
        sourceDataCurve(trig_pop.TS.onsets, {'velocity'}, 'xName','time_s')
     'Sup7A: offset-aligned velocity (LME mean, SEM)', ...
        sourceDataCurve(trig_pop.TS.offsets, {'velocity'}, 'xName','time_s')}];
for k = 1:size(cmp,1)
    blocks(end+1,:) = {sprintf('Sup7%s: (dSPN - iSPN) %s - %s by track position (LME mean, SEM, Holm p)', cmp{k,:}), ...
        sourceDataCurve(cmp_dat.(cmp{k,1}), {'diff_DI_exp01'}, ...
            'labels',{sprintf('%s_minus_%s_DIdiff', cmp{k,2}, cmp{k,3})}, 'pval',true, 'xName','position_bin')};
end
blocks = [blocks
    {'Sample size per mouse: Track 1 bouts (A-F)', n_bout.perMouse
     'Sample size per mouse: Track 1 trials (K-M)', n_trial.perMouse}];
writeSourceData(statsFile, fig, blocks);

end
