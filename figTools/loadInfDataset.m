function out = loadInfDataset(varargin)
% LOADINFDATASET  Infinite-track data for Fig 5 and Sup 9. Two datasets are
% needed and they are NOT interchangeable:
%   pop   - all recorded neurons (processedData_all.mat): population panels
%           (Fig5 A, B, F-I; Sup9 E-G)
%   tuned - tuning-subtype file (tuningSubType_..._inf_newtype.mat), whose
%           avgFc3_<type> columns hold only the neurons of each tuning
%           class: subpopulation panels (Fig5 J-N)
%   sc    - single-cell file: rasters and field statistics (Fig5 D-E, Sup9 B-D)
%
%   dat = loadInfDataset();
%   dat = loadInfDataset('loadTuned',false);   % population panels only
%
% out.pop / out.tuned each contain:
%   tbl, trial_selection, bout_selection, bout (bout_descriptive table)
% out.pop also has binType, the bout-distance bin column made by bin_BoutDist.
% out.info holds the per-session tuning counts, out.sc the single-cell
% variables (with drecnum / irecnum and the "tuned to both" flags).
%
% Name-value options (defaults):
%   popPath, fpath (tuned), scpath, mice_include, exclude_mice,
%   min_boutdist (150 cm), min_boutdur (5 s), bout_rg ([0 150]),
%   fr (31.25), binsize (1.5 cm), scVars, loadTuned / loadSc (true)

ip = inputParser;
ip.addParameter('popPath', '');
ip.addParameter('fpath',   '');
ip.addParameter('scpath',  '');
ip.addParameter('mice_include', {'tdt','a2a'});
ip.addParameter('exclude_mice', {'a2a_262'});
ip.addParameter('min_boutdist', 150);  % cm
ip.addParameter('min_boutdur', 5);     % s
ip.addParameter('bout_rg', [0 150]);   % cm, range binned within a bout
ip.addParameter('fr', 31.25);
ip.addParameter('binsize', 1.5);       % cm, gives 100 bins over bout_rg
ip.addParameter('scVars', {});
ip.addParameter('loadTuned', true);
ip.addParameter('loadSc', true);
ip.parse(varargin{:});
o = ip.Results;

out.ext = {'bDist','track','only_bDist','only_track','nta','nta_track'};
out.fr = o.fr; out.binsize = o.binsize; out.min_boutdist = o.min_boutdist;
out.bout_rg = o.bout_rg;

% ---- population data: all recorded neurons ----
conditions = {'sess_type','inf';'mismatch',0};
pop_tbl = mkTblForPlot(o.popPath, conditions, ...
    'trial_fields',{'trialTimes','trial_skipped'}, ...
    'bout_fields',{'bout_rewarded','peak_vel','bout_length','normalized_dist2rew'});
[pop_tbl, out.pop.bout] = bout_descriptive(pop_tbl, o.fr, o.binsize, false);
% bout-distance bins (same call as paperFigurePopulation_Inf_BoutDist.m)
d2rew = [min(pop_tbl.info_normalized_dist2rew), max(pop_tbl.info_normalized_dist2rew)];
[pop_tbl, ~, binType] = bin_BoutDist(pop_tbl, out.pop.bout, o.bout_rg, ...
    'dist2reward', d2rew, 'min_boutdist', o.min_boutdist, ...
    'min_boutdur', o.min_boutdur, 'BinSize', o.binsize);
inMice = contains(cellstr(pop_tbl.mouse), o.mice_include) & ...
    ~contains(cellstr(pop_tbl.mouse), o.exclude_mice);
out.pop.tbl = pop_tbl;
out.pop.binType = char(binType);
out.pop.trial_selection = ~pop_tbl.skip_repeat;
out.pop.bout_selection = pop_tbl.info_bout_length > o.min_boutdur & ...
    pop_tbl.info_peak_vel > 15 & inMice;
% bouts long enough for the triggered averages (Sup9 E-F)
out.pop.trig_selection = pop_tbl.tot_dist >= o.bout_rg(2) & ...
    pop_tbl.info_peak_vel > 15 & pop_tbl.info_bout_length > o.min_boutdur & inMice;
fprintf('Population: %d bouts binned by distance, %d bouts for triggered averages\n', ...
    numel(unique(pop_tbl.ubid(out.pop.bout_selection & out.pop.trial_selection & ...
    ~isnan(pop_tbl.(binType))))), ...
    numel(unique(pop_tbl.ubid(out.pop.trig_selection & out.pop.trial_selection))));

% ---- tuning-subtype data: one column per tuning class ----
if o.loadTuned
    conditions = {'exp_stage',{'inf'}; 'mismatch',0};
    tun_tbl = mkTblForPlot(o.fpath, conditions, ...
        'trial_fields',{'trialTimes','trial_skipped'}, ...
        'bout_fields',{'bout_rewarded','peak_vel','bout_length'});
    [tun_tbl_all, out.tuned.bout] = bout_descriptive(tun_tbl, o.fr, o.binsize, false);
    inMiceT = contains(cellstr(tun_tbl.mouse), o.mice_include) & ...
        ~contains(cellstr(tun_tbl.mouse), o.exclude_mice);
    out.tuned.tbl = tun_tbl;
    out.tuned.trial_selection = ~tun_tbl.skip_repeat & tun_tbl.exp_stage == 'inf' & inMiceT;
    out.tuned.bout_selection = tun_tbl.info_bout_length > 0 & tun_tbl.info_peak_vel > 15 & ...
        inMiceT & tun_tbl_all.tot_dist > o.min_boutdist;
    S = load(o.fpath,'cat_info');
    info_tbl = struct2table(S.cat_info);
    out.info = info_tbl(matches(info_tbl.exp_stage,'inf'), :);
end

% ---- single-cell data ----
if o.loadSc
    scVars = unique([{'dspn','ispn','dstb','istb','drm_w1','drm_w2','irm_w1','irm_w2'}, ...
        o.scVars(:)'], 'stable');
    sc = load(o.scpath, scVars{:});
    sc.drecnum = cellfun(@(x) x.recNum, {sc.dspn.info});
    sc.irecnum = cellfun(@(x) x.recNum, {sc.ispn.info});
    sc.dboth = sum(sc.dstb,2) == 2;   % tuned to both distance and visual position
    sc.iboth = sum(sc.istb,2) == 2;
    out.sc = sc;
end
end
