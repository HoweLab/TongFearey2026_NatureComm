function rows = reportCrossCorr(avgout, varargin)
% REPORTCROSSCORR  Report rows for the population cross-correlation output
% of crosscor_DI_behav (avgout), in the statsRow schema.
%   rows = reportCrossCorr(avgout_1w_onset, 'Figure','Fig2', 'Panel','E middle', ...
%              'Measure','onset-aligned', 'n', n_bout, 'lagScale', 1/31, 'lagUnit','s')
%
% One row per behavior field of avgout (e.g. acceleration, velocity):
%   Estimate = peak normalized cross-correlation (max(c), as avgout.maxlg)
%   X        = lag at the peak, multiplied by lagScale (lag > 0: DI lags behavior)
%   pValue   = Spearman p at that lag (avgout.spear_tho_p), only with 'withP'
%   Note     = c at lag 0 (+ circular-shift permutation result with 'withP')
% N is the number of overlapping bins at the peak lag; the N_* sample-size
% columns (from 'n') describe the data behind the LME mean traces.
%
% Name-value options: Figure, Panel, Measure, Note, n, lagScale (default 1),
% lagUnit (default 'bin').

ip = inputParser;
ip.addParameter('Figure', "");
ip.addParameter('Panel', "");
ip.addParameter('Measure', "");
ip.addParameter('Note', "");
ip.addParameter('n', []);
ip.addParameter('lagScale', 1);
ip.addParameter('lagUnit', "bin");
ip.addParameter('withP', false); % Spearman p at the peak lag (not shown in the figure)
ip.parse(varargin{:});
o = ip.Results;

rows = table();
for b = fieldnames(avgout)'
    a = avgout.(b{1});
    c = a.xcorr_c(:)';
    lags = a.xcorr_lags(:)';
    [cPeak, k] = max(c);
    k0 = find(lags == 0, 1);
    test = ""; p = NaN; note = "";
    if o.withP
        test = "Spearman correlation at peak lag (LME mean traces)";
        p = a.spear_tho_p(k);
        permSig = "n/a";
        if ~isnan(a.perm_p(k)), permSig = string(logical(a.perm_p(k))); end
        note = sprintf('Permutation (circular shift, p<0.05) significant at peak lag: %s; ', permSig);
    end
    rows = [rows; statsRow('Figure',o.Figure, 'Panel',o.Panel, ...
        'Measure', strtrim(sprintf('%s dSPN-iSPN vs %s cross-correlation', o.Measure, b{1})), ...
        'Group', sprintf('DI diff vs %s', b{1}), 'Test', test, ...
        'X', lags(k) * o.lagScale, 'Estimate', cPeak, ...
        'EstimateType', "peak normalized xcorr", ...
        'StatName', "lag at peak (" + string(o.lagUnit) + ")", 'Stat', lags(k) * o.lagScale, ...
        'pValue', p, 'nStruct', o.n, ...
        'N', a.n_overlap(k), 'N_unit', "overlapping bins", ...
        'Note', strtrim(sprintf('%sxcorr at lag 0 = %.3g; lag in samples = %d. %s', ...
        note, c(k0), lags(k), o.Note)))]; %#ok<AGROW>
end
end
