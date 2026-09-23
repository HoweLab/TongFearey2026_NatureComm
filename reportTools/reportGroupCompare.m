function rows = reportGroupCompare(groups, data, varargin)
% REPORTGROUPCOMPARE  Report rows for a box/violin comparison of groups
% (e.g. boxViolin_colorSig), in the statsRow schema.
%   rows = reportGroupCompare({'dSPN','iSPN'}, {rho_d, rho_i}, ...
%              'pvals', {p_d, p_i}, 'pThresh', 0.001, ...
%              'Figure','Fig1', 'Panel','G', 'Measure','Spearman rho', ...
%              'n', {n_d, n_i})
%
% Output rows:
%   per group : median (IQR, mean +- SEM in Note); signrank vs 0 only if
%               'testVsZero' is true (not shown in the figures)
%   per group : % of elements with individual p < pThresh (if 'pvals'
%               given), with % positive / negative among the significant.
%               Denominator = all elements, as in boxViolin_colorSig.
%   between   : ranksum for 2 groups; pairwise ranksum with Holm
%               correction for >2 groups
%
% Name-value options:
%   Figure, Panel, Measure, Note - text copied to every row
%   pvals   - cell of per-element p-values, same size as data
%   pThresh - threshold for per-element significance (default 0.05)
%   n       - cell of sample-size structs, one per group (or one struct
%             used for all groups)
%   testVsZero - add a Wilcoxon signed-rank test vs 0 per group (default false)

ip = inputParser;
ip.addParameter('Figure', "");
ip.addParameter('Panel', "");
ip.addParameter('Measure', "");
ip.addParameter('Note', "");
ip.addParameter('pvals', {});
ip.addParameter('pThresh', 0.05);
ip.addParameter('n', {});
ip.addParameter('testVsZero', false);
ip.addParameter('paired', false); % paired points: signed-rank instead of rank-sum
ip.parse(varargin{:});
o = ip.Results;

groups = cellstr(groups);
nG = numel(groups);
nStructs = o.n;
if isstruct(nStructs), nStructs = repmat({nStructs}, 1, nG); end
if isempty(nStructs), nStructs = cell(1, nG); end
common = {'Figure',o.Figure, 'Panel',o.Panel};

rows = table();
for g = 1:nG
    x = data{g}(:);
    xv = x(~isnan(x));
    nx = numel(xv);
    q = prctile(xv, [25 50 75]);
    test = ""; pSR = NaN; statSR = NaN; statName = "";
    if o.testVsZero && nx > 0
        [pSR, ~, st] = signrank(xv);
        statSR = st.signedrank;
        test = "Wilcoxon signed-rank vs 0"; statName = "signed rank";
    end
    rows = [rows; statsRow(common{:}, 'Measure', o.Measure, 'Group', groups{g}, ...
        'Test', test, ...
        'Estimate', q(2), 'EstimateType', "median (box plot)", ...
        'StatName', statName, 'Stat', statSR, 'pValue', pSR, ...
        'nStruct', nStructs{g}, 'N', nx, ...
        'Note', strtrim(sprintf('IQR = [%.4g, %.4g]; mean +- SEM = %.4g +- %.4g; %d NaN excluded. %s', ...
        q(1), q(3), mean(xv), std(xv) / sqrt(nx), numel(x) - nx, o.Note)))]; %#ok<AGROW>

    if ~isempty(o.pvals)
        pv = o.pvals{g}(:);
        isSig = pv < o.pThresh;
        nSig = sum(isSig);
        nPos = sum(x(isSig) > 0);
        nNeg = sum(x(isSig) < 0);
        rows = [rows; statsRow(common{:}, ...
            'Measure', sprintf('%s: fraction individually significant (p<%g)', o.Measure, o.pThresh), ...
            'Group', groups{g}, 'Test', "per-element test (count)", ...
            'Estimate', 100 * nSig / numel(x), 'EstimateType', "% significant", ...
            'StatName', "n significant", 'Stat', nSig, ...
            'nStruct', nStructs{g}, 'N', numel(x), ...
            'Note', sprintf('%d/%d significant; of these %.1f%% positive (%d), %.1f%% negative (%d)', ...
            nSig, numel(x), 100 * nPos / nSig, nPos, 100 * nNeg / nSig, nNeg))]; %#ok<AGROW>
    end
end

% between-group comparisons
pairs = nchoosek(1:nG, 2);
pBetween = nan(size(pairs,1), 1);
stat = nan(size(pairs,1), 1);
z = nan(size(pairs,1), 1);
for k = 1:size(pairs,1)
    a = data{pairs(k,1)}; b = data{pairs(k,2)};
    if o.paired % same sessions / fields in both groups
        ok = ~isnan(a(:)) & ~isnan(b(:));
        [pBetween(k), ~, st] = signrank(a(ok), b(ok));
        stat(k) = st.signedrank;
    else
        [pBetween(k), ~, st] = ranksum(a(~isnan(a)), b(~isnan(b)));
        stat(k) = st.ranksum;
    end
    if isfield(st, 'zval'), z(k) = st.zval; end
end
if size(pairs,1) > 1
    pAdj = bonf_holm(pBetween, 0.05);
    pType = "Holm-adjusted across pairs";
else
    pAdj = pBetween;
    pType = "";
end
if o.paired
    testName = "Wilcoxon signed-rank (paired)"; statName = "signed rank";
else
    testName = "Wilcoxon rank-sum"; statName = "rank sum W";
end
for k = 1:size(pairs,1)
    a = data{pairs(k,1)}; b = data{pairs(k,2)};
    rows = [rows; statsRow(common{:}, 'Measure', o.Measure, ...
        'Group', sprintf('%s vs %s', groups{pairs(k,1)}, groups{pairs(k,2)}), ...
        'Test', testName, ...
        'Estimate', median(a,'omitnan') - median(b,'omitnan'), ...
        'EstimateType', "difference of medians", ...
        'StatName', statName, 'Stat', stat(k), ...
        'pValue', pAdj(k), 'pType', pType, ...
        'Note', strtrim(sprintf('n = %d vs %d; z = %.3g. %s', ...
        sum(~isnan(a)), sum(~isnan(b)), z(k), o.Note)))]; %#ok<AGROW>
end
end
