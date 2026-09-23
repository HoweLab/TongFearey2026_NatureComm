function [p, pctSig] = boxViolin_colorSig(pThresh, varargin)
% boxViolin_colorSig(pThresh, data1, p1, data2, p2, ..., figureID)
%
% Example:
%   [p, pctSig] = boxViolin_colorSig(0.05, data1, pvals1, data2, pvals2, 5)
%
% Inputs:
%   pThresh   - significance threshold, e.g. 0.05
%   data1     - data vector for group 1
%   p1        - p-value vector for group 1 (same length as data1)
%   data2     - data vector for group 2
%   p2        - p-value vector for group 2
%   ...
%   figureID  - figure number (last input)
%
% Outputs:
%   p       - ranksum p-value between group 1 and 2 (if n >= 2)
%   pctSig  - percentage significant in each group
%
% Notes:
%   - Significant points (p < pThresh) are colored red
%   - Non-significant points are colored gray
    rng('default') 
    p = NaN;

    fid = varargin{end};
    figure(fid);
    clf;
    varargin(end) = [];

    if mod(length(varargin), 2) ~= 0
        error('Inputs after pThresh must be data/p-value pairs, followed by figureID.');
    end

    n = length(varargin) / 2;
    xp = 1:n;

    boxDat = [];
    boxG   = [];
    pctSig = nan(1, n);
    pctPositive = nan(1,n);
    hold on;

    for i = 1:n
        data  = varargin{2*i - 1};
        pvals = varargin{2*i};

        data  = data(:);
        pvals = pvals(:);

        if length(data) ~= length(pvals)
            error('For group %d, data and p-value vectors must have the same length.', i);
        end

        sigIdx = pvals < pThresh;
        pctSig(i) = 100 * mean(sigIdx);
        sigDat = data(sigIdx);
        posIdx = sigDat > 0;  % positive correlation for sig cors
        pctPositive(i) = 100 * mean(posIdx);

        xt = repelem(xp(i), length(data))';

        ptColor = repmat([0.5 0.5 0.5], length(data), 1);
        ptColor(sigIdx, :) = repmat([1 0 0], sum(sigIdx), 1);

        swarmchart(xt, data, 10, ptColor, 'filled');

        boxDat = [boxDat; data];
        boxG   = [boxG; repelem(i, length(data))'];
    end

    boxplot(boxDat, boxG, 'plotstyle', 'compact');
    delete(findobj(gca, 'Tag', 'Outliers')); % remove outlier symbols
    axis square;
    hold on;

    if n >= 2
        p = ranksum(boxDat(boxG == 1), boxDat(boxG == 2));
        title(sprintf('G1 and G2 significant, p = %.3g', p));
    end

    % get axis limits after plotting everything
    yl = ylim;
    yr = yl(2) - yl(1);

    % optionally add a bit of headroom for labels
    ylim([yl(1), yl(2) + 0.08 * yr]);
    yl = ylim;
    yr = yl(2) - yl(1);

    % place percentage labels above each group
    yText = yl(2) - 0.04 * yr;
    for i = 1:n
        text(xp(i), yText, sprintf('%.1f%% sig; %.1f%% sig pos', pctSig(i), pctPositive(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontSize', 10);
    end
end