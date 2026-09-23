function match2figs_wSubplots(figs)
% MATCH2FIGS_WSUBPLOTS Match y-limits of corresponding subplots across figures.
%
% INPUT
%   figs : n-by-m cell array
%          each element can be:
%             - a figure handle
%             - or a 1x1 cell containing a figure handle
%
% BEHAVIOR
%   For each row of figs:
%       - find all plotting axes in each figure
%       - sort them in subplot order
%       - check that all figures in that row have the same number of axes
%       - match y-limits for corresponding axes across figures
%
% NOTES
%   - Legends and colorbars are ignored
%   - Supports yyaxis left/right

    [nRow, nCol] = size(figs);
    ax_cell = cell(nRow, nCol);

    % ---------- collect and sort axes for each figure ----------
    for r = 1:nRow
        for c = 1:nCol
            f = figs{r,c};

            % allow either direct figure handle or nested {figHandle}
            if iscell(f)
                f = f{1};
            end

            if ~ishandle(f) || ~strcmp(get(f, 'Type'), 'figure')
                error('figs{%d,%d} is not a valid figure handle.', r, c);
            end

            ax = findall(f, 'Type', 'axes');

            % remove legends / colorbars
            is_good = arrayfun(@(x) ...
                ~isa(x, 'matlab.graphics.illustration.Legend') && ...
                ~isa(x, 'matlab.graphics.illustration.ColorBar'), ax);
            ax = ax(is_good);

            % sort axes in subplot order: top-to-bottom, then left-to-right
            ax = sort_axes_subplot_order(ax);

            ax_cell{r,c} = ax;
        end
    end

    % ---------- match corresponding axes within each row ----------
    for r = 1:nRow
        nAx_ref = numel(ax_cell{r,1});

        for c = 2:nCol
            if numel(ax_cell{r,c}) ~= nAx_ref
                error(['Row %d does not have the same number of subplot axes ' ...
                       'across figures.'], r);
            end
        end

        % for each subplot index, match across all figures in this row
        for k = 1:nAx_ref
            all_axes = gobjects(1, nCol);
            for c = 1:nCol
                all_axes(c) = ax_cell{r,c}(k);
            end
            match_yaxis(all_axes);
        end
    end
end


function ax_sorted = sort_axes_subplot_order(ax)
% Sort axes by subplot position:
% top row first, and within each row left to right

    if isempty(ax)
        ax_sorted = ax;
        return
    end

    pos = vertcat(ax.Position);   % [left bottom width height]
    left = pos(:,1);
    bottom = pos(:,2);

    % sort by descending bottom (top-to-bottom), then ascending left
    sortMat = [-bottom, left];
    [~, idx] = sortrows(sortMat, [1 2]);

    ax_sorted = ax(idx);
end


function match_yaxis(all_axes)
% Match y-limits across a set of axes.
% Supports both normal axes and yyaxis left/right.

    if isempty(all_axes)
        return
    end

    % number of YAxis objects must match across all axes
    nY = numel(all_axes(1).YAxis);
    for k = 2:numel(all_axes)
        if numel(all_axes(k).YAxis) ~= nY
            error('Axes do not have the same number of YAxis objects.');
        end
    end

    % initialize from first axis
    ylimits = zeros(nY, 2);
    for i = 1:nY
        ylimits(i,:) = all_axes(1).YAxis(i).Limits;
    end

    % first pass: find overall min/max for each y-axis
    for ax = all_axes
        for i = 1:nY
            curlims = ax.YAxis(i).Limits;
            ylimits(i,1) = min(ylimits(i,1), curlims(1));
            ylimits(i,2) = max(ylimits(i,2), curlims(2));
        end
    end

    % second pass: assign matched limits
    for ax = all_axes
        for i = 1:nY
            ax.YAxis(i).Limits = ylimits(i,:);
        end
    end
end