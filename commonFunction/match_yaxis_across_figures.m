function match_yaxis_across_figures(figCells)
% match_yaxis_across_figures(figCells)
%
% Input:
%   figCells - 1 x n cell array of figure handles
%
% Purpose:
%   If each figure has the same number/order of subplots, this function
%   matches the y-axis limits of corresponding subplots across figures.
%
% Notes:
%   - Supports normal y-axis and yyaxis left/right.
%   - Axes are matched by their plotting order within each figure.
%
% Example:
%   match_yaxis_across_figures({fig1, fig2, fig3})

    if ~iscell(figCells)
        error('Input must be a cell array of figure handles.');
    end

    nFig = numel(figCells);

    % ------------------------------------------------------------
    % Collect axes from each figure
    % ------------------------------------------------------------
    axesCell = cell(1, nFig);

    for f = 1:nFig
        fig = figCells{f};

        if ~isgraphics(fig, 'figure')
            error('Element %d is not a valid figure handle.', f);
        end

        % Find regular axes only, excluding legends/colorbars if present
        ax = findall(fig, 'Type', 'axes');

        % Remove axes used by legends/colorbars, just in case
        ax = ax(~arrayfun(@(a) isa(a, 'matlab.graphics.illustration.ColorBar') || ...
                             isa(a, 'matlab.graphics.illustration.Legend'), ax));

        % MATLAB returns axes in reverse stacking order, so flip it
        ax = flipud(ax(:));

        axesCell{f} = ax;
    end

    % ------------------------------------------------------------
    % Check that all figures have the same number of axes
    % ------------------------------------------------------------
    nAxes = numel(axesCell{1});

    for f = 2:nFig
        if numel(axesCell{f}) ~= nAxes
            error('All figures must have the same number of subplots/axes.');
        end
    end

    % ------------------------------------------------------------
    % Match corresponding axes across figures
    % ------------------------------------------------------------
    for a = 1:nAxes

        % Number of y-axes, usually 1 or 2 for yyaxis
        nYAxes = numel(axesCell{1}(a).YAxis);

        % Check matching yyaxis structure
        for f = 2:nFig
            if numel(axesCell{f}(a).YAxis) ~= nYAxes
                error('Corresponding subplot %d does not have the same number of y-axes across figures.', a);
            end
        end

        % Initialize target limits for each y-axis
        targetLims = nan(nYAxes, 2);

        for y = 1:nYAxes
            curMin = inf;
            curMax = -inf;

            for f = 1:nFig
                ax = axesCell{f}(a);
                lims = ax.YAxis(y).Limits;

                curMin = min(curMin, lims(1));
                curMax = max(curMax, lims(2));
            end

            targetLims(y, :) = [curMin, curMax];
        end

        % Apply matched limits back to corresponding axes
        for f = 1:nFig
            ax = axesCell{f}(a);

            for y = 1:nYAxes
                ax.YAxis(y).Limits = targetLims(y, :);
            end
        end
    end
end