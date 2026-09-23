function [bin_out,bin_out_struct] = addSingleMouseLinePlot(curaxs, vuBin_tbl, binName, var2avg,name_id, ...
    ifsingletrial, numPattern, varargin)

% Average within mouse/session/trial/bin/celltype first
ip = inputParser;
ip.addParameter('trial_unit', 'trial') % bout for GUI
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

grp_mouse = groupsummary( vuBin_tbl, ...
    {'mouse','sess','field', trial_unit, binName,'isD','isI'},  "mean", var2avg);
allmice = unique(grp_mouse.mouse);
allbins = unique(grp_mouse.(binName));
allbins = allbins(~isnan(allbins));
if isempty(curaxs)
    f1 = figure(randi(1000000)); clf
    curaxs = axes('Parent', f1);   
end
hold(curaxs, 'on');
% -------- plot one line per mouse on the provided axes --------
bin_out = cell(1, numel(allmice));
bin_out_struct = struct;

% distinct colors
colors = distinct_colors(numPattern);

% optional line-style variation
lineStyles = {'-', '--', ':', '-.'};

plotIdx = 0;

for i = 1:numel(allmice)

    curDat = grp_mouse( ...
        grp_mouse.mouse == allmice(i) & ...
        grp_mouse.isD == 1 & ...
        ~isnan(grp_mouse.(binName)), :);

    if isempty(curDat)
        continue
    end

    plotIdx = plotIdx + 1;

    % average across trial/session within each bin for this mouse
    bin_out{i} = bin_data_colvec( ...
        curDat.(['mean_' var2avg]), ...
        curDat.(binName), ...
        allbins);

    bin_out_struct.(string(allmice(i))) = bin_out{i};

    % cycle color and line style
    colorIdx = mod(plotIdx - 1, numPattern) + 1;
    styleIdx = mod(floor((plotIdx - 1) / numPattern), numel(lineStyles)) + 1;

    plot(curaxs, allbins, mean(bin_out{i}, 1, 'omitnan'), ...
        'Color', colors(colorIdx,:), ...
        'LineStyle', lineStyles{styleIdx}, ...
        'LineWidth', 1.5, ...
        'DisplayName', char(string(allmice(i))));
end
xlim([min(allbins), max(allbins)])    

legend(curaxs, 'show', 'Interpreter', 'none');
% plot cross animal average
if isempty(curaxs)
    u_all = mean(cell2mat(bin_out'), 1 , 'omitnan');
    plot(curaxs, allbins, u_all,'LineWidth',2,'Color','blue');
    
end
% -------- set y-limits on the correct axes --------
allvel = grp_mouse(~isnan(grp_mouse.(binName)), :).(['mean_' var2avg]);
ylim(curaxs, [min(allvel, [], 'omitnan'), max(allvel, [], 'omitnan')]);
title(['mouse mean_' var2avg '_' name_id], 'Interpreter', 'none')
% -------- optional: plot individual trials per mouse --------
if ifsingletrial

    nMouse = numel(allmice);
    ncol = ceil(sqrt(nMouse));
    nrow = ceil(nMouse / ncol);

    % Figure 1: velocity
    figure(randi(1000000)); clf
    sgtitle('Single-mouse velocity trials', 'Interpreter', 'none');

    for i = 1:nMouse
        curDat = grp_mouse( ...
            grp_mouse.mouse == allmice(i) & ...
            grp_mouse.isD == 1 & ...
            ~isnan(grp_mouse.(binName)), :);

        subplot(nrow, ncol, i);
        hold on

        if ~isempty(curDat)
            % one row per trial
            bin_v = bin_data_colvec(curDat.(['mean_' var2avg]), curDat.(binName), allbins);
            plot(allbins, bin_v', 'Color', [0.7 0.7 0.7]);
            plot(allbins,mean(bin_v,1,'omitnan'),'LineWidth',2,'Color','blue')
        end
        xline(0, '--')
        title(string(allmice(i)), 'Interpreter', 'none');
        xlabel(binName, 'Interpreter', 'none');
        ylabel(var2avg);
        xlim([min(allbins), max(allbins)])
    end
    sgtitle([var2avg '_' name_id], 'Interpreter', 'none')
end
end

function cmap = distinct_colors(n)
% Generate n visually distinct RGB colors for categorical line plots.

    if n <= 0
        cmap = [];
        return
    end

    % evenly spaced hues
    hues = linspace(0, 1, n + 1)';
    hues(end) = [];

    % reorder hues so adjacent colors are more separated
    step = max(1, floor(n / 2));
    idx = mod((0:n-1) * step, n) + 1;

    % if the step causes repeats for even n, fall back to simple ordering
    if numel(unique(idx)) < n
        idx = 1:n;
    end

    hues = hues(idx);

    sat = 0.75 * ones(n, 1);
    val = 0.85 * ones(n, 1);

    cmap = hsv2rgb([hues, sat, val]);
end