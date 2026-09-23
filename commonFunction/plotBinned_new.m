function plotBinned_new(binned_data, trial_info, varargin)
% fig = plotBinnedData(binned_data, trial_info ...)
% There are two fixed input:
%   binned_data: is a struct with each field in the struct is a timeseries
%   data binned by location. Fields in the struct are m*n where m is the
%   number of bins and n is the number of trials. (IMPORTANT: must contain
%   field velocityBin
%   trial_info: is a struct with trial by trial infomation such as trial
%   time, world skip etc.
%   traces: a cell array specify which traces to plot except velocity
%   (which will be plotted on yyaxis right)


% parse input to setup optional input.
ip = inputParser;
ip.addParameter('traces',{'avgDSPNBin','avgISPNBin'})
ip.addParameter('left_traces',{'velocity'})
ip.addParameter('fig',gcf); % plot on current figure unless specified
ip.addParameter('plot_idx',[]); % any specified index to plot, assumed to be boolean index
ip.addParameter('compare_idx',[]);
ip.addParameter('x_label','Track Position (virmen Units)');
ip.addParameter('y_label','dFF');
ip.addParameter('left_y_label','Velocity (m/s)')
ip.addParameter('line_colors', {});
ip.addParameter('plot_title','');
ip.addParameter('plot_vel', false);
ip.addParameter('fromLME',false);
ip.addParameter('plot_sig',false);
ip.addParameter('ax_Tag','traces')
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

if isfield(trial_info,'startLoc')
    startLoc = trial_info.startLoc;
    rewardLoc = trial_info.rewardLoc;
else
    startLoc = -65;
    rewardLoc = 92;
end

% skip any trials that were incomplete (specified as trial_skipped in
% trial_info).
if ~fromLME
    ind = ~trial_info.trial_skipped;
    % see if additional index is given
    if isempty(plot_idx)
        plot_idx = true(size(binned_data.(traces{1}),1),1);
    end
    if ~islogical(plot_idx)
        temp = plot_idx;
        plot_idx = false(size(binned_data.(traces{1}),1),1);
        plot_idx(temp) = true;
        clear temp
    end
    ind = ind & plot_idx;

end

% check if there exist tiledlayout on current figure, if so, makesure it's
% on 'flow' and add new tiles to it.
if ~isempty(get(fig,'Children'))
    if isa(fig.Children,'matlab.graphics.layout.TiledChartLayout') && matches(fig.Children.TileArrangement,'flow')
        tile = fig.Children;
    else
        tile = tiledlayout(fig,'flow');
    end
else
    tile = tiledlayout(fig,'flow');
end

% By default I will plot two subplots for each set of input.
% first, plot each traces seperately, and then plot the one-by-one
% difference for each pair of traces given.
plot_diff = false;
plot_comp = false;
ax = nexttile(tile); % initiate the subplots first
ax.Tag = ax_Tag;
if numel(traces) == 2 & ~fromLME
    plot_diff = true;
end
% another condition we might plot different is to plot compare, just do a
% separate check and wrap plotting in an if statement later
if (~isempty(compare_idx) && numel(traces) == 1)
    plot_comp = true;
end

% before plotting, do some setting first.

if isempty(line_colors)
    line_colors = {
        [0.8, 0.2, 0.2]; % Red
        [0.2, 0.8, 0.2]; % Green
        [0.8, 0.8, 0];   % Yellow
        [0.2, 0.8, 0.8]; % Cyan
        [0.8, 0.2, 0.8]; % Magenta
        [0.8, 0.4, 0];   % Orange
        [0.8, 0.6, 0.64] % Pink
        };
end
num_c = numel(line_colors);
hold(ax,'on')
lines = gobjects(1,numel(traces));
line_names = replace(traces,'Bin','');
y_lims = [0,0];
% loop through each traces to plot, get traces and plot them on ax.
for t = 1:numel(traces)
    if fromLME
        xdata = binned_data.xdata;
        mu = binned_data.(traces{t}).mu;
        sem = binned_data.(traces{t}).sem;
    else
        binSize = binned_data.edges(2) - binned_data.edges(1);
        xdata = binned_data.edges(1:end-1) + 0.5*binSize;
        data = binned_data.(traces{t})(ind,:);
        mu = mean(data,1,'omitnan');
        sem = std(data,1,"omitnan")/sqrt(size(data,1));
    end
    plot_bins = find(xdata >= startLoc & xdata <= rewardLoc);
    xdata = xdata(plot_bins);
    mu = mu(plot_bins);
    sem = sem(plot_bins);
    y_lims(1) = min(y_lims(1),min(mu-sem));
    y_lims(2) = max(y_lims(2),max(mu+sem));
    h = shadedErrorBar(xdata, mu, sem,...
        'lineProps',{'-','Color',line_colors{rem(t-1,num_c)+1}},...
        'plotAxes',ax);
    lines(t) = h.mainLine;
    %    yline(0,'Color','r')
    axis(ax,'tight')
    ylabel(ax,y_label)
    xlabel(ax,x_label)
    if plot_sig && isfield(binned_data.(traces{t}),'pval')
        pval = binned_data.(traces{t}).pval(plot_bins);
        if ~isempty(xdata(pval<0.05))
            plot(xdata(pval < 0.05),0,'*', ...
                'Color',line_colors{rem(t-1,num_c)+1}, ...
                'tag','sigMarkers');
            y_lims(1) = min(y_lims(1),min(mu-sem)) * 1.05;
            y_lims(2) = max(y_lims(2),max(mu+sem)) * 1.05;
        end
    end
end
yline(ax,0,'k--')
ax.XLim = [min(xdata) max(xdata)];
ax.YLim = y_lims;
mks = findobj(ax,'Tag','sigMarkers');
arrayfun(@(x) set(x,'YData', y_lims(2) * 0.95),mks)


% plot velocity for this subplot
if plot_vel

    vel_ylims = [inf -inf];

    yyaxis(ax, 'right');

    for t = 1:numel(left_traces)

        if fromLME
            vel_mu  = binned_data.(left_traces{t}).mu;
            vel_sem = binned_data.(left_traces{t}).sem;
        else
            velBin  = binned_data.velocityBin(ind,:);
            vel_mu  = mean(velBin, 1, 'omitnan');
            vel_sem = std(velBin, 1, 'omitnan') ./ sqrt(size(velBin,1));
        end

        h = shadedErrorBar( ...
            xdata, ...
            vel_mu(plot_bins), ...
            vel_sem(plot_bins), ...
            'lineProps', 'b', ...
            'plotAxes', ax);

        lines(end+1) = h.mainLine;
        line_names(end+1) = left_traces(t);

        % collect velocity data limits only within plotted bins
        cur_low  = min(vel_mu(plot_bins) - vel_sem(plot_bins), [], 'omitnan');
        cur_high = max(vel_mu(plot_bins) + vel_sem(plot_bins), [], 'omitnan');

        vel_ylims(1) = min(vel_ylims(1), cur_low);
        vel_ylims(2) = max(vel_ylims(2), cur_high);
    end

    ylabel(ax, left_y_label)
    set(ax, 'YColor', 'b')

    % ------------------------------------------------------------
    % Force left and right yyaxis zero to align
    % ------------------------------------------------------------

    % get left axis limits
    yyaxis(ax, 'left');
    left_ylims = ylim(ax);

    % make sure left axis includes zero
    left_ylims(1) = min(left_ylims(1), 0);
    left_ylims(2) = max(left_ylims(2), 0);

    % optional padding on left axis
    left_pad = 0.05 * range(left_ylims);
    left_ylims = [left_ylims(1) - left_pad, left_ylims(2) + left_pad];

    % zero position on left axis, as a fraction from bottom to top
    zero_frac = (0 - left_ylims(1)) / range(left_ylims);

    % make sure right axis includes zero
    vel_ylims(1) = min(vel_ylims(1), 0);
    vel_ylims(2) = max(vel_ylims(2), 0);

    % add padding to velocity data range
    vel_pad = 0.05 * range(vel_ylims);
    vel_ylims = [vel_ylims(1) - vel_pad, vel_ylims(2) + vel_pad];

    % data requirements
    right_min_data = vel_ylims(1);
    right_max_data = vel_ylims(2);

    % choose right-axis range needed so zero lands at same zero_frac
    range_needed_lower = abs(right_min_data) / zero_frac;
    range_needed_upper = abs(right_max_data) / (1 - zero_frac);

    right_range = max(range_needed_lower, range_needed_upper);

    % construct right limits with aligned zero
    right_ylims = [
        -zero_frac * right_range, ...
         (1 - zero_frac) * right_range
    ];

    % apply final limits
    yyaxis(ax, 'left');
    ylim(ax, left_ylims);

    yyaxis(ax, 'right');
    ylim(ax, right_ylims);
end

line_names = replace(line_names, '_', ' ');
legend(ax, lines, line_names, 'Location', 'northeastoutside')
title(ax, plot_title)




% if there are only two traces do t-test and plot signifcants.
if plot_diff
    % save binned_data in variable here for easier reference.
    data_a = binned_data.(traces{1})(ind,:);
    data_b = binned_data.(traces{2})(ind,:);
    
    % make data struct for the diff traces
    diff = data_a - data_b;
    diff_dat.diff.mu = mean(diff,3,'omitnan')';
    diff_dat.diff.sem = std(diff,0,3,'omitnan')'/sqrt(size(diff,1));
    [~,pval,~] = ttest2(data_a,data_b,'dim',3,'tail','both');
    diff_dat.diff.pval = pval;
    % call self to add a new tile.
    plotBinned_new(diff_dat,'traces',{'diff'},'fromLME',true,'fig',fig, ...
            'plot_sig',true,'line_colors',{[0 0 0]}, ...
            'plot_vel',plot_vel,'left_traces',left_traces,'ax_Tag','diff');
elseif plot_comp 
    % this is bad coding but I will just clear the axes and replot if we 
    % want to compare trials for a single field.
    cla(ax,'reset')
    hold(ax,'on')
    % split data by condition
    cond1 = (ind & compare_idx == 1);
    cond2 = (ind & compare_idx == 2);
    line_names = {[line_names{1} ' cond1'],[line_names{1} ' cond2']};
    % save binned_data in variable here for easier reference.
    data_a = binned_data.(traces{1})(cond1,plot_bins);
    data_b = binned_data.(traces{1})(cond2,plot_bins);
    % calculate mu and sem for each splitted data, and plot.
    mu = mean(data_a,1,'omitnan');
    sem = std(data_a,1,"omitnan")/sqrt(size(data_a,1));
    h = shadedErrorBar(x_data, mu, sem,...
        'lineProps',{'-','Color',line_colors{1}},...
        'plotAxes',ax);
    lines(1) = h.mainLine;
    mu = mean(data_b,1,'omitnan');
    sem = std(data_b,1,"omitnan")/sqrt(size(data_b,1));
    h = shadedErrorBar(x_data, mu, sem,...
        'lineProps',{'-','Color',line_colors{2}},...
        'plotAxes',ax);
    lines(2) = h.mainLine;
    axis(ax,'tight');
    % do TTest and plot significance on previous plot
    b_pref = ttest2(data_a,data_b,'dim',1,'tail','left') == 1;
    a_pref = ttest2(data_a,data_b,'dim',1,'tail','right') == 1;
    % calculate a location to plot sig stars.
    sig_location = max([mean(data_a,1,'omitnan');mean(data_b,1,'omitnan')], [], 'all') * 1.1;
    if any(a_pref)
        plot(ax, x_data(a_pref), sig_location, '*',...
            'MarkerEdgeColor', line_colors{1}, 'tag','sigMarkers');
    end
    if any(b_pref)
        plot(ax, x_data(b_pref), sig_location, '*',...
            'MarkerEdgeColor', line_colors{2}, 'tag','sigMarkers');
    end
    cur_ylims = ylim;
    ylim(ax, [cur_ylims(1), sig_location*1.1])
end


end
