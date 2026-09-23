function plotTriggered_new(trig_data,varargin)
% I will write a more generalized version of plotTriggeredData, somewhat
% like plotBinned. So that the programs allows more free modification by
% vriable inputs and take in minimal amount of fixed input.

ip = inputParser;
% ip.addParameter('isGUI',false);
ip.addParameter('fromLME',false)
ip.addParameter('trig_events',{'onsets','peaks','offsets'})
ip.addParameter('traces',{'avgDSPN','avgISPN'});
ip.addParameter('fig',gcf); % plot on current figure unless specified
ip.addParameter('plot_idx',[]); % any specified index to plot, assumed to be boolean index
% program will decide if plot_idx are trials or movebout depending on the
% value for isGUI.
ip.addParameter('x_label','sec');
ip.addParameter('y_label','dFF');
ip.addParameter('line_colors', {});
ip.addParameter('plot_title','');
ip.addParameter('plot_vel', true);
ip.addParameter('plot_sig',false);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% if not given a plot index then plot all trig_events
if ~fromLME
    if isempty(plot_idx)
        plot_idx = true(size(trig_data.(trig_events{1}).(traces{1}).events));
    end
    if ~islogical(plot_idx)
        temp = plot_idx;
        plot_idx = false(size(trig_data.(trig_events{1}).(traces{1}).events));
        plot_idx(temp) = true;
        clear temp
    end
    num_events = sum(plot_idx);
    disp(['Plotting ' num2str(num_events) ' movebouts out of '...
        num2str(numel(plot_idx)) ' total movebouts'])
end



% Plot settings like linecolors
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
lines = gobjects(1,numel(traces));
line_names = traces;
% try to loop through things so I can plot for varying number of axes.
tile = tiledlayout('horizontal','TileSpacing','tight');
all_axes = gobjects(numel(trig_events),1);
ax_names = trig_events;
y_lims = [0,0];
x_lims = zeros(numel(trig_events),2);
yvel_lims = [0,0];
for a = 1:numel(trig_events)
    all_axes(a) = nexttile(tile);
    hold(all_axes(a),'on')
    cur_dat = trig_data.(trig_events{a});
    % get mu and sem of current data for plotting.
    if fromLME
        xdata = cur_dat.xdata;
    else
        xdata = cur_dat.(traces{1}).windowIdx ./ 31';
    end
    x_lims(a,:) = [min(xdata),max(xdata)];
    for t = 1:numel(traces)
        if fromLME
            mu = cur_dat.(traces{t}).mu;
            sem = cur_dat.(traces{t}).sem;
        else
            activity = cur_dat.(traces{t}).activity(:,:,plot_idx);
            mu = mean(activity,3,'omitnan')';
            sem = std(activity,0,3,'omitnan')'/sqrt(num_events);
        end
        h = shadedErrorBar(xdata,mu,sem,'lineProps', ...
            {'-','Color',line_colors{rem(t-1,num_c)+1}}, ...
            'plotAxes',all_axes(a));
        lines(t) = h.mainLine;
        % add sig points of the same color as plot line if input specified
        % plot sig
        y_lims(1) = min(y_lims(1),min(mu-sem));
        y_lims(2) = max(y_lims(2),max(mu+sem));
        if plot_sig && isfield(cur_dat.(traces{t}),'pval')
            pval = cur_dat.(traces{t}).pval;
            if any(pval< 0.05)
                plot(xdata(pval < 0.05),0,'*', ...
                    'Color',line_colors{rem(t-1,num_c)+1}, ...
                    'tag','sigMarkers');
                y_lims(1) = min(y_lims(1),min(mu-sem)) * 1.05;
                y_lims(2) = max(y_lims(2),max(mu+sem)) * 1.05;
            end
        end
    end
    % add xline and yline
    xline(all_axes(a),0,'k--')
    yline(all_axes(a),0,'k--')
    
    % after plotting all traces, plot velocity on the right axes
    if plot_vel
        if fromLME
            mu_vel = cur_dat.velocity.mu;
            sem_vel = cur_dat.velocity.sem;
        else
            velocity = cur_dat.velocity.activity(:,:,plot_idx);
            mu_vel = mean(velocity,3,'omitnan')';
            sem_vel = std(velocity,0,3,'omitnan')'/sqrt(num_events);
        end
        yyaxis(all_axes(a),'right')
        h = shadedErrorBar(xdata,mu_vel,sem_vel,'lineProps','-b', ...
            'plotAxes',all_axes(a));
        vel_line = h.mainLine;
        yvel_lims(1) = min(yvel_lims(1),min(mu_vel-sem_vel));
        yvel_lims(2) = max(yvel_lims(2),max(mu_vel+sem_vel));
    end
end

if plot_vel
    lines(end+1) = vel_line;
    line_names{end+1} = 'velocity';
    line_names = replace(line_names,'_', ' ');

    % ------------------------------------------------------------
    % Force left and right y-axes to share the same y = 0 line
    % ------------------------------------------------------------

    % make sure both global limits include zero
    y_lims(1)    = min(y_lims(1), 0);
    y_lims(2)    = max(y_lims(2), 0);
    yvel_lims(1) = min(yvel_lims(1), 0);
    yvel_lims(2) = max(yvel_lims(2), 0);

    % add small padding
    left_pad = 0.05 * range(y_lims);
    if left_pad == 0
        left_pad = 1;
    end
    y_lims = [y_lims(1) - left_pad, y_lims(2) + left_pad];

    vel_pad = 0.05 * range(yvel_lims);
    if vel_pad == 0
        vel_pad = 1;
    end
    yvel_lims = [yvel_lims(1) - vel_pad, yvel_lims(2) + vel_pad];

    % fraction of vertical axis where zero appears on the left axis
    zero_frac = (0 - y_lims(1)) / range(y_lims);

    % protect against edge cases
    zero_frac = max(min(zero_frac, 1 - eps), eps);

    % expand velocity limits so zero appears at the same vertical fraction
    vel_range_needed_lower = abs(yvel_lims(1)) / zero_frac;
    vel_range_needed_upper = abs(yvel_lims(2)) / (1 - zero_frac);

    vel_range_final = max(vel_range_needed_lower, vel_range_needed_upper);

    yvel_lims_aligned = [
        -zero_frac * vel_range_final, ...
         (1 - zero_frac) * vel_range_final
    ];

    % apply limits to every subplot
    for a = 1:numel(all_axes)

        yyaxis(all_axes(a),'left')
        all_axes(a).YLim = y_lims;
        all_axes(a).XLim = x_lims(a,:);

        if a ~= 1
            all_axes(a).YAxis(1).Visible = 'off';
        end

        yyaxis(all_axes(a),'right')
        all_axes(a).YLim = yvel_lims_aligned;

        if a ~= numel(trig_events)
            all_axes(a).YAxis(2).Visible = 'off';
        else
            all_axes(a).YAxis(2).Color = 'b';
        end

        yyaxis(all_axes(a),'left')
    end

else
    % no velocity axis: only set left y-limits
    y_lims(1) = min(y_lims(1), 0);
    y_lims(2) = max(y_lims(2), 0);

    left_pad = 0.05 * range(y_lims);
    if left_pad == 0
        left_pad = 1;
    end
    y_lims = [y_lims(1) - left_pad, y_lims(2) + left_pad];

    for a = 1:numel(all_axes)
        % no yyaxis call here: it would add an empty right-hand y axis
        all_axes(a).YLim = y_lims;
        all_axes(a).XLim = x_lims(a,:);

        if a ~= 1
            all_axes(a).YAxis(1).Visible = 'off';
        end
    end
end

% move significance markers to top of the final left axis limits
mks = findobj(fig,'Tag','sigMarkers');
arrayfun(@(x) set(x,'YData', y_lims(2) * 0.95), mks);



% if plot_vel
%     lines(end+1) = vel_line;
%     line_names{end+1} = 'velolcity';
%     line_names = replace(line_names,'_', ' ');
%     for a = 1:numel(all_axes)
%         % set ax limit for the velocity side
%         yyaxis(all_axes(a),'right')
%         all_axes(a).YLim = yvel_lims;
%         if a ~= numel(trig_events)
%             all_axes(a).YAxis(2).Visible = 'off';
%         else
%             all_axes(a).YAxis(2).Color = 'b';
%         end
%         yyaxis(all_axes(a),'left')
%     end
% end
% % change plot limits
% % get the range of y_limits to set, then set ylim for all existing axis
% for a = 1:numel(all_axes)
%     % set ax limit for the signal side
%     all_axes(a).YLim = y_lims;
%     all_axes(a).XLim = x_lims(a,:);
%     if a ~= 1
%         all_axes(a).YAxis(1).Visible = 'off';
%     end
% end
mks = findobj(fig,'Tag','sigMarkers');
arrayfun(@(x) set(x,'YData', y_lims(2) * 0.95),mks)


% plot settings, prepare axes
if numel(traces) == 2 && ~fromLME
    plot_diff = true;
else
    plot_diff = false;
end
if plot_diff
    fig = figure(); % do difference plot on a different figure.
    % For this segement I need to first add a graphic object place holder
    % at the end of my array of lines. 
    diff_dat = struct();
    for a = 1:numel(all_axes)
        cur_dat = trig_data.(trig_events{a});
        data_a = cur_dat.(traces{1}).activity(:,:,plot_idx);
        data_b = cur_dat.(traces{2}).activity(:,:,plot_idx);
        % 
        % b_pref = ttest2(data_a,data_b,'dim',3,'tail','left') == 1;
        % a_pref = ttest2(data_a,data_b,'dim',3,'tail','right') == 1;

        % plot the black DvI line
        DvI = data_a - data_b;
        diff_dat.(trig_events{a}).diff.mu = mean(DvI,3,'omitnan')';
        diff_dat.(trig_events{a}).diff.sem = std(DvI,0,3,'omitnan')'/sqrt(num_events);
        [~,pval,~] = ttest2(data_a,data_b,'dim',3,'tail','both');
        diff_dat.(trig_events{a}).diff.pval = pval;
           
        % will make a recursive call, compute other neccesary input to do it.
        diff_dat.(trig_events{a}).xdata = cur_dat.(traces{1}).windowIdx ./ 31';
        velocity = cur_dat.velocity.activity(:,:,plot_idx);
        diff_dat.(trig_events{a}).velocity.mu = mean(velocity,3,'omitnan')';
        diff_dat.(trig_events{a}).velocity.sem = std(velocity,0,3,'omitnan')'/sqrt(num_events);
    end
    plotTriggered_new(diff_dat,'traces',{'diff'},'fromLME',true,'fig',fig, ...
        'plot_sig',true,'line_colors',{[0 0 0]});
end

legend(lines,line_names,'Location','northeastoutside')
sgtitle(tile, plot_title);

end

