function varargout  = plot_subtable(subtbl, y_vars, tlts, varargin)
    ip = inputParser;
    ip.addParameter('xdata', [1:100; linspace(1,100,100)])
    ip.addParameter('fake_info', struct('startLoc',1,'rewardLoc',100))
    ip.addParameter('vel_y', 'mean_velocity')
    ip.addParameter('vel_title', {})
    ip.addParameter('bin_unit','vuBin')
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
 
    figs = cell(numel(subtbl), numel(y_vars)); 
    figv = cell(1, numel(subtbl));
    binned_dat = cell(numel(subtbl), numel(subtbl));
    for t = 1:numel(subtbl) % if multiple tables to plot
        for i = 1:numel(y_vars)
            try
                binned_dat{t,i} = getPlotDatFromTbl(subtbl{t},bin_unit,'y',y_vars{i}, ...
                    'vel_y',vel_y,'xdata',xdata);
                figs{t,i} = figure('Position',[900 100 900 700]);
                plotBinned_new(binned_dat{t,i},fake_info, ...
                    'fromLME',true,'traces',{'dSPN','iSPN'},'plot_vel',false)
                plotBinned_new(binned_dat{t,i},fake_info, ...
                    'fromLME',true,'traces',{'diff_di'}, 'plot_sig',true, ...
                    'line_colors',{[0 0 0]},'plot_vel',false);
                sgtitle(tlts{t,i},'Interpreter', 'none')
            catch
                disp('Not plotted')
            end
        end
        if ~isempty(vel_title) % plot velocity
            figv{t} = figure('Position',[900 100 900 700]);
            plotBinned_new(binned_dat{t,i},fake_info, ...
                'fromLME',true,'traces',{'velocity'}, 'plot_sig',false, ...
                'line_colors',{[0 0 0]},'plot_vel',false);
            sgtitle(vel_title{t},'Interpreter', 'none')
        end
    end
    varargout{1} = figs;
    if ~isempty(vel_title)
        varargout{2} = figv;
    end
    if nargout == 3
        varargout{3} = binned_dat;
    end
end
