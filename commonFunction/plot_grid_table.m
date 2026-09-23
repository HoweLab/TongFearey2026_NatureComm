function [figs, figv, binned_dat, v_binned] = plot_grid_table(subtbl, y_vars, tlts, varargin)

    ip = inputParser;
    ip.addParameter('xdata', [1:100; linspace(1,100,100)])
    ip.addParameter('fake_info', struct('startLoc',1,'rewardLoc',100))
    ip.addParameter('vel_y', 'mean_velocity')
    ip.addParameter('vel_title', {})
    ip.addParameter('bin_unit','vuBin')
    ip.addParameter('randEffects', '(1|mouse) + (1|mouse:sessID)')
    ip.addParameter('savefpath', [])
    ip.addParameter('expname', {'GUI','LinearTrack'})
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    %% timewindow to plot

    figs = cell(numel(subtbl), numel(y_vars)); 
    figv = cell(1, numel(subtbl));
    binned_dat = cell(numel(subtbl), numel(y_vars));
    v_binned = cell(1, numel(subtbl));
    for t = 1:numel(subtbl) % if multiple tables to plot
        for i = 1:numel(y_vars)
            subtbl{t}.sessty(subtbl{t}.sessty ==expname{1}) = 'gui';
            subtbl{t}.sessty(subtbl{t}.sessty ==expname{2}) = 'track'; % to fit function
            binned_dat{t,i} = getPlotDatFromTbl(subtbl{t},bin_unit,'y',y_vars{i}, ...
                'vel_y',vel_y,'xdata',xdata, 'randomEffects',randEffects, ...
                'compare','sesstype','sessTypeVarName','sessty');
            figs{t,i} = figure('Position',[900 100 900 700]);
            plotBinned_new(binned_dat{t,i},fake_info, ...
                'fromLME',true,'traces',{'track', 'gui'},'plot_vel',false)
            
            plotBinned_new(binned_dat{t,i},fake_info, ...
                'fromLME',true,'traces',{'diff_tg'}, 'plot_sig',true, ...
                'line_colors',{[0 0 0]},'plot_vel',false);
            sgtitle(tlts{t,i},'Interpreter', 'none')
            if ~isempty(savefpath)
                save_img(figs{t,i}, savefpath, tlts{t,i})
            end
        end
        % get behav data
        v_binned{t} = getPlotDatFromTbl(subtbl{t},bin_unit,'y',vel_y, ...
            'vel_y',[],'xdata',xdata, 'randomEffects',randEffects, ...
            'compare','sesstype','sessTypeVarName','sessty');
        if ~isempty(vel_title) % plot velocity
            figv{t} = figure('Position',[900 100 900 700]);
            plotBinned_new(v_binned{t},fake_info, ...
                'fromLME',true,'traces',{'track', 'gui'},'plot_vel',false)
            plotBinned_new(v_binned{t},fake_info, ...
                'fromLME',true,'traces',{'diff_tg'}, 'plot_sig',true, ...
                'line_colors',{[0 0 0]},'plot_vel',false);
            sgtitle(vel_title{t},'Interpreter', 'none')
            if ~isempty(savefpath)
                save_img(figv{t}, savefpath, vel_title{t})
            end
        end
        fprintf(' - For table %d: Number of mice %d \n', length(unique(subtbl{t}.mouse)))
        [G, mouseID] = findgroups(subtbl{t}.mouse);
        numSessPerMouse = splitapply(@(s) numel(unique(s)), subtbl{t}.sessID, G);
        result = table(mouseID, numSessPerMouse);
        fprintf(' - Number of sessions used %d \n', sum(result.numSessPerMouse));
   
    end




end


