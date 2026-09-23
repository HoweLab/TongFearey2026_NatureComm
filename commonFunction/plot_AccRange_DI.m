function [trig_dat] = plot_AccRange_DI(out,figdir, varargin)

    % set optional input for 1w vs. gui spatial bin 
    ip = inputParser;
    % for getPlotDatFromTbl function
    ip.addParameter('vel_y', 'mean_acceleration')
    ip.addParameter('bin_unit', {'vuBin','progressBin'}) % match table for spBin; match window for trig
    ip.addParameter('y','mean_avgFc3')
    % stats used to select trials/bouts
    ip.addParameter('stats_vars',{'acc_peak','dec_peak'}) % match number of windows
    % for naming/title
    ip.addParameter('exp_name', {'1w','gui'})
    ip.addParameter('binName', 'spBin') % or Trig
    ip.addParameter('wind_name', {'start bins','end bins'})
    % if want significance between levels
    ip.addParameter('ifgetsig', 0) % get significance between levels 
    ip.addParameter('levlplotted', 'all') % plot all levels; otherwise input index [1,3]
    ip.addParameter('ComparisonTp', 'bt_level') % comparison level: bt_exp or bt_level
    % DI or DA data
    ip.addParameter('dataType', 'DI') % D and I data
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    % plot results from function select_VelAccRange_tbl
    [numtbl,numwind]  = size(out.bin_tbl);
    if strcmp(levlplotted, 'all')
        numRange = size(out.trial_stats{1,1}.rangeVal,1);
        levlplotted = 1:numRange;
    end
    
    numRange = length(levlplotted); % like [1,3] if wants only low and high
    if strcmp(dataType, 'DI') & ~ifgetsig % DI but plot seperate range
        [figs, trig_dat, tlts] =  deal(cell(numtbl, numwind, numRange));
    elseif ~strcmp(dataType, 'DI') % DA 
        [figs, trig_dat, tlts,v_bin] =  deal(cell( numwind, numRange));
    elseif ifgetsig & strcmp(dataType, 'DI') & strcmp(ComparisonTp, 'bt_level')   % sig comparison between level
        pairComb = nchoosek(levlplotted, 2);
        [trig_dat] =  deal(cell(numtbl, numwind, size(pairComb, 1)));
    elseif ifgetsig & strcmp(dataType, 'DI') & ~strcmp(ComparisonTp, 'bt_level')  % sig comparison between exp
        [trig_dat] =  deal(cell(numwind, numRange));
    end

    for w = 1:numwind
        x = [out.bin_idx{w}; linspace(min(out.bin_idx{w}), max(out.bin_idx{w}),length(out.bin_idx{w}))];
        fake_info = struct('startLoc',min(x(2,:)),'rewardLoc',max(x(2,:)));
        curbin = bin_unit{w};
        % Method1: plot DI at each level for each experiment
        if strcmp(dataType, 'DI')  & ~ifgetsig 
            for t = 1:numtbl 
                curtbl = out.bin_tbl{t, w};
                if strcmp(binName, 'spBin')
                    curbin = bin_unit{t};
                end
                for i = 1:numRange
                    curlevel = levlplotted(i);
                    tlts{t, w,i} = {strcat(exp_name{t}, '_', dataType,'_', binName,'_', wind_name{w}, '_', mat2str(curlevel))};
                    tblinput = curtbl(curtbl.(stats_vars{w}) == curlevel,:);
                    %numSPN{t, w,i} = count_numFc3(tblinput, curbin, x(1,:));
                    [figs{t,w,i},~,trig_dat{t,w,i}] = plot_subtable({tblinput}, ...
                        {y}, tlts{t, w,i}, 'bin_unit',curbin,'xdata', x, 'fake_info', fake_info, ...
                        'vel_y',vel_y); % get but not plot vel
                end                
            end
            match2figs_wSubplots(reshape(squeeze(figs(:,w,:)), 1, []))
            for t = 1:numtbl
                save_img(squeeze(figs(t,w,:))',figdir, squeeze(tlts(t,w,:))')
            end
            close all  
        % Method 2: get significance for each experiment at each level and across levels 
        elseif strcmp(dataType, 'DI')  & ifgetsig
            if strcmp(ComparisonTp, 'bt_level') % compare kinematic level under same experiment
                for i = 1:size(pairComb, 1)
                    for t = 1:numtbl 
                        if strcmp(binName, 'spBin'),curbin = bin_unit{t};end
                        curtbl = out.bin_tbl{t, w};
                        curPair = sort(pairComb(i, :), 2, 'descend');
                        % get DI difference ComparisonTp
                        [trig_dat{t,w,i}, ~] = get_DIdiff_glm_SX( ...
                            curtbl(curtbl.(stats_vars{w}) == curPair(1),:),... 
                            curtbl(curtbl.(stats_vars{w}) == curPair(2),:), ...
                            'y',y,'fake_info', fake_info, 'field',curbin,'xdata', x, 'plotornot', 0, 'vel_y', vel_y);    
                    end
                end
            else % same kinematics level but across experiment
                for i = 1:numRange
                    curlevel = levlplotted(i);
                    if numtbl == 2 
                        tbl1 = out.bin_tbl{1, w}; % track 
                        tbl2 = out.bin_tbl{2, w}; % GUI 
                        [trig_dat{w,i}, ~] = get_DIdiff_glm_SX( ...
                            tbl1(tbl1.(stats_vars{w}) ==curlevel,:),... 
                            tbl2(tbl2.(stats_vars{w}) == curlevel,:), ...
                            'y',y,'fake_info', fake_info, 'field', bin_unit{w}, 'xdata', x, 'plotornot', 0, 'vel_y', vel_y);  
                    else
                        disp('More than two tables comparison hasnt been written; come back later')
                        exit
                    end
                end
            end
              
        % Method 3 for DA data that compare things differently    
        elseif ~strcmp(dataType, 'DI')
            curtbl = out.bin_tbl(:, w);
            curtbl = vertcat(curtbl{:}); % merge the same window
            curbin = bin_unit{w};
            for i = 1:numRange
                tlts{w,i} = {strcat( dataType,'_', binName,'_', wind_name{w}, '_', mat2str(i))};
                tblinput = curtbl(curtbl.(stats_vars{w}) == i,:);
                [figs{w,i},~,trig_dat{w,i},v_bin{w,i}] = plot_grid_table({tblinput}, ...
                    {y}, tlts{w,i}, 'bin_unit',curbin,'xdata', x, 'fake_info', fake_info, ...
                    'vel_y',vel_y); % get but not plot vel                
            end
            match2figs_wSubplots(figs(w,:))
            save_img(figs(w,:)',figdir, tlts(w,:)')
        end
        
    end
    


    % plot accl outside of the loop above
    if ~ifgetsig % only plot when not getting significance
        f_acc = cell(numtbl, numwind);
        for t = 1:numtbl
            for w = 1:numwind
                accl_name = [exp_name{t} '_' binName '_accl_' wind_name{w}];
                f_acc{t,w}{1} = figure('Name',['Accl_' wind_name{w}]); 
                hold on; colors = lines(numRange);
                for i = 1:numRange
                    if strcmp(dataType, 'DI')
                        plot_error(gca(f_acc{t,w}{1}), [out.bin_idx{w}], trig_dat{t,w,i}{1}.velocity.mu, ...
                            trig_dat{t,w,i}{1}.velocity.sem,'Color',colors(i,:));
                    else
                        plot_error(gca(f_acc{t,w}{1}), [out.bin_idx{w}], v_bin{w,i}{1}.(exp_name{t}).mu, ...
                            v_bin{w,i}{1}.(exp_name{t}).sem,'Color',colors(i,:));
                    end   
                end
                sgtitle(accl_name, 'Interpreter', 'none')
                save_img(f_acc(t,w),figdir, {accl_name})
            end
        end
        close all
    end
end

