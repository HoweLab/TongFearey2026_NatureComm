%% Main function
function [binned_data, figs] = get_DIdiff_glm_SX(tbl_exp0, tbl_exp1, varargin)

    % use glme to get the experimental level difference of DI-diff
    
    ip = inputParser;
    %  change naming for concat: size = [n_vars_to_change x 2]
    ip.addParameter('matchingVar',{'trial','bout'; 'vuBin','progressBin'}) 
    ip.addParameter('y','mean_avgFc3') % response var
    ip.addParameter('field','vuBin')   % for virmen bins
    ip.addParameter('xdata', [1:100; linspace(1,100,100)])
    ip.addParameter('fake_info', struct('startLoc',1,'rewardLoc',100))
    ip.addParameter('plotornot', 1)
    ip.addParameter('vel_y', 'mean_velocity')
    ip.addParameter('randomEffects','(1|mouse) + (1|mouse:sess) + (1|mouse:sess:field)')
    ip.addParameter('compare','DI_diff_exp') % estimate diff of DI diff
    ip.addParameter('plot_v',0)
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    
    var0 = tbl_exp0.Properties.VariableNames; % track
    var1 = tbl_exp1.Properties.VariableNames; % gui
    
    if any(~matches(var0, var1))   % if var name not match exactly
        warning('table input has diff vars, matching variables')
        for v = 1:size(matchingVar,1) % change vars in tbl_exp1 to match exp0
            if contains(matchingVar(v,1), var0)
                tbl_exp1.(matchingVar{v,1}) = tbl_exp1.(matchingVar{v,2});
            end  
        end
        vars_both = intersect(var0, tbl_exp1.Properties.VariableNames);
        tbl_exp00 = tbl_exp0(:, vars_both);
        tbl_exp11 = tbl_exp1(:, vars_both);
    else
        tbl_exp00 = tbl_exp0;
        tbl_exp11 = tbl_exp1;
    end
    % Binary variable for experimental condition
    tbl_exp00.exp(:) = 0;
    tbl_exp11.exp(:) = 1; 
    tbl_in = vertcat(tbl_exp00, tbl_exp11);
    % Binary variable for two cell types
    % tbl_in.isD(tbl_in.celltype == 'dSPN') = 1;
    % tbl_in.isI(tbl_in.celltype == 'iSPN') = 1;
    tbl_in.exp0(tbl_in.exp == 0) = 1;
    tbl_in.exp1(tbl_in.exp == 1) = 1;

    binned_data = getPlotDatFromTbl(tbl_in,field,'y',y, ...
        'vel_y',vel_y,'xdata',xdata ,'compare', compare);
    figs = struct;
    % Plot
    if plotornot
        figs.DI_exp0 = figure('Position',[900 100 900 700]);
        plotBinned_new(binned_data,fake_info, ...
            'fromLME',true,'traces',{'dSPN','iSPN'}, ...
            'plot_vel',false)
        plotBinned_new(binned_data,fake_info, ...
            'fromLME',true,'traces',{'diff_di'}, ...
            'plot_sig',true,'line_colors',{[0 0 0]});
        
        figs.figExp01 = figure('Position',[900 100 900 700]);
        plotBinned_new(binned_data,fake_info, ...
            'fromLME',true,'traces',{'diff_DI_exp01'}, ...
            'plot_sig',true)
        if plot_v
            figs.v_exp0 = figure('Position',[900 100 900 700]);
            plot(binned_data.velocity.mu(1:3,:)','LineWidth',2)
            legend('v_exp0','v_exp1','v_diff')
        end
    end
end