function plot_crosscor(avgout, sessout)
    
    behav = fieldnames(avgout);
    x = avgout.(behav{1}).xcorr_lags;
    [~,lg_1] = max(abs(avgout.(behav{1}).xcorr_c)); % lag at max cor
    r_lg = avgout.(behav{1}).xcorr_c(lg_1); % max cor
    [~,lg_2] = max(abs(avgout.(behav{2}).xcorr_c)); % vel
    r_lg_2 = avgout.(behav{2}).xcorr_c(lg_2); % max cor
    max_rg = max(r_lg, r_lg_2);
    min_rg = min([avgout.(behav{1}).xcorr_c(:); avgout.(behav{2}).xcorr_c(:)]);
    figure; 
    
    subplot(1,2,1)
    xlim([min(x), max(x)]); ylim([min_rg, max_rg])
    hold on
    xline(0, 'LineStyle','--','LineWidth',1)
    title([behav{1} ' vs. DI max lag ' num2str(x(lg_1)) ', R=' num2str(r_lg)], 'Interpreter','none')
    
    plot(x, avgout.(behav{1}).xcorr_c', 'LineWidth', 2)
    
    % % significance markers
    % plot(x(avgout.(behav{1}).spear_tho_p < 0.001), ...
    %      r_lg + 0.05, '*', 'Color', 'b', 'Tag', 'sigMarkers')
    
    % highlight one point in red
    idx = find(x == x(lg_1), 1);   
    if ~isempty(idx)
        scatter(x(idx), avgout.(behav{1}).xcorr_c(idx), 40, 'r', 'filled')
    end

    subplot(1,2,2)
    xlim([min(x), max(x)]); ylim([min_rg, max_rg])
    hold on
    xline(0, 'LineStyle','--','LineWidth',1)
    title([behav{2} ' vs. DI max lag ' num2str(x(lg_2)) ', R=' num2str(r_lg_2)],  'Interpreter','none')
    
    plot(x, avgout.(behav{2}).xcorr_c', 'LineWidth', 2)
    % 
    % % significance markers
    % plot(x(logical(avgout.(behav{2}).perm_p)), ...
    %      r_lg_2 + 0.05, '*', 'Color', 'b', 'Tag', 'sigMarkers')
    
    % highlight one point in red
    idx = find(x == x(lg_2), 1);   % if lg_2 is exactly one of the x values
    if ~isempty(idx)
        scatter(x(idx), avgout.(behav{2}).xcorr_c(idx), 40, 'r', 'filled')
    end
    


    % Per session
    if ~isempty(sessout{1,1})
        % max lag distribution
        figure; 
        subplot(1,2,1)
        sessout = sessout(cellfun(@(x) ~isempty(x),sessout ));
        boxViolin(cellfun(@(x) x.(behav{1}).maxlg(2), sessout), ...
            cellfun(@(x) x.(behav{2}).maxlg(2), sessout),get(gcf,'Number'))
        subplot(1,2,2)
        boxViolin(cellfun(@(x) x.(behav{1}).maxlg(1), sessout), ...
            cellfun(@(x) x.(behav{2}).maxlg(1), sessout),get(gcf,'Number'))
        sgtitle('max lag and r-val (at max lag) distribution')
        % mean cor at each lag
        
        acc_cor = cell2mat(cellfun(@(x) x.(behav{1}).cor_lags(1,:),...
            sessout, 'UniformOutput', false));
        vel_cor = cell2mat(cellfun(@(x) x.(behav{2}).cor_lags(1,:),...
            sessout, 'UniformOutput', false));
        figure; xlim([min(x), max(x)])
        plot_error(gca, x, mean(acc_cor), std(acc_cor)./sqrt(size(acc_cor,1)),'Color', [1,0,0]); 
        hold on;  plot_error(gca, x, mean(vel_cor), std(vel_cor)./sqrt(size(vel_cor,1)), 'Color', [0,0,1] ); 
        title('Sess-average correlation at each lag'); 
        % figure; plot(x, acc_cor); title('acc vs DI cor per session')
        % figure; plot(x, vel_cor); title('vel vs DI cor per session')
    end
end
