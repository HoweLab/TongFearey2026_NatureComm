function plot_SigOnTrace(dat, sig, field2plot, sigfield, bins, titleN, fpath, fr)
% dat is the binned_dat that has all the figures overlaid

    allfields = fieldnames(dat);
    x = bins / fr;

    % one color per field
    cmap = lines(numel(allfields));

    for i = 1:numel(field2plot)

        figure(randi(1000000));
        ax = gca;
        hold(ax, 'on');

        hLines = gobjects(numel(allfields), 1);

        for j = 1:numel(allfields)   % put low, high or data of same type tgt

            mu  = dat.(allfields{j}).(field2plot{i}).mu;
            sem = dat.(allfields{j}).(field2plot{i}).sem;

            % plot mean +/- sem with assigned color
            hLines(j) = plot_error_colored(ax, x, mu, sem, cmap(j,:));

        end

        % significance markers
        pval = sig.(sigfield{i}).pval;
        if ~isempty(x(pval < 0.05))
            plot(ax, x(pval < 0.05), zeros(sum(pval < 0.05),1), '*', ...
                'Color', 'r', 'Tag', 'sigMarkers');
        end

        xlim(ax, [min(x), max(x)]);
        xline(ax, 0, '--');
        yline(ax, 0, '--');

        % legend is the clearest way to show color-field mapping
        legend(ax, hLines, allfields, 'Interpreter', 'none', 'Location', 'best');

        % title text listing fields in the same order as the plotted colors
        colorMapText = strjoin(allfields', ', ');
        title(ax, sprintf('%s%s | Color order: %s', titleN, field2plot{i}, colorMapText), ...
            'Interpreter', 'none');

        save_img(gcf, fpath, [titleN, field2plot{i}]);
        close(gcf);

    end
end


function h = plot_error_colored(ax, x, mu, sem, c)
% plots mean with shaded SEM in color c

    % make sure row vectors
    x   = x(:)';
    mu  = mu(:)';
    sem = sem(:)';

    % shaded error region
    fill(ax, [x fliplr(x)], [mu-sem fliplr(mu+sem)], c, ...
        'FaceAlpha', 0.20, 'EdgeColor', 'none');

    % mean trace
    h = plot(ax, x, mu, 'Color', c, 'LineWidth', 1.5);
end