function [inf_tbl_all, bout] = bout_descriptive(inf_tbl, fr, binsize, ifplot)
    % bout traveled distance distribution: integrate v.
    bout_id = double(inf_tbl.bout);
    vel = double(inf_tbl.velocity); % cm/s
    % - find continouse segment of same boutid
    bout = RunningSeg(bout_id); % maybe not necessary since ubid col
    dt = 1/fr;
    incre_dist = arrayfun(@(s,e) dt * cumsum(vel(s:e),'omitnan'), ...
                        bout.StartIdx, bout.EndIdx, ...
                        'UniformOutput', false); % cm
    totldist = cellfun(@(b) b(end), incre_dist);
    bout = horzcat(bout, table(totldist),table(incre_dist)); 

    % - add a column to inf_tbl with cumsum dist
    [tot_dist, bout_dist, bout_allbin] = deal(zeros(height(inf_tbl),1));
    alledge = 0:binsize:max(bout.totldist);
    pkvel = [];
    for i = 1:height(bout)
        bout_dist(bout.StartIdx(i):bout.EndIdx(i)) = bout.incre_dist{i};
        tot_dist(bout.StartIdx(i):bout.EndIdx(i)) =  bout.totldist(i);
        % bin index for all bin range with fixed bin size
        bout_allbin(bout.StartIdx(i):bout.EndIdx(i)) = discretize(bout.incre_dist{i},alledge);  
        rewarded(i) = any(inf_tbl.info_bout_rewarded(bout.StartIdx(i):bout.EndIdx(i)));
        pkvel(i) = mean(inf_tbl.info_peak_vel(bout.StartIdx(i):bout.EndIdx(i)));
        duration(i) =  mean(inf_tbl.info_bout_length(bout.StartIdx(i):bout.EndIdx(i)))./fr;
    end 
    bout.rewarded(:) = rewarded;
    bout.pkvel(:) = pkvel;
    bout.duration(:) = duration;
    if ifplot
        % bout distance distribution (cm)
        figure(1)
        histogram(totldist,50); xlabel('bout distance (cm)')
        hold on
        xline(mean(totldist),'LineWidth',2,'Color','r','DisplayName','Mean bout dist')
        title(sprintf('median bout dist = %.1f cm, mean bout dist =  %.1f cm', ...
            median(totldist), mean(totldist)))
        % bout length distribution (s)
        bout_t = inf_tbl.info_bout_length/fr; 
        bout_t = unique(bout_t(~isnan(bout_t)));
        figure(2)
        histogram(bout_t, 30); xlabel('bout duration (s)')
        title(sprintf('median bout dur = %.1f s',median(bout_t)))
        % bout distance for reward vs. unrewarded
        figure(3); subplot(1,2,1)
        rew_dist = totldist(rewarded==1);
        histogram(rew_dist,50); xlabel('bout distance (cm)')
        hold on
        xline(mean(rew_dist),'LineWidth',2,'Color','r','DisplayName','Mean bout dist')
        title(sprintf('Rewarded: %d, median bout dist = %.1f cm, mean bout dist =  %.1f cm', ...
            length(rew_dist), median(rew_dist), mean(rew_dist)))
        subplot(1,2,2)
        rew_dist = totldist(rewarded==0);
        histogram(rew_dist,50); xlabel('bout distance (cm)')
        hold on
        xline(mean(rew_dist),'LineWidth',2,'Color','r','DisplayName','Mean bout dist')
        title(sprintf('Not Rewarded: %d, median bout dist = %.1f cm, mean bout dist =  %.1f cm', ...
            length(rew_dist),median(rew_dist), mean(rew_dist)))

        % boutpeak velocity distribution (cm)
        figure(randi(10000))
        histogram(pkvel,50); xlabel('bout peak velocity (cm/s)')
    end
    if ~any(matches({'bout_dist', 'tot_dist'}, fieldnames(inf_tbl)))
        inf_tbl_all = [inf_tbl, table(bout_dist),table(tot_dist), table(bout_allbin)];
    else
        inf_tbl_all = inf_tbl;
    end
end

function segments = RunningSeg(x)
% Returns a table of [StartIdx, EndIdx, Value] for contiguous runs of
% identical (non-NaN) values in a vector x. NaNs are treated as separators.

    if ~isvector(x), error('x must be a vector'); end
    x = x(:).';                      % row vector
    idx = find(~isnan(x));           % positions of valid values
    if isempty(idx)
        segments = table([],[],[], 'VariableNames',{'StartIdx','EndIdx','Value'});
        return
    end

    v = x(idx);                      % valid values
    % A new run starts when value changes OR the original indices are not consecutive
    isNew = [true, (v(2:end) ~= v(1:end-1)) | (idx(2:end) ~= idx(1:end-1)+1)];
    runStarts = find(isNew);
    runEnds   = [runStarts(2:end)-1, numel(v)];

    StartIdx = idx(runStarts).';
    EndIdx   = idx(runEnds).';
    Value    = v(runStarts).';

    segments = table(StartIdx, EndIdx, Value);
end

