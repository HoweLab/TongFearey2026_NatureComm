function [sessout,avgout] = crosscor_DI_behav(tbl, varargin)
    % corrcoef(x,y) computes the Pearson correlation at zero lag after mean subtraction.
    % parse input
    ip = inputParser;
    ip.addParameter('BinOrTrig', 'bin')
    ip.addParameter('x_data', [1:100; linspace(1,100,100)])
    % avg bin
    ip.addParameter('grp_vars', {'mouse','field','sess','trial'})
    ip.addParameter('bin_unit', 'vuBin') % onset for trig avg
    % trig avg
    ip.addParameter('whichFc3', 'avgFc3')
    ip.addParameter('whichBehavs', {'acceleration','velocity'})
    ip.addParameter('cor_level',{'mouse','field','sess'}) % sess level cor
    ip.addParameter('onlyAvg', 1)
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    % cross-corr of mean velocity/acceleration with DI diff
    if strcmp(BinOrTrig, 'bin')
        vuBin_tbl = groupsummary(tbl,[grp_vars(:)' bin_unit {'isD','isI','celltype','numFc3'}], ...
            'mean',{whichFc3,'velocity','acceleration'}); % calculate mean for avgFc3 and velocity
        [grps, ~, sessID] = sort_tbl(vuBin_tbl,  cor_level,'ascend'); % tid is mouse ID
        Fc3 = ['mean_',whichFc3];
        behav_name = strcat('mean_', whichBehavs);
    elseif strcmp(BinOrTrig, 'trig')
        vuBin_tbl = tbl;
        [grps, ~, sessID] = sort_tbl(vuBin_tbl,  cor_level, 'ascend'); % tid is mouse ID
        Fc3 = whichFc3;
        behav_name = whichBehavs;
    end
    sessout = cell(height(sessID),1 );
    if ~onlyAvg % no session-level
        % on the session level       
        x = x_data;
        bin = bin_unit;
        for i = 1:height(sessID) 
            fprintf('session %d \n', i)
            curtbl = vuBin_tbl(grps == i,:);
            if ~all(ismember(x(1,:), unique(curtbl.(bin))))
                warning('skip currrent session, not full rank')
                continue
            end
            sessout{i}.sessName = sessID(i,:);
            for b = 1:numel(behav_name)
                binned = groupsummary(curtbl, {bin, 'isD', 'isI'}, 'mean',{behav_name{b}, Fc3});
                idx = ismember(binned.(bin), x(1,:));
                behav = binned.(['mean_',behav_name{b}])((idx & binned.isD == 1));
                DIdiff = binned.(['mean_',Fc3])(idx & binned.isD == 1) - binned.(['mean_',Fc3])(idx & binned.isI == 1) ;
                [c,lags] = xcorr(DIdiff', behav', 'normalized');
                sessout{i}.(behav_name{b}).maxlg = [max(c), lags(c == max(c))];
                sessout{i}.(behav_name{b}).cor_lags = [c;lags];     
                sessout{i}.(behav_name{b}).behavdat = behav; % array
                sessout{i}.(behav_name{b}).DIdat = DIdiff;            
            end
        end
    end

    % All mice
    avgout = struct;
    binned_DI = getPlotDatFromTbl(vuBin_tbl, bin_unit, 'y', Fc3, ...
        'vel_y',[], 'xdata', x_data);

    for b = 1:numel(behav_name)

        binned_beh = getPlotDatFromTbl(vuBin_tbl, bin_unit, 'y', [], ...
            'vel_y', behav_name{b}, 'xdata', x_data);

        x = binned_DI.diff_di.mu(:);          % DI trace
        y = binned_beh.velocity.mu(:);     % behavior trace

        % Ensure same length
        n = min(numel(x), numel(y));
        x = x(1:n);
        y = y(1:n);

        % Compute lagged correlations and p-values
        maxlag = round(numel(x)/2);

        % Keep xcorr-style output too if you want
        [c, lags_xcorr] = xcorr(x, y,maxlag, 'normalized');
        n_overlap = nan(1, numel(lags_xcorr));
        spearp_lags = nan(1, numel(lags_xcorr));
        p_sig = nan(1, numel(lags_xcorr));
        for i = 1:numel(lags_xcorr)
            [xs, ys] = local_align_by_lag(x, y, lags_xcorr(i));
            valid = ~(isnan(xs) | isnan(ys));
            xs = xs(valid);
            ys = ys(valid);
            n_overlap(i) = numel(xs);

            if numel(xs) >= 3
                [s_tho, spearp_lags(i)] = corr(xs, ys, 'type', 'Spearman');
                rnull = nan(1,1000);
                for p = 1:1000
                    s = randi(numel(ys)-1);
                    yperm = circshift(ys, s);
                    rnull(p) = corr(xs, yperm, 'type', 'Spearman');
                end
                p_sig(i) = (sum(abs(rnull) >= abs(s_tho)) + 1) / (numel(rnull) + 1) < 0.05;
            end
        end



        avgout.(behav_name{b}).maxlg =[max(c), lags_xcorr(c == max(c))];
        avgout.(behav_name{b}).spear_tho_p   = spearp_lags;
        avgout.(behav_name{b}).xcorr_c    = c; % correlation at each lag
        avgout.(behav_name{b}).xcorr_lags = lags_xcorr;
        avgout.(behav_name{b}).n_overlap   = n_overlap;
        avgout.(behav_name{b}).perm_p   = p_sig;
        % raw dat
        avgout.(behav_name{b}).behavdat = binned_beh.velocity; % structure
        avgout.(behav_name{b}).DIdat    = binned_DI.diff_di;
    end
end

%% helper


function [xs, ys] = local_align_by_lag(x, y, lag)
% Return overlapping vectors for a given lag.
%
% lag > 0:
%   y is effectively shifted to the right, so compare x(1:end-lag) with y(1+lag:end)
%
% lag < 0:
%   compare x(1-lag:end) with y(1:end+lag)

    n = numel(x);
    if lag >= 0
        xs = x(1:n-lag);
        ys = y(1+lag:n);
    else
        k = -lag;
        xs = x(1+k:n);
        ys = y(1:n-k);
    end
end
