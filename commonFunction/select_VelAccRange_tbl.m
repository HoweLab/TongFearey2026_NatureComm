function out = select_VelAccRange_tbl(tbls, varargin)
% Select comparable bouts between 1w vs. gui 
% grid search for threshoding vel/acc parameters for onset/offset events   
% to find minimal diff score between two exps 

    ip = inputParser;
    ip.addParameter('grp_vars', repmat({'mouse','field','sess','bout'},numel(tbls),1))
    ip.addParameter('bin_unit', repmat({'progressBin'},numel(tbls),1))
    ip.addParameter('BinOrTrig', 'Bin') 
    % Trig related values
    ip.addParameter('win_events', {'onsets','offsets'})
    ip.addParameter('wind', {-31:3*31, 3*-31:31})
    ip.addParameter('min_bout_length', 5 * 31)  
    % location bin: to the max accel at the begining and deccl at the end
    ip.addParameter('loc_bin_idx', {1:20, 81:100}) % in terms of bin_unit
    %ip.addParameter('pk_windpw', {1:20, 1:20})
    ip.addParameter('numRange', 3) % tertile
    ip.addParameter('dataType', 'DI') % D and I data
    ip.addParameter('y_val', {'avgFc3'})
    ip.addParameter('behav_val', {'velocity','acceleration'})
    ip.addParameter('velfd2seg','velocity')
    ip.addParameter('acclfd2seg', 'acceleration')
    % How far the range label (acc_peak, v_mean, ...) is written back into
    % bin_tbl, i.e. which rows a caller filtering on <stat> == i gets back:
    %   'window' old behavior: label kinematic range only on the wind 
    %   'event'  every row of the event the window belongs to, an event being
    %           a run of the win_events column that does not restart. Lets the
    %           kinematics be matched on e.g. 0:2 s while the neural data is
    %           read over -1:3 s
    ip.addParameter('label_scope', 'window')
    ip.addParameter('ifplot', 0)
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end

    % get bout/trial-wise behavioral stats
    if ~iscell(tbls) & isscalar(tbls)
        tbls = {tbls};
    end
    if ~iscell(bin_unit) | ischar(bin_unit)
        bin_unit = {bin_unit};
    end

    if strcmp(BinOrTrig, 'Bin')
        bin_idx = loc_bin_idx;
    elseif strcmp(BinOrTrig, 'Trig')
        bin_idx = wind;
    end
    % initiate
    [win_v, win_a, idxes, bin_tbls] = deal(cell(numel(tbls),numel(bin_idx)));
    [expID] = cell(1, numel(tbls));
    stats_out = cell(1, numel(bin_idx));
    for w = 1:numel(bin_idx) % pos index to be assesse
        curbin = bin_idx{w};
        for t = 1:numel(tbls)
            if strcmp(BinOrTrig, 'Bin')
                cur_vars = grp_vars(t,:);
                if strcmp(dataType, 'DI') % DI data
                    bin_tbl = groupsummary(tbls{t}, [cur_vars(:)' bin_unit{t} {'isD','isI','celltype','numFc3'}], ...
                                'mean', [y_val(:)',behav_val(:)']);  
                    bin_tbl = bin_tbl(~isnan(bin_tbl.(bin_unit{t})) ,:); 
                    % have all isD at the top and all isI at the bottom
                    bin_tbl = [bin_tbl(bin_tbl.isD == 1,:); bin_tbl(bin_tbl.isI == 1,:)];
                else % Dopamine data; potentially y_val is a cell with >1 filbers
                    bin_tbl = groupsummary(tbls{t}, [cur_vars(:)' bin_unit{t}], ...
                                'mean', [y_val(:)',behav_val(:)']);  
                end
                if numRange > 1
                    win_idx = strfind(bin_tbl.(bin_unit{t})',curbin) + (0:length(curbin)-1)';
                    win_v{t,w} =  bin_tbl.mean_velocity(win_idx); % [window size x numTrials]
                    win_a{t,w} =  bin_tbl.mean_acceleration(win_idx);
                end
                bin_tbls{t,w} = bin_tbl; % save for output
            elseif strcmp(BinOrTrig, 'Trig')
                curtbl = tbls{t};
                if contains(curtbl.Properties.VariableNames, 'info_bout_length')
                    curtbl = curtbl(curtbl.info_bout_length > min_bout_length,:);
                end
                if strcmp(dataType, 'DI') % so it works for DA data
                    curtbl = [curtbl(curtbl.isD == 1,:); curtbl(curtbl.isI == 1,:)];
                end
                if numRange > 1 % by pass all the middle thing and output bin window
                    win_idx = strfind(curtbl.(win_events{w})',curbin) + (0:length(curbin)-1)';
                    win_v{t,w} =  curtbl.(velfd2seg)(win_idx); % preserve the size
                    win_a{t,w} = curtbl.(acclfd2seg)(win_idx);
                end
                bin_tbls{t,w} = curtbl; 
            end
            if numRange > 1
            idxes{t,w} = win_idx;
            expID{t} = repmat(t, 1, size(win_v{t,w}, 2));
            end
        end
        % stats for each window event
        statsName = {'v_mean', 'a_mean', 'acc_peak', 'dec_peak', 'pos_a_mean', 'neg_a_mean'} ;
        if numRange > 1
            v = cell2mat(win_v(:,w)');
            accl = cell2mat( win_a(:,w)');
            exp = cell2mat(expID)';
            v_mean = mean(v, 1, 'omitnan')';
            a_mean = mean(accl, 1, 'omitnan')';
            acc_peak = max(max(accl,0), [], 1)'; % max accel (a>0, for track beginning)
            dec_peak = abs(min(min(accl,0), [], 1))'; % max decel (a<0, for track term) 
            pos_a_mean = mean(max(accl,0), 1, 'omitnan')'; % mean of acc segs
            neg_a_mean = abs(mean(min(accl,0), 1, 'omitnan'))'; % abs mean of deceleration segs
            
            
            T = table(v_mean, a_mean,acc_peak, dec_peak, ...
                pos_a_mean, neg_a_mean,  exp);
            
            stats_out{w} = define_statsRange(T, numRange,ifplot); % tertile
          
            % put the range idx back to the table
            for t = 1:numel(tbls)   
                curID = stats_out{w}.rangeID(T.exp == t,:);
                % which rows each matched window's label is written to
                scope_id = scope_labels(bin_tbls{t,w}, win_events{w}, label_scope);
                for s = 1:numel(stats_out{w}.statsName) % loop through each stat
                    col = nan(height(bin_tbls{t,w}),1);  
                    if isempty(scope_id)
                        % rep each val windowSize times
                        col(idxes{t,w}(:)) = repelem(curID(:,s), size(idxes{t,w},1));
                    else
                        % same value on every row of the window's event/bout
                        col = spread_to_scope(col, scope_id, idxes{t,w}(1,:), curID(:,s));
                    end
                    bin_tbls{t,w}.(stats_out{w}.statsName{s}) = col;
                end
            end
        else
            % put the range idx back to the table
            for t = 1:numel(tbls)   
                for s = 1:numel(statsName) % loop through each stat  
                    bin_tbls{t,w}.(statsName{s})(:) = 1;
                end
            end

        end
        
    end
    % output
    out.trial_stats = stats_out;
    out.bin_tbl = bin_tbls;
    out.bin_idx = bin_idx;

end


%% helper
function id = scope_labels(tbl, ev, label_scope)
% Row -> event/bout id, or [] when the label stays on the matched window.
    id = [];
    switch lower(label_scope)
        case 'window'
            return
        case 'event'
            % An event starts wherever the event column becomes non-NaN or
            % stops increasing, so two onsets whose windows abut, and the
            % dSPN/iSPN halves of the table, are never merged.
            v = double(tbl.(ev)); ok = ~isnan(v);
            starts = ok & ([true; ~ok(1:end-1)] | [true; v(2:end) <= v(1:end-1)]);
            id = cumsum(starts);
            id(~ok) = NaN;
        otherwise
            error('select_VelAccRange_tbl:badScope', ...
                'label_scope must be ''window'', ''event''.');
    end
end

function col = spread_to_scope(col, scope_id, win_start_rows, vals)
% Give every row sharing an event/bout id with a matched window that
% window's value. Later windows win if two land in the same event.
    key = scope_id(win_start_rows);
    ok = ~isnan(key);
    if ~any(ok), return; end
    lut = nan(max(scope_id(~isnan(scope_id))), 1);
    lut(key(ok)) = vals(ok);
    have = ~isnan(scope_id);
    col(have) = lut(scope_id(have));
end

function out = define_statsRange(T, numRange, ifplot)
% make sliding window of kinematics across trials
% numRange = 4, quartile; numRange= 3, tertile
% 2 cases:
%   Case1: if only one experiment, just divide itself
%   Case2: if two exps, find range that match bt the two 

    X = T{:,1:end-1};
    expv = T{:,end};
    varNames = T.Properties.VariableNames;
    statNames = varNames(1:end-1); % kinematic statics 
    prct = linspace(100/numRange, 100, numRange);
    if numel(unique(expv)) == 1 % Case 1
        p = prctile(X, prct);
        binID = assign_binID(X, p);
        p1 = []; p2 = [];
    elseif numel(unique(expv)) == 2 % when need to get comparable stats for two exps

        % get larger cut-off values bt two distribution through balanced
        % selection
        X1 = X(T.exp == 1,:); 
        X2 = X(T.exp == 2,:); 
        p1 = prctile(X1, prct);
        p2 = prctile(X2, prct);
        p_eq = nan(size(p1));
        for col = 1:size(X1,2)
            p_eq(:, col) = shared_quantile_cutoffs(X1(:, col), X2(:, col), prct);
        end
        
        % beta1 = 0.1; beta2 = 0.0;
        % p = p_eq;
        % % push the lower boundary and lower the upper
        % % boundary of the last tile
        % p(2,:) = beta1 * p_eq(2,:) + (1-beta1) * max(p1(2,:), p2(2,:)); 
        % p(3,:) = beta2 * p_eq(3,:) + (1-beta2) * min(p1(3,:), p2(3,:)); 
        % idx = p1 < p2;
        % p = nan(size(p1));
        % p(idx) = p1(idx);
        % p(idx == 0) = p2(p1 > p2);
        p = p1; % p1 is 1w
        binID1 = assign_binID(X1, p);
        binID2 = assign_binID(X2, p);
        binID = nan(size(X));
        binID(T.exp == 1,:) = binID1;
        binID(T.exp == 2,:) = binID2;
        if ifplot
            figure; histogram(binID1(:,3));hold on; histogram(binID2(:,3))
            figure; histogram(binID1(:,4)); hold on; histogram(binID2(:,4))
        end
    end
    % pack output
    out.rangeVal = p;
    out.rangeVal_p1 = p1;
    out.rangeVal_p2 = p2;
    out.rangeID = binID;
    out.statsName = statNames;
    out.T_stats = T;
end


function [binID] = assign_binID(X, p)
    numRange = size(p,1);
    [m,n] = size(X);
    % avoid looping to assign each value to a range
    binID = squeeze(sum(reshape(X,m,1,n) > reshape(p,1,numRange,n), 2)) + 1;
    %binID(binID > numRange) = nan; % exceed 
end

function p = shared_quantile_cutoffs(X1, X2, prct)
% prct in percent, e.g. [25 50 75]

    X1 = X1(:);
    X2 = X2(:);
    X1 = X1(~isnan(X1));
    X2 = X2(~isnan(X2));

    vals = sort(unique([X1; X2]));
    n1 = numel(X1);
    n2 = numel(X2);

    p = nan(size(prct));

    for k = 1:numel(prct)
        q = prct(k) / 100;

        % empirical CDFs evaluated on shared support
        F1 = arrayfun(@(v) sum(X1 <= v) / n1, vals);
        F2 = arrayfun(@(v) sum(X2 <= v) / n2, vals);

        Fmix = 0.5 * F1 + 0.5 * F2;

        idx = find(Fmix >= q, 1, 'first');
        p(k) = vals(idx);
    end
    % take the middle value for the largest tile
    X1p = prctile(X1, prct);
    X2p = prctile(X2, prct);
    p(end) = mean([X1p(end), X2p(end)]);
end