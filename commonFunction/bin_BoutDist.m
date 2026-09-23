% simple binning bout by traveled distance
function [inf_tbl, numbin, nm] = bin_BoutDist(inf_tbl, bout_tbl, bt_range, varargin)
    % 10/16/25: SX
    % inf_tbl: from mkTblForPlot with all data for each time point
    % bout_tbl: discriptive table containing index for each bout
    % bt_range: distance range traveled within a bout to be assessed

    ip = inputParser;
    ip.FunctionName = 'bin_BoutDist';
    ip.addRequired('inf_tbl', @(x) istable(x));
    ip.addRequired('bout_tbl', @(x) istable(x));
    ip.addRequired('bt_range', @(x) isnumeric(x) && numel(x) == 2);
    
    ip.addParameter('sr', 31.25, @(x) isnumeric(x) && isscalar(x));        % sample rate of the 2p camera
    ip.addParameter('min_boutdist', 150, @(x) isnumeric(x) && isscalar(x)); % cm
    ip.addParameter('max_boutdist', inf, @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('min_boutdur', 5, @(x) isnumeric(x) && isscalar(x));    % sec
    ip.addParameter('max_boutdur', inf, @(x) isnumeric(x) && isscalar(x));
    ip.addParameter('animal2exclude', nan, @(x) iscellstr(x) || isstring(x));
    ip.addParameter('dist2reward',[0,1], @(x) isnumeric(x) && numel(x) == 2);
    ip.addParameter('BinSize', 1, @(x) isnumeric(x) && isscalar(x));        % cm

    ip.parse(inf_tbl, bout_tbl, bt_range, varargin{:});
    for j = string(fieldnames(ip.Results))'
        eval(j + " = ip.Results." + j + ";");
    end

    numBout = height(bout_tbl); % this is double the real bouts num
    if ~(contains('ubid', fieldnames(inf_tbl)) && ...
        max(inf_tbl.ubid) == numBout/2)
        error('Number of bouts not matching bt. inf and bout table')
    end

    % Filter bouts
    ok_bout = false(1,numBout);
    b_idx = false(height(inf_tbl),1);
    if contains({'StartIdx','EndIdx','totldist','incre_dist'},fieldnames(bout_tbl))
        dur = (bout_tbl.EndIdx - bout_tbl.StartIdx)/sr;
        dist = bout_tbl.totldist;
        ok_bout (dur>min_boutdur & dur<max_boutdur &...
             dist>min_boutdist &  dist<max_boutdist) = true;
        for i = find(ok_bout)
            b_idx(bout_tbl.StartIdx(i):bout_tbl.EndIdx(i)) = true;
            b_idx(bout_tbl.StartIdx(i):bout_tbl.EndIdx(i)) =  true;
        end
        if ~isnan(animal2exclude)
            excluded = inf_tbl.info_normalized_dist2rew < dist2reward(1) | ...
            inf_tbl.info_normalized_dist2rew > dist2reward(2) |...
             ~contains(cellstr(inf_tbl.mouse),animal2exclude); % not indexed based on other factors
        else
            excluded = inf_tbl.info_normalized_dist2rew < dist2reward(1) | ...
            inf_tbl.info_normalized_dist2rew > dist2reward(2); % not indexed based on other factors
        end
        b_idx(excluded) = false;
    else
        error('run RunningSeg to get bout table with correct variables')
    end
    
    % Edge with set absolute bin size
    b = bt_range(1); e= bt_range(2);
    edges = b:BinSize:e;
    numbin = length(edges);

    % Bin time series data
    d_csum = inf_tbl.bout_dist;
    subtbl_idx = d_csum>b & d_csum<e & b_idx;
    bins = discretize(inf_tbl(subtbl_idx,:).bout_dist, edges);

    bin_tbl = nan(height(inf_tbl),1);
    bin_tbl(subtbl_idx,1) =  bins;
    nm = strcat('Bout_DistBin_', string(b),'_',string(e));
    inf_tbl.(nm) = bin_tbl;
    
end