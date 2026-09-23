function T = sourceDataPoints(groups, data, varargin)
% SOURCEDATAPOINTS  Every point of a box / violin / paired-line plot.
%   T = sourceDataPoints({'track-sensitive','track-insensitive'}, {rho_ts, rho_ti})
%   T = sourceDataPoints({'dSPN','iSPN'}, {d, i}, 'paired',true, 'valueName','pct_sensitive')
%
% Columns: Group, point (index within the group; the pair/session index when
% 'paired' is true, so the connecting lines can be rebuilt) and the value.
%
% Name-value options:
%   valueName - name of the value column (default 'value')
%   paired    - true when the groups are paired point-by-point (default false)
%   ids       - labels for the points (one per point, same for all groups)

ip = inputParser;
ip.addParameter('valueName', 'value');
ip.addParameter('paired', false);
ip.addParameter('ids', []);
ip.parse(varargin{:});
o = ip.Results;

groups = cellstr(groups);
if o.paired && numel(unique(cellfun(@numel, data))) > 1
    error('sourceDataPoints: paired groups must have the same number of points');
end

T = table();
for g = 1:numel(groups)
    x = data{g}(:);
    idx = (1:numel(x))';
    if isempty(o.ids)
        pointName = idx;
    else
        pointName = string(o.ids(:));
    end
    Tg = table(repmat(string(groups{g}), numel(x), 1), pointName, x, ...
        'VariableNames', {'Group','point',o.valueName});
    T = [T; Tg]; %#ok<AGROW>
end
end
