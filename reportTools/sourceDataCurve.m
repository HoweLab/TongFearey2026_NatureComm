function T = sourceDataCurve(dat, traces, varargin)
% SOURCEDATACURVE  Plotted values of LME curves from getPlotDatFromTbl
% (what plotBinned_new / plotTriggered_new draw with 'fromLME', true).
%   T = sourceDataCurve(binned_data, {'dSPN','iSPN'}, 'xName','position_bin')
%   T = sourceDataCurve(vu_binned, {'diff_di'}, 'labels',{'dSPN_minus_iSPN'}, 'pval',true)
%
% One row per x point: x, then <label>_mean and <label>_SEM for each trace
% (the line and shaded region), plus <label>_p_Holm when 'pval' is true
% (only for plots that mark significance, e.g. 'plot_sig', true).
%
% Name-value options:
%   labels - column name for each trace (default: trace names)
%   row    - row of mu/sem/pval to use (default 1)
%   pval   - add Holm-adjusted p per point (default false)
%   xName  - name of the x column (default 'x')
%   x      - x values to write instead of dat.xdata (e.g. dat.bins/31)

ip = inputParser;
ip.addParameter('labels', traces);
ip.addParameter('row', 1);
ip.addParameter('pval', false);
ip.addParameter('xName', 'x');
ip.addParameter('x', []);
ip.parse(varargin{:});
o = ip.Results;

x = o.x;
if isempty(x), x = dat.xdata; end
T = table(x(:), 'VariableNames', {o.xName});
for t = 1:numel(traces)
    d = dat.(traces{t});
    lbl = matlab.lang.makeValidName(o.labels{t});
    T.([lbl '_mean']) = d.mu(o.row, :)';
    T.([lbl '_SEM']) = d.sem(o.row, :)';
    if o.pval
        T.([lbl '_p_Holm']) = d.pval(o.row, :)';
    end
end
end
