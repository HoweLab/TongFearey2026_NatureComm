function T = sourceDataHist(h)
% SOURCEDATAHIST  Bars of a plotted histogram (histogram object).
%   h = histogram(bout.duration, 50);
%   T = sourceDataHist(h)   % call before the figure is closed
% Columns: bin_left, bin_right, count (as drawn: h.BinEdges, h.Values).

T = table(h.BinEdges(1:end-1)', h.BinEdges(2:end)', h.Values', ...
    'VariableNames', {'bin_left','bin_right','count'});
end
