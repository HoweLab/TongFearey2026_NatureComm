function T = sourceDataPerMouse(bin_out_struct, x, varargin)
% SOURCEDATAPERMOUSE  Per-mouse average curves as plotted by
% addSingleMouseLinePlot (one line per mouse = mean over that mouse's trials).
%   [~, acc_mouse] = addSingleMouseLinePlot(gca, vuBin_tbl, 'vuBin', 'mean_acceleration', ...);
%   T = sourceDataPerMouse(acc_mouse, mouse_bins, 'xName','position_bin')
%
% bin_out_struct.(mouse) is [nTrial x nBin] (second output of
% addSingleMouseLinePlot); x are the bins it used (unique non-NaN bins).
% Output: x column, then one column per mouse.

ip = inputParser;
ip.addParameter('xName', 'x');
ip.parse(varargin{:});

T = table(x(:), 'VariableNames', {ip.Results.xName});
for m = fieldnames(bin_out_struct)'
    T.(m{1}) = mean(bin_out_struct.(m{1}), 1, 'omitnan')';
end
end
