function [bin_out, bin_out_struct] = addSessionLinePlot(curaxs, tbl, binName, var2avg, varargin)
% ADDSESSIONLINEPLOT  One line per session on top of a binned plot, the
% session-wise counterpart of addSingleMouseLinePlot (gray traces in the
% figures: "individual sessions").
%   [~, vel_sess] = addSessionLinePlot(gca, switch_tbl(fam,:), 'vuBin', 'velocity')
%
% Each line is that session's mean over trials, per bin (dSPN rows only, as
% in addSingleMouseLinePlot: the two cell-type halves hold the same behaviour).
% Returns the per-session [nTrial x nBin] matrices, as a cell and as a
% struct keyed by mouse_session (for sourceDataPerMouse).
%
% Name-value options:
%   Color - line color (default gray [0.7 0.7 0.7])
%   sessVars - grouping variables defining a session (default {'mouse','sess'})

ip = inputParser;
ip.addParameter('Color', [0.7 0.7 0.7]);
ip.addParameter('sessVars', {'mouse','sess'});
ip.parse(varargin{:});
o = ip.Results;

grp = groupsummary(tbl, [o.sessVars, {'field','trial',binName,'isD','isI'}], 'mean', var2avg);
grp = grp(grp.isD == 1 & ~isnan(grp.(binName)), :);
allbins = unique(grp.(binName));
[sessIdx, sessID] = findgroups(grp(:, o.sessVars));

if isempty(curaxs)
    curaxs = axes('Parent', figure(randi(1000000)));
end
hold(curaxs, 'on');

bin_out = cell(1, height(sessID));
bin_out_struct = struct;
for s = 1:height(sessID)
    cur = grp(sessIdx == s, :);
    if isempty(cur), continue; end
    bin_out{s} = bin_data_colvec(cur.(['mean_' var2avg]), cur.(binName), allbins);
    name = matlab.lang.makeValidName(strjoin(string(sessID{s,:}), '_'));
    bin_out_struct.(name) = bin_out{s};
    plot(curaxs, allbins, mean(bin_out{s}, 1, 'omitnan'), ...
        'Color', o.Color, 'LineWidth', 1, 'DisplayName', char(name));
end
xlim(curaxs, [min(allbins), max(allbins)])
end
