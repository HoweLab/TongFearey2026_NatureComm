function T = sourceDataBouts(tbl, bout, vars)
% SOURCEDATABOUTS  One row per bout from bout_descriptive, with its
% mouse / session / field / cell-type half and the requested bout values.
%   [~, bout] = bout_descriptive(tblE, 31.25, 1.5, 0);
%   T = sourceDataBouts(tblE, bout, {'duration'})
%
% tbl must be the same table passed to bout_descriptive (bout.StartIdx
% indexes its rows). vars are columns of bout, e.g. 'duration', 'pkvel'.

i0 = bout.StartIdx;
T = table(string(tbl.mouse(i0)), string(tbl.sess(i0)), string(tbl.field(i0)), ...
    string(tbl.celltype(i0)), string(tbl.bout(i0)), ...
    'VariableNames', {'mouse','sess','field','celltype_half','bout'});
for v = cellstr(vars)
    T.(v{1}) = bout.(v{1})(:);
end
end
