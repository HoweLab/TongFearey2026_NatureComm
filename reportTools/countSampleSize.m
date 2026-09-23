function n = countSampleSize(tbl, unitVar)
% COUNTSAMPLESIZE  Sample sizes of a (selected) table from mkTblForPlot.
%   n = countSampleSize(track_tbl(trial_selection,:), 'trial')
%   n = countSampleSize(track_tbl(trial_selection & bout_selection,:), 'bout')
%
% unitVar is the observation unit ('trial' or 'bout'); rows where it is
% undefined are ignored. Returns a struct usable as the 'n' input of
% statsRow:
%   unit, nUnits (= nTrialsOrBouts), nMice, nSessions, nFields, nDSPN, nISPN,
%   perMouse (table of the same counts per mouse)
% Neuron counts use numFc3 (number of cells of that type in a field),
% counted once per mouse/session/field/celltype.

reqVars = {'mouse','sess','field',unitVar,'celltype','numFc3'};
hasCells = ismember('numFc3', tbl.Properties.VariableNames); % some tables (e.g. 2 world
if ~hasCells                                                 % tuning subtypes) have no cell counts
    reqVars = setdiff(reqVars, {'numFc3'}, 'stable');
end
missingVars = reqVars(~ismember(reqVars, tbl.Properties.VariableNames));
if ~isempty(missingVars)
    error('countSampleSize: missing variable(s): %s', strjoin(missingVars, ', '));
end

tbl = tbl(~isundefined(tbl.(unitVar)), reqVars);
n = counts(tbl, unitVar);
n.unit = string(unitVar) + "s";

% per-mouse breakdown
mice = unique(tbl.mouse);
perMouse = table();
for m = 1:numel(mice)
    c = counts(tbl(tbl.mouse == mice(m),:), unitVar);
    perMouse = [perMouse; struct2table(c)]; %#ok<AGROW>
end
perMouse = addvars(perMouse, string(mice), 'Before', 1, 'NewVariableNames', 'mouse');
n.perMouse = removevars(perMouse, 'nMice');
end

function c = counts(tbl, unitVar)
fld = unique(tbl(:, {'mouse','sess','field'}), 'rows');
c.nMice = numel(unique(tbl.mouse));
c.nSessions = height(unique(tbl(:, {'mouse','sess'}), 'rows'));
c.nFields = height(fld);
c.nUnits = height(unique(tbl(:, {'mouse','sess','field',unitVar}), 'rows'));
if ~ismember('numFc3', tbl.Properties.VariableNames)
    c.nTrialsOrBouts = c.nUnits;
    [c.nDSPN, c.nISPN] = deal(NaN); % no cell counts in this table
    return
end
c.nTrialsOrBouts = c.nUnits;
cells = groupsummary(tbl, {'mouse','sess','field','celltype','numFc3'});
c.nDSPN = sum(cells.numFc3(cells.celltype == 'dSPN'));
c.nISPN = sum(cells.numFc3(cells.celltype == 'iSPN'));
end
