function row = statsRow(varargin)
% STATSROW  One row of the common stats-report schema shared by all figure
% scripts. Every summary reporter (reportGroupCompare, reportCrossCorr, ...) builds
% its output with this function so rows from any figure can be vertically
% concatenated and written to Excel as a block with writeSourceData.
%
%   row = statsRow('Figure','Fig1','Panel','E','Measure','dF/F', ...
%                  'Group','dSPN','Estimate',8.1,'SEM',0.4,'pValue',1e-5, ...
%                  'nStruct', n_trial)
%
% Columns not given are left empty (text) or NaN (numeric).
% 'nStruct' (optional) is a sample-size struct from countSampleSize /
% countCorrSampleSize; it fills the N_* columns unless they are given
% explicitly. N_neurons is taken from n.nNeurons if present, otherwise from
% n.nDSPN / n.nISPN depending on whether Group names dSPN or iSPN.
% 'Sig' is filled from pValue (*** <.001, ** <.01, * <.05, n.s.) if not given.

cols = { ... name,            type
    'Figure',            't'
    'Panel',             't'
    'Measure',           't'
    'Group',             't'
    'Test',              't'
    'Window',            't'
    'X',                 'n'
    'Estimate',          'n'
    'EstimateType',      't'
    'SEM',               'n'
    'StatName',          't'
    'Stat',              'n'
    'pValue',            'n'
    'pType',             't'
    'Sig',               't'
    'N',                 'n'
    'N_unit',            't'
    'N_mice',            'n'
    'N_sessions',        'n'
    'N_fields',          'n'
    'N_neurons',         'n'
    'N_trials_or_bouts', 'n'
    'Note',              't'};

ip = inputParser;
for c = 1:size(cols,1)
    if cols{c,2} == 't'
        ip.addParameter(cols{c,1}, "");
    else
        ip.addParameter(cols{c,1}, NaN);
    end
end
ip.addParameter('nStruct', []);
ip.parse(varargin{:});
r = ip.Results;
isDefault = @(name) ismember(name, ip.UsingDefaults);

% fill sample sizes from the n struct unless given explicitly
n = r.nStruct;
if ~isempty(n)
    fillMap = {'N','nUnits'; 'N_unit','unit'; 'N_mice','nMice'; ...
        'N_sessions','nSessions'; 'N_fields','nFields'; ...
        'N_trials_or_bouts','nTrialsOrBouts'};
    for f = 1:size(fillMap,1)
        if isDefault(fillMap{f,1}) && isfield(n, fillMap{f,2})
            r.(fillMap{f,1}) = n.(fillMap{f,2});
        end
    end
    if isDefault('N_neurons')
        grp = string(r.Group);
        if isfield(n,'nNeurons')
            r.N_neurons = n.nNeurons;
        elseif isfield(n,'nDSPN') && contains(grp,'dSPN','IgnoreCase',true) && ~contains(grp,'iSPN','IgnoreCase',true)
            r.N_neurons = n.nDSPN;
        elseif isfield(n,'nISPN') && contains(grp,'iSPN','IgnoreCase',true) && ~contains(grp,'dSPN','IgnoreCase',true)
            r.N_neurons = n.nISPN;
        end
    end
end

if isDefault('Sig') && ~isnan(r.pValue)
    r.Sig = sigStars(r.pValue);
end

row = table();
for c = 1:size(cols,1)
    v = r.(cols{c,1});
    if cols{c,2} == 't'
        row.(cols{c,1}) = string(v);
    else
        row.(cols{c,1}) = double(v);
    end
end
end

function s = sigStars(p)
if p < 0.001
    s = "***";
elseif p < 0.01
    s = "**";
elseif p < 0.05
    s = "*";
else
    s = "n.s.";
end
end
