function out = load2wDataset(varargin)
% LOAD2WDATASET  Two-world table, trial selection and single-cell data with
% the same sessions / mice excluded everywhere (Fig 3, Sup 6, Sup 7).
%   dat = load2wDataset();                       % default paths below
%   dat = load2wDataset('scVars',{'dworldcor','iworldcor'});
%
% Excluded: mice that did not learn, the first two-world session of each
% mouse (from the switch-info sheet) and mismatched sessions.
%
% Output fields:
%   tbl             - two-world table from mkTblForPlot (all rows)
%   trial_selection - logical index into tbl: trials to use
%   fam / nov       - logical index into tbl: Track 1 (familiar, info_trial_novel
%                     == 0) and Track 2 (novel) trials, already selected
%   sc              - single-cell variables, rows already restricted to the
%                     valid dSPNs / iSPNs (drm_w1, dstb, dremap, dspn, ...
%                     plus anything listed in 'scVars'), and iPlot
%   info            - cat_info entries of the sessions actually used; sum a
%                     subtype field over it for cell counts, e.g.
%                     sum([dat.info.rmp_w1_d])
%   exclude_mice, exclude_sess, mismatched - what was dropped
%
% Name-value options (defaults):
%   fpath   - two-world tuning-subtype file
%   scpath  - single-cell (by cell) file
%   Tpath   - novel/familiar switch info sheet
%   exclude_mice - {'TdTG7_07'} (animal did not learn)
%   scVars  - extra variables to read from scpath

ip = inputParser;
ip.addParameter('fpath',  '');
ip.addParameter('scpath', '');
ip.addParameter('Tpath',  '');
ip.addParameter('exclude_mice', {'TdTG7_07'});
ip.addParameter('scVars', {});
ip.parse(varargin{:});
o = ip.Results;

% ---- table and trial selection ----
conditions = {'exp_stage',{'2 world'};'mismatch',0};
out.tbl = mkTblForPlot(o.fpath,conditions);
out.exclude_mice = o.exclude_mice;
out.exclude_sess = findfirstSession(o.Tpath); % first 2-world session of each mouse

sel = true(height(out.tbl),1);
sel(out.tbl.skip_trial | ...
    matches(cellstr(out.tbl.mouse),out.exclude_mice) | ...
    contains(cellstr(out.tbl.fname),out.exclude_sess)) = false;
out.trial_selection = sel;
out.fam = sel & out.tbl.info_trial_novel == 0; % Track 1, familiar world (w1)
out.nov = sel & out.tbl.info_trial_novel == 1; % Track 2, novel world (w2)

% ---- session info: mismatched sessions and the sessions actually used ----
S = load(o.fpath, 'cat_info');
cat_info = S.cat_info;
is2w = contains({cat_info.exp_stage}, '2 world');
mm = cat_info(is2w & logical([cat_info.mismatch]));
out.mismatched = strcat({mm.mouse},'_',{mm.session})';
info_2w = cat_info(is2w & ~logical([cat_info.mismatch]));
info_name = strcat({info_2w.mouse},'_',{info_2w.session},'_',{info_2w.field});
out.info = info_2w(ismember(info_name, unique(cellstr(out.tbl(sel,:).fname))));

% ---- single-cell data, valid neurons only ----
scVars = unique([{'dspn','ispn','dstb','istb','dremap','iremap', ...
    'drm_w1','drm_w2','irm_w1','irm_w2','drecnum','irecnum', ...
    'dFieldWidth','iFieldWidth','iPlot'}, o.scVars(:)'], 'stable');
sc = load(o.scpath, scVars{:});
d_info = [sc.dspn(:).info];
i_info = [sc.ispn(:).info];
drop = [out.exclude_sess; out.exclude_mice(:); out.mismatched];
d_valid = ~contains(string({d_info.recName}), drop)';
i_valid = ~contains(string({i_info.recName}), drop)';
nD = numel(sc.dspn); nI = numel(sc.ispn);

for f = fieldnames(sc)'
    v = sc.(f{1});
    if isstruct(v) && numel(v) == nD && startsWith(f{1},'d')
        sc.(f{1}) = v(d_valid);
    elseif isstruct(v) && numel(v) == nI && startsWith(f{1},'i')
        sc.(f{1}) = v(i_valid);
    elseif ~isstruct(v) && size(v,1) == nD && startsWith(f{1},'d')
        sc.(f{1}) = v(d_valid,:);
    elseif ~isstruct(v) && size(v,1) == nI && startsWith(f{1},'i')
        sc.(f{1}) = v(i_valid,:);
    end
end
sc.d_valid = d_valid;
sc.i_valid = i_valid;
out.sc = sc;
end
