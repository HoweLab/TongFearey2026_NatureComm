function out = selectNovelFamiliar(varargin)
% SELECTNOVELFAMILIAR  Trials from the first exposure to a novel track and
% from the familiar (expert) track, one session per mouse (Fig 4).
%   dat = selectNovelFamiliar();
%   dat = selectNovelFamiliar('numSessNov',1, 'mice_excluded',{'tdt5'});
%
% For each mouse in the switch-info sheet:
%   novel    - the first two-world session(s) with novel-track trials
%   familiar - one-world session(s) before the switch session, taken from
%              the end of training (the second-latest session)
% Mice without both conditions are skipped.
%
% Output fields:
%   tbl       - trials of both conditions stacked (switch_tbl)
%   fam / nov - logical index into tbl (familiar / novel trials)
%   binned    - tbl summarised per mouse/session/field/trial/vuBin/celltype
%               (mean avgFc3, velocity, acceleration) with info_trial_novel
%   afterNov  - trials of the session after the novel one, per mouse
%   perMouse  - what was picked per mouse (sessions, world, trial range, n)
%
% Name-value options (defaults):
%   fpath          - all-sessions track file
%   switchInfoPath - novel/familiar switch info sheet
%   numSessNov / numSessFam - sessions per condition (1)
%   mice_excluded  - mice to drop afterwards ({})
%   verbose        - print the per-mouse summary (true)

ip = inputParser;
ip.addParameter('fpath', '');
ip.addParameter('switchInfoPath', '');
ip.addParameter('numSessNov', 1);
ip.addParameter('numSessFam', 1);
ip.addParameter('mice_excluded', {});
ip.addParameter('verbose', true);
ip.parse(varargin{:});
o = ip.Results;

conditions = {'exp_stage',{'1 world','2 world'};'sess_type','track';'mismatch',0};
track_tbl = mkTblForPlot(o.fpath,conditions);
parsed  = readSwitchInfo(o.switchInfoPath);  % switch session per mouse
allmice = fieldnames(parsed);
groupvars = {'sess'};

switch_tbl = cell2table(cell(0,width(track_tbl)), ...
    'VariableNames',track_tbl.Properties.VariableNames);
afterNov_tbl = switch_tbl;
perMouse = table();

for m = 1:numel(allmice)
    mouse = string(allmice(m));
    fam_sess = categorical(cellstr(parsed.(mouse).fam));
    % novel: two-world trials in the novel track
    tbl_novel = track_tbl(~track_tbl.skip_trial ...
        & track_tbl.exp_stage == '2 world' ...
        & track_tbl.mouse == mouse  ...
        & track_tbl.info_trial_novel==1 ,:);
    % familiar: one-world trials, excluding the switch session
    tbl_fam = track_tbl(~track_tbl.skip_trial ...
        & track_tbl.exp_stage == '1 world' ...
        & track_tbl.mouse == mouse  ...
        & track_tbl.info_trial_novel==0  ...
        & track_tbl.sess ~= fam_sess, :);
    if isempty(tbl_novel) || isempty(tbl_fam)
        if o.verbose, fprintf('Mouse %s has no novel or familar data \n', mouse); end
        continue
    end

    [grpsN, sort_idxN, tidN] = sort_tbl(tbl_novel, groupvars,'ascend');  % earliest novel sessions
    [grpsF, sort_idxF, tidF] = sort_tbl(tbl_fam, groupvars,'descend');   % latest familiar sessions
    Ntrial_idx = sort_idxN(1:min(o.numSessNov,length(sort_idxN)))';
    nov_idx = any(grpsN == Ntrial_idx,2);
    switch_tbl = cat(1,switch_tbl,tbl_novel(nov_idx,:));

    % the session right after the novel session used
    if numel(sort_idxN) >= 2
        Nafter_idx = sort_idxN(2:min(2,length(sort_idxN)))';
        afterNov_tbl = cat(1,afterNov_tbl,tbl_novel(any(grpsN == Nafter_idx,2),:));
    end

    Ftrial_idx = sort_idxF(min(2,length(sort_idxF)): ...
        min(o.numSessFam-1+min(2,length(sort_idxF)),length(sort_idxF)))';
    fam_idx = any(grpsF == Ftrial_idx,2);
    switch_tbl = cat(1,switch_tbl,tbl_fam(fam_idx,:));

    % what was picked for this mouse
    temp_tbl = cat(1, tbl_novel(nov_idx,:), tbl_fam(fam_idx,:));
    cur_nov = temp_tbl.info_trial_novel==1;
    cur_fam = ~cur_nov;
    [~, ~, numNov] = sort_tbl(temp_tbl(cur_nov,:), {'sess','field','trial'},'ascend');
    [~, ~, numFam] = sort_tbl(temp_tbl(cur_fam,:), {'sess','field','trial'},'ascend');
    perMouse = [perMouse; table(mouse, ...
        strjoin(string(unique(tidN.sess(Ntrial_idx))), ' '), ...
        unique(temp_tbl.info_world(cur_nov)), ...
        min(double(temp_tbl.trial(cur_nov))), max(double(temp_tbl.trial(cur_nov))), height(numNov), ...
        strjoin(string(unique(tidF.sess(Ftrial_idx))), ' '), ...
        unique(temp_tbl.info_world(cur_fam)), ...
        min(double(temp_tbl.trial(cur_fam))), max(double(temp_tbl.trial(cur_fam))), height(numFam), ...
        'VariableNames', {'mouse','novel_sess','novel_world','novel_firstTrial', ...
        'novel_lastTrial','novel_nTrials','fam_sess','fam_world','fam_firstTrial', ...
        'fam_lastTrial','fam_nTrials'})]; %#ok<AGROW>
end

keep = ~contains(cellstr(switch_tbl.mouse), o.mice_excluded);
if isempty(o.mice_excluded), keep = true(height(switch_tbl),1); end
out.tbl = switch_tbl(keep,:);
out.fam = out.tbl.info_trial_novel == 0;
out.nov = out.tbl.info_trial_novel == 1;
out.afterNov = afterNov_tbl;
out.perMouse = perMouse;
out.binned = groupsummary(out.tbl, ...
    {'mouse','sess','field','trial','vuBin','isD','isI','celltype','info_trial_novel'}, ...
    'mean',{'avgFc3','velocity','acceleration'});

if o.verbose
    disp(perMouse)
    fprintf('Novel: %d trials; familiar: %d trials (%d mice)\n', ...
        numel(unique(out.tbl.utid(out.nov))), numel(unique(out.tbl.utid(out.fam))), ...
        numel(unique(out.tbl.mouse)));
end
end
