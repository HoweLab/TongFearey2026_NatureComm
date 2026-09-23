function [tbl, info]= mkTblForPlot(varargin)
%MKTBLFORPLOT Summary of this function goes here
%   Detailed explanation goes here
ip = inputParser;
ip.addRequired('fpath',@(x) isfile(x))
ip.addRequired('conditions',@(x) iscell(x) & width(x) == 2 | isempty(x))
% for trial and bout selection should pass in a function that act on the
% trial_info/bout_info struct, and returns a boolean index.
ip.addParameter('skip_trial',@(x) logical(x.trial_skipped));
ip.addParameter('skip_bout',@(x) logical(x.bout_rewarded));
ip.addParameter('skip_repeat',@(x) logical(x.repeat_skipped));
ip.addParameter('trial_fields',{'world','trialTimes','trial_novel'})
ip.addParameter('bout_fields',{'bout_intrack','bout_rewarded','peak_vel','bout_length'})
ip.addParameter('repeat_fields',{'repeatTimes','repeat_skipped'})
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

disp('Making Tables')
tic
% load data from path
all_dat = load(fpath);
mobj_fields = fieldnames(all_dat);
dat_field = mobj_fields{contains(mobj_fields,'data')}; % all data
info = all_dat.(mobj_fields{contains(mobj_fields,'info')}); % all info
dat = all_dat.(dat_field);

% drop the empty rows
info_erow = cellfun(@(x)isempty(x), {info.mouse});
dat_erow = cellfun(@(x)isempty(x), {dat.bout_info});
if all(find(info_erow) == find(dat_erow))
    disp('Drop empty rows')
    info(info_erow) = [];
    dat(dat_erow) = [];
end

% assume we include everything if no conditions are given
ind = 1:numel(info); % nu,ber of sessions
% go through each condition, find indices that satisfy it, and take the
% intersect of all conditions.
for c = 1:length(conditions)
    curcond = conditions{c,2};
    if ischar(curcond) || iscellstr(curcond)
        tblcond = {info(:).(conditions{c,1})};
        tblcond(cellfun(@isempty, tblcond)) = {'NaN'};
        ind = intersect(ind, find(matches(tblcond,curcond)));
    elseif isnumeric(curcond) || islogical(curcond)
        tblcond = [info(:).(conditions{c,1})];
        ind = intersect(ind, find(tblcond == curcond));
    end
end

% start making table from sessions that satsify the condition.
cat_tbl = table();
for i = ind
    %cur_dat = matobj.(dat_field)(1,i);
    cur_dat = dat(1,i);
    tbl = cur_dat.tbl;
    tbl.mouse(:) = categorical({info(i).mouse});
    tbl.sess(:) = categorical({info(i).session});
    tbl.field(:) = categorical({info(i).field});
    tbl.exp_stage(:) = categorical({info(i).exp_stage});
    tbl.fname(:) = {[info(i).mouse '_' info(i).session '_' info(i).field]};
    if isfield(info(i), 'trackLength')
        tbl.trackLgth(:) = [info(i).trackLength];
    end
    % there are some non-integer trials (probabaly due to averaging during
    % split_n_bin). It might make most sense to drop out these values since
    % a quick investigation reveal that other fields (e.g. y_pos and phase)
    % are pretty messed up as well.
    % This will clean the nan trials at the begining and end of a session
    % too.
    if ismember('trial',tbl.Properties.VariableNames)
        % tbl(mod(tbl.trial,1) ~=0,:) =[]; % integer test, remove non in rows
        tbl.trial = round(tbl.trial);
        tbl.skip_trial = ismember(tbl.trial, ...
            cur_dat.trial_info.trial(skip_trial(cur_dat.trial_info)));
        for fld = trial_fields
            tbl_varname = ['info_' fld{1}];
            tbl.(tbl_varname)(:) = nan;
            for t = 1:numel(cur_dat.trial_info.trial)
                idx = tbl.trial == cur_dat.trial_info.trial(t);
                tbl{idx,tbl_varname} = repelem(cur_dat.trial_info.(fld{1})(t), ...
                    sum(idx))';
            end
        end
        tbl = convertvars(tbl,{'trial'},'categorical');
        info(i).trial_info = cur_dat.trial_info;
    else
        tbl.trial(:) = categorical(nan);
        % tbl.skip_trial(:) = false;
        % for fld = trial_fields
        %     tbl_varname = ['info_' fld{1}];
        %     tbl.(tbl_varname)(:) = nan;
        % end
    end
    if ismember('bout',tbl.Properties.VariableNames)
        tbl.skip_bout = ismember(tbl.bout, ...
            cur_dat.bout_info.bout(skip_bout(cur_dat.bout_info)));
        for fld = bout_fields
            tbl_varname = ['info_' fld{1}];
            tbl.(tbl_varname)(:) = nan;
            for b = 1:numel(cur_dat.bout_info.bout)
                idx = tbl.bout == cur_dat.bout_info.bout(b);
                tbl{idx,tbl_varname} = repelem(cur_dat.bout_info.(fld{1})(b), ...
                    sum(idx))';
            end
        end
        tbl = convertvars(tbl,{'bout'},'categorical');
        info(i).bout_info = cur_dat.bout_info;
    else
        tbl.bout(:) = categorical(nan);
        % tbl.skip_bout(:) = false;
        % for fld = bout_fields
        %     tbl_varname = ['info_' fld{1}];
        %     tbl.(tbl_varname)(:) = nan;
        % end
    end
    if ismember('repeat',tbl.Properties.VariableNames)
        tbl.skip_repeat = ismember(tbl.repeat, ...
            cur_dat.trial_info.repeat(skip_repeat(cur_dat.trial_info)));
        for fld = repeat_fields
            tbl_varname = ['info_' fld{1}];
            tbl.(tbl_varname)(:) = nan;
            for b = 1:numel(cur_dat.trial_info.repeat)
                idx = tbl.repeat == cur_dat.trial_info.repeat(b);
                tbl{idx,tbl_varname} = repelem(cur_dat.trial_info.(fld{1})(b), ...
                    sum(idx))';
            end
        end
        tbl = convertvars(tbl,{'repeat'},'categorical');
        info(i).trial_info = cur_dat.trial_info;
    else
        tbl.repeat(:) = categorical(nan);

    end
    % The last step, convert variable to categorical and cat current tbl to
    % the overal cat_tbl.
    cat_tbl(height(cat_tbl)+1:height(cat_tbl)+height(tbl),:) = tbl;
end

% remove info that we didn't werent through
info = info(ind);

% after concatenating data, add a column that gives a unique index of
% bout/trial
mstb = cat_tbl{:,{'mouse','sess','trial','bout'}};
[~, ~, utid] = unique(mstb(~ismissing(cat_tbl.trial),[1 2 3]),'row');
cat_tbl.utid(~ismissing(cat_tbl.trial)) = utid;
[~, ~, ubid] = unique(mstb(~ismissing(cat_tbl.bout),[1 2 4]),'row');
cat_tbl.ubid(~ismissing(cat_tbl.bout)) = ubid;

% reform table, make celltype a column.
varnames = cat_tbl.Properties.VariableNames;
others = varnames(~contains(varnames,'SPN'));
d_vars = varnames(contains(varnames,'DSPN'));
tblD = [cat_tbl(:,d_vars) cat_tbl(:,others)];
tblD = renamevars(tblD, d_vars, replace(d_vars,'DSPN','Fc3'));
tblD.isD(:) = 1;
tblD.isI(:) = 0;
tblD.celltype(:) = categorical({'dSPN'});
i_vars = varnames(contains(varnames,'ISPN'));
tblI = [cat_tbl(:,i_vars) cat_tbl(:,others)];
tblI = renamevars(tblI, i_vars, replace(i_vars,'ISPN','Fc3'));
tblI.isD(:) = 0;
tblI.isI(:) = 1;
tblI.celltype(:) = categorical({'iSPN'});
tbl = [tblD; tblI];


disp(['Time taken ' num2str(toc) 's'])

end

