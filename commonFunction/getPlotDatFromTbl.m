function out = getPlotDatFromTbl(tbl,field,varargin)
%PLOTFROMTBL Summary of this function goes here
%   Detailed explanation goes here
% parse input to setup optional input.
ip = inputParser;
ip.addParameter('y','avgFc3')
ip.addParameter('vel_y','velocity')
% For velocity can only do random effect until the session level. When we
% get to trial level we run in to issue with numeric stability (?).
% Only for binned Location, not sure others.
ip.addParameter('randomEffects','(1|mouse) + (1|mouse:sess) + (1|mouse:sess:field)')
ip.addParameter('compare','celltype')
ip.addParameter('sessTypeVarName','sesstype')
ip.addParameter('skip',{})
ip.addParameter('xdata',[])
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

disp('Getting plot data from the table ..')
tic

if ~isempty(randomEffects)
    randomEffects = [' + ' randomEffects];
end
% num coeff same except for DI_diff_exp
num_coeff = 2;
num_v_coeff = 1;
switch compare
    case 'celltype'
        tbl.isD(tbl.celltype == 'dSPN') = 1;
        tbl.isI(tbl.celltype == 'iSPN') = 1;
        fixedEffects = {'isD', 'isI'};
        % field name for the output struct
        out_vals = {'iSPN','dSPN','diff_di','diff_id'};
        formula = {[y ' ~  1 + ' fixedEffects{1}, randomEffects],...
            [y ' ~  1 + ' fixedEffects{2} , randomEffects]};
        % calculate velocity with mixed model as well?
        tbl.calc_velmu = tbl.isD; % use half of the data/cuz there are repeats
        velmu_formula = {[vel_y ' ~ 1', randomEffects]};
    case 'novelty'
        tbl.novel(tbl.info_trial_novel == 1) = 1;
        tbl.familiar(tbl.info_trial_novel == 0) = 1;
        fixedEffects = {'novel', 'familiar'};
        % field name for the output struct
        out_vals = {'familiar','novel','diff_nf','diff_fn'};
        formula = {[y ' ~  1 + ' fixedEffects{1}, randomEffects],...
            [y ' ~  1 + ' fixedEffects{2} , randomEffects]};
        % calculate velocity with mixed model as well?
        tbl.calc_velmu(:) = 1; 
        velmu_formula = {[vel_y ' ~ 1 ', randomEffects]};
    case 'sesstype'
        tbl.gui(tbl.(sessTypeVarName) == 'gui') = 1;
        tbl.track(tbl.(sessTypeVarName) == 'track') = 1;
        fixedEffects = {'gui','track'};
        % field name for the output struct
        out_vals = {'track','gui','diff_gt','diff_tg'};
        formula = {[y ' ~  1 + ' fixedEffects{1}, randomEffects],...
            [y ' ~  1 + ' fixedEffects{2}, randomEffects]};
        % calculate velocity with mixed model as well?
        tbl.calc_velmu(:) = 1; 
        velmu_formula = {[vel_y ' ~ 1', randomEffects]};
    case 'singleType' % if not get difference
        out_vals = {'exp1'};
        formula = {[y ' ~  1 ', randomEffects]};
        tbl.calc_velmu(:) = 1; 
        velmu_formula = {[vel_y ' ~ 1', randomEffects]};
        num_coeff = 1;
        num_v_coeff = 1;
    case 'DI_diff_exp'
        % avgFc3 ~ 1 + isD*expType + (1|mouse) + (1|mouse:sess) + (1|mouse:sess:field)'
        fixedEffects = {'isD', 'isI'};
        interaction = {'exp0','exp1'}; % two types of exp
        % ispn/dspn = estimate for exp 0; diff_di/diff_id = estimate for
        % DI/ID diff for exp 0; diff_DI_exp01 = DI_diff for exp0 - exp1
        out_vals = {'iSPN','dSPN','diff_di','diff_id',...
                    'diff_I_exp01','diff_D_exp01','diff_DI_exp01','diff_ID_exp01' };
        formula = {[y ' ~  1 + ' fixedEffects{1} '*' interaction{1}, randomEffects],...
            [y ' ~  1 + ' fixedEffects{2} '*' interaction{1}, randomEffects]};
        % calculate velocity with mixed model 
        tbl.calc_velmu(:) = 1; 
        velmu_formula = {[vel_y ' ~ 1 + ' interaction{1}, randomEffects], ...
            [vel_y ' ~ 1 + ' interaction{2}, randomEffects]}; % vel for both exp
        num_coeff = 4;
        num_v_coeff =2;
end
% preallocate
if isempty(xdata)
    start_idx = min(tbl{:,field});
    end_idx = max(tbl{:,field});
    all_bins = start_idx:end_idx;
    xdata = start_idx:end_idx;
else
    all_bins = xdata(1,:);
    xdata = xdata(2,:);
end



% mu, sem, pval size depends on number for model to fit and number of
% time/location bins.
% Shuold look something like
% [out_val{1}.sem(1), out_val{1}.sem(2) .... out_val{1}.sem(n);
%  out_val{2}.sem(1), out_val{2}.sem(2) .... out_val{2}.sem(n);
%  ....
%  out_val{m}.sem(1), out_val{1}.sem(2) .... out_val{m}.sem(n);]\
% if we have m output variable and n timepoints
mu = nan(numel(out_vals),numel(all_bins));
sem = nan(numel(out_vals),numel(all_bins));
pval = nan(numel(out_vals),numel(all_bins));
rnd = cell(numel(formula),numel(all_bins));


% lets initiate separate array for velocity calculation
vmu = nan(num_v_coeff * numel(velmu_formula),numel(all_bins));
vsem = nan(num_v_coeff * numel(velmu_formula),numel(all_bins));
vpval = nan(num_v_coeff * numel(velmu_formula),numel(all_bins));
% do a parallel for loop so that things can be done faster
% This will speed thing up by around 4 fold if I have parallel processing
% on already.
y_var = y;
vel_var = vel_y;
parfor i = 1:numel(all_bins)
    %disp(['Bin number ' num2str(i)])
    cb = all_bins(i);
    subtbl = tbl(tbl.(field) == cb,:);
    if height(subtbl) > 2 & ~isempty(y_var) 
        % fit model for current timepoint.
        [mu(:,i), sem(:,i), pval(:,i), rnd(:,i)] = fit_for_subtbl(subtbl,formula,num_coeff);        
    end
    if height(subtbl) > 2 & ~isempty(vel_var) % only do behavior if vel_y is not empty
        % now do another subtbl to calculate velocity
        subtbl = subtbl(subtbl.calc_velmu == 1,:);        
        [vmu(:,i), vsem(:,i), vpval(:,i)] = fit_for_subtbl(subtbl,velmu_formula,num_v_coeff);
    end
end
% after running the parallel for loop we can sort output to the struct.
out = struct();
for ck = 1:numel(out_vals)
    out.(out_vals{ck}).mu = mu(ck,:);
    out.(out_vals{ck}).sem = sem(ck,:);
    out.(out_vals{ck}).pval = bonf_holm(pval(ck,:),0.05);
    out.(out_vals{ck}).pval_uncorrected = pval(ck,:);
end
% random effect only had 1 value, not sure what that is, give an output
% first.
out.rnd = rnd;

out.velocity.mu = vmu;
out.velocity.sem = vsem;
out.velocity.pval_uncorrected = vpval;
for i = 1:size(vpval, 1)
    out.velocity.pval(i,:) = bonf_holm(vpval(i,:),0.05);
end

out.bins = all_bins;
out.xdata = xdata;
disp(['Time taken ' num2str(toc) 's'])
end


function [mu, sem, pval, rnd] = fit_for_subtbl(subtbl,formula,num_coeff)
% num slope should be something like number of fixed effect -1
nf = numel(formula);
numr = nf * num_coeff; % for 2 coeffs
mu = nan(numr,1);
sem = nan(numr,1);
pval = nan(numr,1);
rnd = cell(nf,1);
for cf = 1:numel(formula)
    %%
    try
        glme = fitglme(subtbl,formula{cf});
        mu(cf) = glme.Coefficients.Estimate(1);
        sem(cf) = glme.Coefficients.SE(1);
        pval(cf) = glme.Coefficients.pValue(1);
        [~, bnames, bstats] = randomEffects(glme);
        rnd{cf} = bstats(matches(bnames.Group,'mouse'),:);
        % allow this function to store vary how many elements in the output
        % array.
        for i = 2:num_coeff % fill the rest
            idx = (i-1)*nf + cf;
            mu(idx) = glme.Coefficients.Estimate(i);
            sem(idx) = glme.Coefficients.SE(i);
            pval(idx) = glme.Coefficients.pValue(i);
        end
    
    catch
        try % if glme not working, do lme
            %disp('- current bin using lme instead of glme')
            glme =  fitlme(subtbl,  formula{cf});
            mu(cf) = glme.Coefficients.Estimate(1);
            sem(cf) = glme.Coefficients.SE(1);
            pval(cf) = glme.Coefficients.pValue(1);
            [~, bnames, bstats] = randomEffects(glme);
            rnd{cf} = bstats(matches(bnames.Group,'mouse'),:);
            % allow this function to store vary how many elements in the output
            % array.
            for i = 2:num_coeff % fill the rest
                idx = (i-1)*nf + cf;
                mu(idx) = glme.Coefficients.Estimate(i);
                sem(idx) = glme.Coefficients.SE(i);
                pval(idx) = glme.Coefficients.pValue(i);
            end

        catch ME
            %disp(ME.message)
            % if glme didn't work then we fill where suppose to have value with
            % nan.
            mu(cf) = nan;
            sem(cf) = nan;
            pval(cf) = nan;
            for i = 2:num_coeff % fill the rest
                idx = (i-1)*nf + cf;
                mu(idx) = nan;
                sem(idx) = nan;
                pval(idx) = nan;
            end
        end
    end
end
end

