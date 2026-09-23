function [avg_sc, avg_sess] = get_singleCell_sess_dff(sc, rmfield, xvec)

    %% get in track 
    %rmfield = 'rm_w1' familiar trials
    sc_rm = cellfun(@(x) x.(rmfield)(:,xvec), {sc.rms}, UniformOutput=false);
    sessID = cellfun(@(x) x.recName, {sc.info}, UniformOutput=false); 
    [unique_idx, unique_sess] = findgroups(sessID);
    
    % get average activity per sess for each cell
    avg_sc = cell2mat(cellfun(@(x) mean(x(:, xvec), 'all','omitnan'), sc_rm, 'UniformOutput', false)); 
    avg_sess = NaN(1, numel(unique_sess));
    % for i = 1:numel(unique_sess)
    %     cur_sc = sc_rm(unique_idx == i)';
    %     sesstemp = mean(permute(cat(3, cur_sc{:}), [3 1 2]), 1, "omitnan");   % calc n_trial x n_bin
    % end
        


end




