
function [grpsN, sort_idx, tidN] = sort_tbl(tbl_novel, groupvars,ascendORdescend)
    
    [grpsN,tidN] = findgroups(tbl_novel(:,groupvars));
    sort_mat = str2double(replace(string(tidN{:,:}),'fld','')); 
    [~, sort_idx] = sortrows(sort_mat,ascendORdescend);
end