%% SX
function exclude_sess = findfirstSession(Tpath)
    % find the first recorded session of the 2w sessions 
    T = readtable(Tpath);
    T_2w = T(contains(T.Behavior, '2 worlds'),:);
    fileN_2w = T_2w(:,1:2);
    allfirst2w = fileN_2w.FileNumber;
    firstPart = cellfun(@(x) split(x,'_'), allfirst2w, 'UniformOutput', false);
    exclude_sess = cellfun(@(x) x{2}, firstPart, 'UniformOutput', false);
    exclude_sess = strcat(fileN_2w.Mouse, '_', exclude_sess);

end