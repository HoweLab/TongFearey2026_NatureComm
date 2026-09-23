function parsed  = readSwitchInfo(Tpath)
% Tpath: path to the table containing swicth information
    T = readtable(Tpath);
    mouse = unique(T.Mouse);
    fileN = T.FileNumber;
    parsed = struct;
    nov_files = matches(T.Behavior, '2 worlds');
    fam_files = matches(T.Behavior, '1 world');
    % Get for each mouse
    for i = 1:numel(mouse)
        tempIdx = matches(T.Mouse, mouse{i});
        novF = fileN(nov_files & tempIdx);
        partstemp = strsplit(novF{:}, '_');
        parsed.(mouse{i}).nov = partstemp{2};
        famF= fileN(fam_files & tempIdx);
        partstemp = strsplit(famF{:}, '_');
        parsed.(mouse{i}).fam = partstemp{2};
    end
end