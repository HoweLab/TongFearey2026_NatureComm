function save_img(fig,figDir,fname)
% A function to quickly save all three figure type we need for illustrator
% work with the same filename different extension. 
% Three img type are: 
%   1) Original matlab .fig, store the figure object for easy edit.
%   2) Vector form .eps file, the format editable in illustrator
%   3) A low res .png file, for quicker layout, and preview.

% Changed by Samme: to save a bunch of images
if ischar(fname) | isstring(fname)
    fname = {fname};
end
    if numel(fig) ~= numel(fname)
        error('number of figures not matching the name')
    end
    if isscalar(fig) && ~iscell(fig)
        fig = {fig};
    end
    for i = 1:numel(fig)
        outpath = fullfile(figDir,fname{i});
        outpath = char(outpath);
        if iscell(fig{i})
            curfig = fig{i}{1};
        else
            curfig = fig{i};
        end
        set(curfig,'Renderer','painters')
        saveas(curfig,outpath,'epsc')
        saveas(curfig,outpath,'fig')
        exportgraphics(curfig,[outpath '.png'],'Resolution',150);
    end
end