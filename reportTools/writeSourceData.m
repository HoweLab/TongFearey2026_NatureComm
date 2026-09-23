function writeSourceData(xlsxFile, sheetName, blocks)
% WRITESOURCEDATA  Write several tables to one Excel tab, side by side.
%   blocks = {'Fig1D top: acceleration (LME mean, SEM)', T1;
%             'Fig1D top: acceleration per mouse',       T2; ...};
%   writeSourceData(fullfile(save2where,'PaperStats_report.xlsx'), 'Fig1', blocks)
%
% Each block gets its own columns: row 1 = block title, row 2 = column
% names, data below; one empty column between blocks. The tab is cleared
% first, other tabs (other figures) are kept.

col = 1;
for b = 1:size(blocks, 1)
    T = blocks{b, 2};
    titleCell = sprintf('%s1', colLetter(col));
    if b == 1 && isfile(xlsxFile)
        writecell(blocks(b,1), xlsxFile, 'Sheet', sheetName, 'Range', titleCell, ...
            'WriteMode', 'overwritesheet');
    else
        writecell(blocks(b,1), xlsxFile, 'Sheet', sheetName, 'Range', titleCell);
    end
    writetable(T, xlsxFile, 'Sheet', sheetName, 'Range', sprintf('%s2', colLetter(col)));
    col = col + width(T) + 1;
end
fprintf('Source data written: %s [%s], %d blocks\n', xlsxFile, sheetName, size(blocks,1));
end

function s = colLetter(n)
s = '';
while n > 0
    r = mod(n - 1, 26);
    s = [char(65 + r) s]; %#ok<AGROW>
    n = floor((n - 1) / 26);
end
end
