function [cm] = corMat(dat1,dat2)
    cm = nan(size(dat1,2),size(dat2,2));
    for i = 1 : size(dat1,2) % two for loop should have same size
        for j = 1 : size(dat2,2)
            cm(i,j) = corr(dat1(:,i),dat2(:,j),'type','Spearman','rows','complete');
        end
    end


