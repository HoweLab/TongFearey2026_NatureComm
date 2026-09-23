function [] = boxViolin(varargin) 
%boxScatter(input1,input2,input3,inputN,figureID)

%Inputs: 
%n comma separated data vectors (i.e. boxScatter(input1,input2,input3,inputN,figureID)
% last input is an integer for figure # (useful when subplotting)
% e.g. if want to plot 2 inputs on figure 5, call boxScatter(input1,input2,5)


fid = varargin{end};
figure(fid);
varargin(end) = [];
n = length(varargin);

% Define x locations of plot
xp = [1];
if n>1
    for i = 2 : n
        xp = [xp, xp(i-1)+1];
    end
end

% % Generate scatter around those xp locations for each data input
    boxDat = []; boxG = [];
    for i = 1 : n
        clear xt;
        xt = repelem(xp(i),length(varargin{i}))';
        figure(fid); %subplot(1,5,1); 
        hold on;
        % jw = 0.5 * min(diff(unique(varargin{i})));
        swarmchart(xt,varargin{i},10,[.2 .2 .2],'filled')%,'XJitterWidth',jw
        % scatter(xt, varargin{i},10,repmat([.2 .2 .2],length(xt),1),'filled')
        
        boxDat = [boxDat; varargin{i}];
        boxG = [boxG; ones(length(varargin{i}),1).*i];
        
    end
    boxplot(boxDat,boxG,'plotstyle','compact'); axis square;
    p = ranksum(boxDat(boxG==1),boxDat(boxG==2));
    title(['G1 and G2 significant, p = ' num2str(p)])

    delete(findobj(gca,'Tag','Outliers')) % to not show outliers since we are already plotting scatters in the swarmchart
    hold on
    y = axis;



end

