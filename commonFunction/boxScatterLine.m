function [] = boxScatterLine(varargin) 

fid = varargin{end};
figure(fid); clf(fid)
varargin(end) = [];
n = length(varargin);

% Define x locations of plot
xpLookUp = [1.25 1.75 3.25 3.75 5.25 5.75 7.25 7.75 9.25 9.75];
xp = xpLookUp(1:n);   

% % Generate scatter around those xp locations for each data input
    boxDat = []; boxG = [];
    for i = 1 : n
        clear xt;
        xt = xp(i);
        figure(fid); hold on;
        scatter(repmat(xt,length(varargin{i}),1), varargin{i},10,repmat([.2 .2 .2],length(xt),1),'filled')
        
        boxDat = [boxDat; varargin{i}];
        boxG = [boxG; ones(length(varargin{i}),1).*i];
        
    end
    boxplot(boxDat,boxG,'plotstyle','compact'); axis square;
    hold on
    y = axis;
    
    for j = 1 : 2: n-(mod(n,2))
    plot([repmat(xp(j),length(varargin{j}),1) repmat(xp(j+1),length(varargin{j+1}),1)]',...
        [varargin{j} varargin{j+1}]','color',[.5 .5 .5])
    end
    
    % plot([y(1) y(2)],[0 0],'r-')


% if ifLines == 1
%     for i = 1 : n
%         clear xt;
%         xt = xp(i) + .02.*randn(length(varargin{i}),1);
%         figure(99); hold on;
% %         scatter(xt, varargin{i},10,repmat([.2 .2 .2],length(xt),1),'filled')
%         plot(xt, varargin{i},'color',[.5 .5 .5])
% 
%         
% %         boxDat = [boxDat; varargin{i}];
% %         boxG = [boxG; ones(length(varargin{i}),1).*i];
%         
%     end
% end




end
