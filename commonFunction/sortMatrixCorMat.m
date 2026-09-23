function [out] = sortMatrixCorMat(varargin)
% SX Note: corr run between two matrix or with itself after scaled by max
% values within each trial/bout/Rep

% maxn = data matrix after max scaling
% sortin = index to sort the rows of nerons according to which bin max activity
%       happen from small bin id to large bin (from 1 to 100)
% cm = correlation coefficient of size [numBin, numBin]
% err = distance between index of max correlation with diagonal 
            out.maxn = [];
            out.sortind = [];
            out.cm = [];
            out.err = [];
            
    switch nargin % based on # of input
        case 2
            
            dat = varargin{1};
            ip = varargin{2};
            dat_max = dat./nanmax(dat,[],2); % max scaling
            [~,ind] = nanmax(dat_max,[],2);
            [~,ind2] = sort(ind);
            
            cm = corMat(dat_max,dat_max);
            [~,icm] = max(cm,[],2); % max index of bin correspond to max correlation
            err = abs(icm-[1:length(icm)]); % max correlation wrt diagnoal
            % erro shoul dbe 0 if correlate with itself
            if ip == 1
                figure; subplot(1,2,1);
                imagesc(dat_max(ind2,:));
                caxis([0 1])

                subplot(1,2,2);

                imagesc(cm);
                caxis([0 1])

                hold on
                plot(icm,[1:length(icm)],'r','LineWidth',2);
                hold on
                plot(1:length(icm),1:length(icm),'k--')
                colormap(parula)
                caxis([0 .5])
            end
            
            out.maxn = dat_max;
            out.sortind = ind2;
            out.cm = cm;
            out.err = err;
            
            
        case 3
            
            dat = varargin{1}; % [numBouts/Trials/Spatial reps, numBins(100)]
            dat2 = varargin{2};
            ip = varargin{3};
            % seems like minmax scaling: dat normalized by max of each bout, maxVal = [numbout, 1]
            maxval = nanmax([dat,dat2],[],2);
            dat_max = dat./maxval;
            dat2_max = dat2./maxval;
            
            [~,ind] = nanmax(dat_max,[],2); % max position of each neuron
            [~,ind2] = sort(ind); % Sort the rows (bouts/trial) based on max position (early bin to late)
            % so here it calculated a cross corr of normlized dat?
            % dat_max = [numbouts/trials,nBins]
            % Two input should have same size
            cm = corMat(dat_max,dat2_max);
            [~,icm] = max(cm,[],2); % max for each row
            err = abs(icm'-[1:length(icm)]); %max correlation wrt diagnoal
            
            if ip == 1
                figure; subplot(1,3,1);
                imagesc(dat_max(ind2,:)); % this is sorting based on index of max value
                caxis([0 1])

                subplot(1,3,2);
                imagesc(dat2_max(ind2,:));
                caxis([0 1])

                subplot(1,3,3);

                imagesc(cm);
                
                hold on
                plot(icm,[1:length(icm)],'r','LineWidth',2);
                hold on
                plot(1:length(icm),1:length(icm),'k--')
                colormap(parula)
                clim([0 1]) % no negative corr. coefficient
            end            
                       
            out.maxn = dat_max;
            out.maxn2 = dat2_max;
            out.sortind = ind2;
            out.cm = cm;
            out.err = err;
            
            
    end

    % basically 
% function [cm] = corMat(dat1,dat2)
% 
%     for i = 1 : size(dat1,2)
%         for j = 1 : size(dat2,2)
%             cm(i,j) = corr(dat1(:,i),dat2(:,j),'type','Spearman','rows','complete');
%         end
%     end
    