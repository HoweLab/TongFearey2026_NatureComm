function [cat_dat, grp_idx, tbl_name] = find_cells_summaryROIs(tbl, fpath, tracetp, skiplogical, varargin)
%% fpath: directory to summaryROI
    
    ip = inputParser;
    ip.addParameter('skipHighFreq',1)
    ip.addParameter('sr', 31)
    ip.parse(varargin{:});
    for j=fields(ip.Results)'
        eval([j{1} '=ip.Results.' j{1} ';']);
    end
    %fpath  = 'D:\Lily_ProcessedData\SummaryROI\linearTrackRecording_24to25_1w_deconv';
    filelist = cat(1,dir(fullfile(fpath,'*summaryROIs.mat')));
    roi_names = {filelist.name}';
    [grp_idx, tbl_name] = findgroups(tbl(:, {'fname'}));
    dat_file = unique(tbl.fname);
    base2 = regexprep(roi_names, '_[^_]+$', '');
    [tf, ~] = ismember(base2, dat_file);
    roi_names = roi_names(tf);
    base2 = base2(tf);
    if ~any(matches(base2, dat_file))
        error('file not match, double check')
    end
    
    cat_dat = cell(height(tbl),1); % size of the table
    for f = 1:numel(dat_file)
        fprintf('Processing %s \n', base2{f});
        dat = load(fullfile(fpath, roi_names{f})); % should open
        if skiplogical % if exclude cells
            if isfield(dat.dSPN, 'skipLogical')
                didx = dat.dSPN.skipLogical;
                iidx = dat.iSPN.skipLogical;
            else
                didx = true(size(dat.dSPN.(tracetp),2),1)'; 
                iidx = true(size(dat.iSPN.(tracetp),2),1)';                 
            end
        else
            didx = true(size(dat.dSPN.(tracetp),2),1)'; 
            iidx = true(size(dat.iSPN.(tracetp),2),1)'; 
        end
        
        dSPN = dat.dSPN.(tracetp)(:,didx);
        iSPN = dat.iSPN.(tracetp)(:,iidx);
        if skipHighFreq % skip high frequency firing neurons
            % out_dspn = get_event(dat.dSPN.Fc, dat.dSPN.Fc3', sr);
            % out_ispn = get_event(dat.iSPN.Fc, dat.iSPN.Fc3, sr);
        end
        if ~sum(grp_idx == f & tbl.isD)== size(dSPN,1)
            warning('Frames number not matching for file %s',base2{f});
        end
        % Fc3 
        cat_dat(grp_idx == f & tbl.isD,:) = num2cell(dSPN, 2);
        cat_dat(grp_idx == f & tbl.isI,:) = num2cell(iSPN, 2);    
    end

end % end of function

function [out] = get_event(dffTrace, eventTrace, fs)

% dffTrace: original dF/F trace, 1 x nFrames
% eventTrace: thresholded event trace, 0 where no event
% fs: frame rate
% eventID: which detected event to fit
for i = 1:size(dffTrace, 2)
    eventFrames = eventTrace ~= 0;
    d = diff([false eventFrames false]);
    
    onsets = find(d == 1);
    offsets = find(d == -1) - 1;
    
    onset = onsets;
    offset = offsets;
    
    % Find peak within detected event using original dF/F
    [~, localPeak] = max(dffTrace(onset:offset));
    peakFrame = onset + localPeak - 1;
    
    % Fit from peak to either event offset or a bit beyond
    fitEnd = min(numel(dffTrace), peakFrame + round(2 * fs)); % 2-s window
    fitFrames = peakFrame:fitEnd;
    
    y = dffTrace(fitFrames)';
    t = ((0:numel(fitFrames)-1)' ./ fs);
    
    % Optional: remove NaNs
    valid = ~isnan(y);
    y = y(valid);
    t = t(valid);
    
    % Fit y = A*exp(-t/tau) + c
    ft = fittype('A*exp(-t/tau) + c', ...
        'independent', 't', ...
        'coefficients', {'A','tau','c'});
    
    A0 = max(y) - min(y);
    tau0 = 0.7;
    c0 = min(y);
    
    opts = fitoptions(ft);
    opts.StartPoint = [A0, tau0, c0];
    opts.Lower = [0, 0.05, -Inf];
    opts.Upper = [Inf, 10, Inf];
    
    fitObj = fit(t, y, ft, opts);
    
    tau = fitObj.tau;
    
    % event width
    evt_with = (offsets - onsets)./fs;
    evt_freq = length(offsets)./length(onsets)./fs;
    prct_evt = mean(eventFrames, 'omitnan');
    
    out.tau = tau;
    out.event_width = evt_with;
    out.event_freq = evt_freq;
    out.prct_evt = prct_evt;
    
end
end

