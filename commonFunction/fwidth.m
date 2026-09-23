function [fieldind,thres] = fwidth(temp)


temp = smoothdata(temp,2);
[max_val,peak] = nanmax(temp);
min_val = nanmin(temp);

% try the 20% value between min and max of the average curve
thres = (max_val-min_val).* 0.2 + min_val;
% added older method here, just find below 20% of the peak value, uncomment
% to use the new threshold.
% thres = max_val.*.2;

first = find(fliplr(temp(1:peak-1))<thres,1);
last = find(temp(peak+1:end)<thres,1);

if isempty(first)
    fieldind = [1,peak+last];
elseif isempty(last)
    fieldind = [peak-first,length(temp)];
else
    fieldind = [peak-first, peak+last];
end
