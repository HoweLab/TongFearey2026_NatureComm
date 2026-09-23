function T = sourceDataXcorr(avgout, lagScale, lagName)
% SOURCEDATAXCORR  Cross-correlation curves as plotted by plot_crosscor
% (avgout from crosscor_DI_behav): one row per lag, one column per behavior.
%   T = sourceDataXcorr(avgout_1w_onset, 1/31, 'lag_s')
%   T = sourceDataXcorr(avgout_1w, 1, 'lag_bin')

behav = fieldnames(avgout);
lags = avgout.(behav{1}).xcorr_lags(:);
T = table(lags, lags * lagScale, 'VariableNames', {'lag_samples', lagName});
for b = 1:numel(behav)
    T.([behav{b} '_xcorr']) = avgout.(behav{b}).xcorr_c(:);
end
end
