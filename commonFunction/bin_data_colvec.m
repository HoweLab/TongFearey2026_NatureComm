function dat_mat = bin_data_colvec(dat, cur_bin, bin_vals)
% Inputs:
%   dat      : vector or matrix with first dimension as samples
%   cur_bin  : vector of bin labels, same length as dat along first dimension
%   bin_vals : bin values to search for, e.g. 1:100
%
% Output:
%   dat_mat  : [nTrial x nBin]

    cur_bin = cur_bin(:);
    bin_vals = bin_vals(:)';

    if size(dat,1) ~= numel(cur_bin)
        error('dat must have the same number of rows as cur_bin has elements.');
    end

    nBin = numel(bin_vals);

    % Detect new trial when bin sequence resets/decreases
    newTrial = [true; diff(cur_bin) < 0];

    % Assign trial ID to each sample
    trial_id = cumsum(newTrial);
    nTrial = max(trial_id);

    % Keep only valid bins and non-NaN data
    [isValidBin, bin_col] = ismember(cur_bin, bin_vals);

    if isvector(dat)
        dat = dat(:);
        good = isValidBin & ~isnan(dat);

        dat_mat = accumarray( ...
            [trial_id(good), bin_col(good)], ...
            dat(good), ...
            [nTrial, nBin], ...
            @mean, ...
            NaN);

    else
        % If dat has multiple columns, output is [nTrial x nBin x nFeature]
        nFeature = size(dat,2);
        dat_mat = NaN(nTrial, nBin, nFeature);

        for f = 1:nFeature
            x = dat(:,f);
            good = isValidBin & ~isnan(x);

            dat_mat(:,:,f) = accumarray( ...
                [trial_id(good), bin_col(good)], ...
                x(good), ...
                [nTrial, nBin], ...
                @mean, ...
                NaN);
        end
    end
end