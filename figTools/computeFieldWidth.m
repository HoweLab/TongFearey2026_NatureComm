function [fieldWidth, bounds] = computeFieldWidth(sc, varargin)
% COMPUTEFIELDWIDTH  Place-field width of each neuron, from its trial-averaged
% rate maps (as in paperFigScripts/paperFigureSingleCell_Inf.m).
%   dFieldWidth = computeFieldWidth(dat.sc.dspn);          % [nCell x 2]
%   [fw, bd] = computeFieldWidth(ispn, 'scale', [1.5 1]);
%
% For every neuron and every map, the field bounds come from fwidth (first
% bins on each side of the peak that fall below 20% of the peak-to-trough
% range of the smoothed map); the width is their difference. Neurons whose
% map has no usable field are left as NaN.
%
% Column 1 is the first map (bout distance/time for the infinite track,
% Track 1 for the two-world task), column 2 the second (visual position /
% Track 2), matching dstb / istb.
%
% Name-value options:
%   maps  - rate-map fields of sc(i).rms (default {'rm_w1','rm_w2'})
%   scale - cm per bin for each map (default [1.5 1], as in the infinite
%           track scripts: distance maps are 1.5 cm bins, the visual-position
%           maps are left in bins)

ip = inputParser;
ip.addParameter('maps', {'rm_w1','rm_w2'});
ip.addParameter('scale', [1.5 1]);
ip.parse(varargin{:});
o = ip.Results;

nCell = numel(sc);
nMap = numel(o.maps);
bounds = nan(nCell, 2, nMap);
for m = 1:nMap
    fld = o.maps{m};
    bd = nan(nCell, 2);
    parfor n = 1:nCell
        temp = nanmean(sc(n).rms.(fld)); %#ok<NANMEAN> % mean over trials/bouts
        try
            bd(n,:) = fwidth(temp);
        catch
            % no field found for this neuron: stays NaN
        end
    end
    bounds(:,:,m) = bd;
end

fieldWidth = nan(nCell, nMap);
for m = 1:nMap
    fieldWidth(:,m) = (bounds(:,2,m) - bounds(:,1,m)) * o.scale(m);
end
end
