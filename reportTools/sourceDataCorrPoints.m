function T = sourceDataCorrPoints(cor, varname, group, pThresh)
% SOURCEDATACORRPOINTS  Every point of a boxViolin_colorSig correlation plot.
%   T = [sourceDataCorrPoints(track_corr.dspn_cor, 'binned_v_eachN_cor', 'dSPN', 0.001)
%        sourceDataCorrPoints(track_corr.ispn_cor, 'binned_v_eachN_cor', 'iSPN', 0.001)];
%
% cor is a struct array from encoding_dff_vel (one element per imaging
% field). Rows are in the same order as grab_celldat concatenates them.
% Columns: Group, Field (mouse_session_field; "session_<i>" if the save has
% no sessName), Neuron (index within the
% field; NaN for per-field correlations), rho, p, Significant (p < pThresh,
% plotted in red; otherwise gray).

[fieldName, neuron, rho, p] = deal({});
for i = 1:numel(cor)
    tho = cor(i).(varname).tho(:);
    n = numel(tho);
    if isfield(cor, 'sessName')
        nm = string(cor(i).sessName);
        nm = nm(1);
    else % older saves without sessName: index of the struct element
        nm = "session_" + i;
    end
    fieldName{i} = repmat(nm, n, 1);
    if contains(varname, 'eachN')   % per-neuron correlation
        neuron{i} = (1:n)';
    else                            % one population value per field
        neuron{i} = nan(n, 1);
    end
    rho{i} = tho;
    p{i} = cor(i).(varname).p(:);
end
rho = vertcat(rho{:});
p = vertcat(p{:});
T = table(repmat(string(group), numel(rho), 1), vertcat(fieldName{:}), ...
    vertcat(neuron{:}), rho, p, p < pThresh, ...
    'VariableNames', {'Group','Field','Neuron','rho','p', ...
    sprintf('Significant_p_lt_%g', pThresh)});
T.Properties.VariableNames{end} = matlab.lang.makeValidName(T.Properties.VariableNames{end});
end
