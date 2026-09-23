function n = countCorrSampleSize(cor, varname, tbl)
    % COUNTCORRSAMPLESIZE  Sample sizes for a correlation struct array from
    % encoding_dff_vel (one element per imaging field, e.g. dspn_cor).
    %   n = countCorrSampleSize(track_corr.dspn_cor, 'binned_v_eachN_cor', ...
    %                           track_tbl(trial_selection,:))
    %
    % cor(i).(varname).tho holds one value per neuron (per-neuron correlation)
    % or one value per field (population correlation). tbl is the table the
    % correlations were computed from; it is only used to map each field
    % (cor(i).sessName = mouse_session_field) back to its mouse and session.
    % Returns a struct usable as the 'n' input of statsRow / reportGroupCompare.
    
    tho = arrayfun(@(c) c.(varname).tho(:), cor, 'UniformOutput', false);
    nPerField = cellfun(@(x) sum(~isnan(x)), tho);
    used = nPerField > 0;
    
    isPerNeuron = any(cellfun(@numel, tho) > 1);
    n.nUnits = sum(nPerField);
    if isPerNeuron
        n.unit = "neurons";
        n.nNeurons = n.nUnits;
    else
        n.unit = "fields";
    end
    n.nFields = sum(used);
    
    % map field names back to mouse / session
    if isfield(cor, 'sessName')
        fn = arrayfun(@(c) firstString(c.sessName), cor(used));
        key = unique(tbl(:, {'mouse','sess','field'}), 'rows');
        keyName = string(key.mouse) + "_" + string(key.sess) + "_" + string(key.field);
        [found, loc] = ismember(fn, keyName);
        if ~all(found)
            warning('countCorrSampleSize: %d field(s) not found in tbl; mouse/session counts exclude them.', sum(~found));
        end
        key = key(loc(found), :);
        n.nMice = numel(unique(key.mouse));
        n.nSessions = height(unique(key(:, {'mouse','sess'}), 'rows'));
    else % older saves without sessName: one struct element per session/field
        n.nSessions = numel(cor);
        n.nMice = NaN;
    end
end

function s = firstString(x)
s = string(x);
s = s(1);
end
