function plot_error(cur_axis, x,y,err,varargin)

    % use the fill function to plot erorr cloud
    % credit to chatgpt 
    p = inputParser;
    addParameter(p,'Color',[0 0 1])
    addParameter(p,'Alpha',0.3)
    addParameter(p,'LineWidth',2)
    addParameter(p,'LineStyle','-')
    addParameter(p,'fadeFactor', 0.5) % 0 = no change, 1 = fully white
    parse(p, varargin{:});
    c  = p.Results.Color;
    a  = p.Results.Alpha;
    lw = p.Results.LineWidth;
    ls = p.Results.LineStyle;
    fadeFactor = p.Results.fadeFactor;
    
    % -------------------- Ensure column vectors --------------------
    if isvector(y),   y   = y(:);   end
    if isvector(err), err = err(:); end

    [Ny, Ky] = size(y);
    [Ne, Ke] = size(err);

    % Allow err to be N x 1 and y to be N x K (replicate err across columns)
    if Ne ~= Ny
        error('plot_error:SizeMismatch', 'y and err must have the same number of rows.');
    end
    if Ke == 1 && Ky > 1
        err = repmat(err, 1, Ky);
        Ke = Ky;
    end
    if Ky ~= Ke
        error('plot_error:SizeMismatch', ...
            'y and err must have the same number of columns (each column is one trace).');
    end

    K = Ky;

    % -------------------- Validate / reshape x --------------------
    if isvector(x)
        x = x(:);  % N x 1
        if numel(x) ~= Ny
            error('plot_error:SizeMismatch', 'Length of x must match number of rows in y.');
        end
        x = repmat(x, 1, K); % make N x K for uniform handling
    else
        if size(x,1) ~= Ny
            error('plot_error:SizeMismatch', 'x must have the same number of rows as y.');
        end
        if size(x,2) == 1 && K > 1
            x = repmat(x, 1, K);
        elseif size(x,2) ~= K
            error('plot_error:SizeMismatch', 'x must be N x 1 or N x K to match y.');
        end
    end

    % -------------------- Handle colors --------------------
    if size(c,1) == 1 & K == 1
        c = repmat(c, K, 1); 
    elseif size(c,1) == 1 & K ~= 1
        c = lines(K);
    elseif size(c,1) ~= K || size(c,2) ~= 3
        error('plot_error:ColorSize', ...
            'Color must be 1x3 or Kx3, where K is number of traces (columns).');
    end
    patchColors = c + (1 - c) * fadeFactor;
    
    % -------------------- Plot --------------------
    axes(cur_axis);
    hold(cur_axis, 'on');

    h = struct('patch', cell(1,K), 'line', cell(1,K));

    for k = 1:K
        xk = x(:,k);
        yk = y(:,k);
        ek = err(:,k);

        upper = yk + ek;
        lower = yk - ek;
        h(k).patch = fill(cur_axis, ...
            [xk; flipud(xk)], ...
            [upper; flipud(lower)], ...
            patchColors(k,:), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', a);

        h(k).line = plot(cur_axis, xk, yk, ...
            'Color', c(k,:), ...
            'LineWidth', lw, ...
            'LineStyle', ls);
    end
end