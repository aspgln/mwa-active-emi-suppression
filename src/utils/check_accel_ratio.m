function [R_accel, fill_line] = check_accel_ratio(kdata, dim)
%CHECK_ACCEL_RATIO  Determine acceleration ratio & first fill line in k-space
%   [R_accel, fill_line] = check_accel_ratio(kdata)
%   [R_accel, fill_line] = check_accel_ratio(kdata, dim)
%
%   kdata – N-D k-space array (e.g., [Nfe, Nch, Npe, ...])
%   dim   – (optional) the dimension along which to check zeros
%           (default = 3)
%
%   R_accel   – number of leading zero lines before first non-zero
%   fill_line – which line to begin filling: 1 or 2

    % default dimension
    if nargin < 2 || isempty(dim)
        dim = 3;
    end

    % sanity check
    if dim < 1 || dim > ndims(kdata)
        error('check_accel_ratio:InvalidDim', ...
              'dim must be between 1 and %d.', ndims(kdata));
    end

    % build indexing: 1 for every dim, ':' at the chosen one
    N = ndims(kdata);
    subs = repmat({1}, 1, N);
    subs{dim} = ':';

    % extract the 1-D patch
    patch = kdata(subs{:});

    % initialize outputs
    R_accel  = 1;
    fill_line = 1;

    % check first two lines along that dim
    if patch(1) == 0
        % first line zero → find first non-zero
        while patch(R_accel) == 0
            R_accel = R_accel + 1;
        end
        fill_line = 2;
    elseif numel(patch) >= 2 && patch(2) == 0
        % second line zero → find first non-zero starting from 2
        while patch(1 + R_accel) == 0
            R_accel = R_accel + 1;
        end
        fill_line = 1;
    end
end
