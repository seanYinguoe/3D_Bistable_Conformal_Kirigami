function t = assign_thickness(bistability, t_min, t_max, thr)
%ASSIGN_THICKNESS Assign ligament thickness from bistability.
%   t = assign_thickness(bistability, t_min, t_max)
%   t = assign_thickness(bistability, t_min, t_max, thr)
%
% Rule:
%   bistability >= thr -> t_max
%   bistability <  thr -> t_min
%   bistability = NaN  -> t_min

if nargin < 4 || isempty(thr)
    thr = 0.4;
end

if ~isscalar(t_min) || ~isscalar(t_max) || ~isscalar(thr)
    error('t_min, t_max, and thr must be scalar values.');
end
if ~isfinite(t_min) || ~isfinite(t_max) || ~isfinite(thr)
    error('t_min, t_max, and thr must be finite.');
end
if t_min < 0 || t_max < 0
    error('Negative ligament thickness is not allowed.');
end
if t_min > t_max
    error('t_min must be less than or equal to t_max.');
end

if istable(bistability)
    if width(bistability) ~= 1
        error('Table input for bistability must have exactly one column.');
    end
    bistability = bistability{:,1};
end

if ~isnumeric(bistability)
    error('bistability must be a numeric vector or a single table column.');
end

t = t_min * ones(size(bistability));
mask = isfinite(bistability) & (bistability >= thr);
t(mask) = t_max;
end
