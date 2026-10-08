function [q1, q2, q3, is_valid, reason] = scale_facs_to_q(lambda, edgeLen)
%SCALE_FACS_TO_Q Convert one unit's scale factors to target triangle points.
%   [q1, q2, q3, is_valid, reason] = scale_facs_to_q(lambda, edgeLen)
%
% Inputs
%   lambda  : 1x3 edge scale factors [lambda12, lambda23, lambda31]
%   edgeLen : reference unit edge length
%
% Outputs
%   q1, q2, q3 : 1x3 coordinates of target unit triangle vertices
%   is_valid   : true if conversion is valid (finite + triangle inequality)
%   reason     : '' when valid, otherwise failure reason string

q1 = [NaN, NaN, NaN];
q2 = [NaN, NaN, NaN];
q3 = [NaN, NaN, NaN];
is_valid = false;
reason = '';

if ~isnumeric(lambda) || numel(lambda) ~= 3 || any(~isfinite(lambda))
    reason = 'nonfinite_scale';
    return;
end

lambda = reshape(lambda, 1, 3);
e12 = lambda(1) * edgeLen;
e23 = lambda(2) * edgeLen;
e31 = lambda(3) * edgeLen;

if any([e12, e23, e31] <= 0)
    reason = 'nonpositive_edge';
    return;
end

% Triangle inequality
if (e12 + e23 <= e31) || (e23 + e31 <= e12) || (e31 + e12 <= e23)
    reason = 'invalid_target_triangle';
    return;
end

% Build target unit coordinates
q1 = [0, 0, 0];
q2 = [0, -e12, 0];
q3_y = (e23^2 - e31^2 - e12^2) / (2 * e12);
q3_x = -sqrt(max(e31^2 - q3_y^2, 0));
q3 = [q3_x, q3_y, 0];

is_valid = true;
end
