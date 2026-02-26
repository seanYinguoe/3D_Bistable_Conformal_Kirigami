function [anisotropy_level, lambda_sorted, sort_idx] = scale_facs_to_angles(scale_facs)
%SCALE_FACS_TO_ANGLES Convert triangle edge scale factors to anisotropy list.
%   [anisotropy_level, lambda_sorted, sort_idx] = scale_facs_to_angles(scale_facs)
%
% Input:
%   scale_facs : Mx3 edge scale factors per triangle
%                [edge(n1,n2), edge(n2,n3), edge(n3,n1)]
%
% Output:
%   anisotropy_level : Mx4 list [alpha1 alpha2 alpha3 eps_bist], where
%                      alpha1 >= alpha2 >= alpha3 and eps_bist = lambda3
%   lambda_sorted: Mx3 sorted scale factors [lambda1 lambda2 lambda3],
%                  lambda1 >= lambda2 >= lambda3
%   sort_idx     : Mx3 indices mapping sorted lambdas to original columns
%
% Notes:
%   - Uses triangle geometry (law of cosines) on side lengths proportional to
%     lambda values. This is equivalent to law-of-sines consistency used in
%     your workflow.
%   - Rows that violate triangle inequality or contain invalid values return NaN.

if size(scale_facs,2) ~= 3
    error('scale_facs must be Mx3.');
end

% Sort lambdas descending per triangle: lambda1 >= lambda2 >= lambda3
[lambda_sorted, sort_idx] = sort(scale_facs, 2, 'descend');

l1 = lambda_sorted(:,1);
l2 = lambda_sorted(:,2);
l3 = lambda_sorted(:,3);

% Initialize outputs
M = size(scale_facs,1);
anisotropy_level = nan(M,4, 'like', scale_facs);

% Valid rows: finite, positive, and satisfy triangle inequality
valid = all(isfinite(lambda_sorted),2) & all(lambda_sorted > 0,2) & ...
        (l1 < l2 + l3) & (l2 < l1 + l3) & (l3 < l1 + l2);

if any(valid)
    lv1 = l1(valid); lv2 = l2(valid); lv3 = l3(valid);

    % Law of cosines (angles opposite corresponding side lengths)
    c1 = (lv2.^2 + lv3.^2 - lv1.^2) ./ (2*lv2.*lv3);
    c2 = (lv1.^2 + lv3.^2 - lv2.^2) ./ (2*lv1.*lv3);

    % Numeric safety
    c1 = min(max(c1, -1), 1);
    c2 = min(max(c2, -1), 1);

    a1 = acos(c1);
    a2 = acos(c2);
    a3 = pi - a1 - a2;

    % eps_bist requested as lambda3
    eps_bist = lv3;
    anisotropy_level(valid,:) = [a1, a2, a3, eps_bist];
end

%{
% Example:
% sf = [1.60 1.45 1.20; 1.20 1.20 1.20];
% [anisotropy_level, lambda_sorted] = scale_facs_to_angles(sf)
%}

end
