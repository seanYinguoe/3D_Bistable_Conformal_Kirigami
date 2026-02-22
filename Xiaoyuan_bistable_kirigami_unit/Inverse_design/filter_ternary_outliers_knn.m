function T_out = filter_ternary_outliers_knn(T_in, k_nn, outlier_strength)
%FILTER_TERNARY_OUTLIERS_KNN Remove isolated ternary islands while preserving dense core.
%   T_out = filter_ternary_outliers_knn(T_in, k_nn, outlier_strength)
%
% Recommended defaults for regular ternary lattices:
%   k_nn ~ 6 to 10, outlier_strength ~ 2 to 3.
% These values suppress sparse islands while keeping the central dense cluster.

if nargin < 2 || isempty(k_nn)
    k_nn = 6;
end
if nargin < 3 || isempty(outlier_strength)
    outlier_strength = 3;
end

T_out = T_in;
if isempty(T_in)
    return;
end

requiredVars = {'a1','a2','a3'};
if ~all(ismember(requiredVars, T_in.Properties.VariableNames))
    error('Table must contain variables: a1, a2, a3.');
end

n_pts = height(T_in);
if n_pts <= 3
    return;
end

% Ternary -> Cartesian
[x, y] = ternary_to_xy(T_in.a1, T_in.a2, T_in.a3);
X = [x, y];

% Pairwise Euclidean distances (toolbox-free)
G = sum(X.^2, 2);
D2 = max(G + G' - 2*(X*X'), 0);
D = sqrt(D2);
D(1:n_pts+1:end) = inf;

% kNN mean distance per point
k_eff = min(max(1, round(k_nn)), n_pts - 1);
D_sorted = sort(D, 2, 'ascend');
knn_mean_dist = mean(D_sorted(:, 1:k_eff), 2);

% Robust base radius from local density
d_med = median(knn_mean_dist);
d_mad = median(abs(knn_mean_dist - d_med));
if d_mad <= eps
    r_base = d_med;
else
    r_base = d_med + outlier_strength * 1.4826 * d_mad;
end

% Clamp with q75 lower bound so dense core doesn't fragment
q75 = local_quantile(knn_mean_dist, 0.75);
r = max(r_base, q75);
if ~isfinite(r) || r <= 0
    r = max(q75, max(knn_mean_dist));
end
if ~isfinite(r) || r <= 0
    r = 1e-12;
end

% Build largest connected component with adaptive relaxation
min_frac = 0.30;
bestMask = false(n_pts,1);
for attempt = 1:4
    A = (D <= r);
    A(1:n_pts+1:end) = false;

    compId = connected_components_bfs(A);
    compSizes = accumarray(compId, 1);
    [maxSize, maxComp] = max(compSizes);

    bestMask = (compId == maxComp);
    if (maxSize / n_pts) >= min_frac
        break;
    end
    r = r * 1.2;
end

% Optional center safeguard: always keep isotropic-near points
x0 = 0.5;
y0 = sqrt(3)/6;
if ~isfinite(d_med) || d_med <= 0
    d_med = median(knn_mean_dist(knn_mean_dist > 0));
    if isempty(d_med) || ~isfinite(d_med)
        d_med = 1e-12;
    end
end
r0 = 2 * d_med;
d0 = sqrt((x - x0).^2 + (y - y0).^2);
keepCenter = (d0 <= r0);

% If center points connect to largest CC within r, they are already kept.
% Merge anyway to guarantee center is never removed.
keepMask = bestMask | keepCenter;

% Fallback to original robust-distance removal if CC unexpectedly empty
if ~any(keepMask)
    if d_mad <= eps
        thr_robust = inf;
    else
        thr_robust = d_med + outlier_strength * 1.4826 * d_mad;
    end
    d_sorted = sort(knn_mean_dist);
    idx99 = max(1, ceil(0.99 * numel(d_sorted)));
    thr_p99 = d_sorted(idx99);
    is_outlier = (knn_mean_dist > thr_robust) | (knn_mean_dist > thr_p99);
    keepMask = (~is_outlier) | keepCenter;
end

T_out = T_in(keepMask, :);

end

function q = local_quantile(v, p)
% Toolbox-free quantile estimator using sorted index.
v = sort(v(:));
n = numel(v);
if n == 0
    q = NaN;
    return;
end
p = min(max(p, 0), 1);
idx = max(1, ceil(p * n));
q = v(idx);
end

function compId = connected_components_bfs(A)
% Connected components via BFS on logical adjacency matrix.
n = size(A,1);
compId = zeros(n,1);
comp = 0;
for i = 1:n
    if compId(i) ~= 0
        continue;
    end
    comp = comp + 1;
    compId(i) = comp;

    queue = zeros(n,1);
    qHead = 1;
    qTail = 1;
    queue(qTail) = i;

    while qHead <= qTail
        u = queue(qHead);
        qHead = qHead + 1;

        nbr = find(A(u,:));
        for k = 1:numel(nbr)
            v = nbr(k);
            if compId(v) == 0
                compId(v) = comp;
                qTail = qTail + 1;
                queue(qTail) = v;
            end
        end
    end
end
end

function [x, y] = ternary_to_xy(a1, a2, a3)
% Local ternary-to-Cartesian mapping to keep this file self-contained.
s = a1 + a2 + a3;
B = a2 ./ s;
C = a3 ./ s;
x = B + 0.5 * C;
y = (sqrt(3) / 2) * C;
end
