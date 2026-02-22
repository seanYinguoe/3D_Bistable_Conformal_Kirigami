function T_out = filter_ternary_outliers_knn(T_in, k_nn, outlier_strength)
%FILTER_TERNARY_OUTLIERS_KNN Remove isolated ternary outliers via kNN distance.
%   T_out = filter_ternary_outliers_knn(T_in, k_nn, outlier_strength)

T_out = T_in;
n_pts = height(T_in);
if n_pts <= 3
    return;
end

% Normalize ternary coordinates and map to 2D Cartesian.
s = T_in.a1 + T_in.a2 + T_in.a3;
B = T_in.a2 ./ s;
C = T_in.a3 ./ s;
x = B + 0.5 * C;
y = (sqrt(3)/2) * C;
X = [x, y];

% Pairwise distances and kNN mean distances.
G = sum(X.^2, 2);
D2 = max(G + G' - 2*(X*X'), 0);
D = sqrt(D2);
D(1:n_pts+1:end) = inf;

k_eff = min(k_nn, n_pts - 1);
if k_eff < 1
    return;
end

D_sorted = sort(D, 2, 'ascend');
knn_mean_dist = mean(D_sorted(:, 1:k_eff), 2);

% Robust threshold: median + outlier_strength * 1.4826*MAD
d_med = median(knn_mean_dist);
d_mad = median(abs(knn_mean_dist - d_med));
if d_mad <= eps
    thr_robust = inf;
else
    thr_robust = d_med + outlier_strength * 1.4826 * d_mad;
end

% Secondary threshold: 99th percentile.
d_sorted = sort(knn_mean_dist);
idx99 = max(1, ceil(0.99 * numel(d_sorted)));
thr_p99 = d_sorted(idx99);

is_outlier = (knn_mean_dist > thr_robust) | (knn_mean_dist > thr_p99);
T_out = T_in(~is_outlier, :);
end
