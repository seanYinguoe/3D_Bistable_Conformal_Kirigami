%% Sensitivity study over internal angles (alpha1, alpha2, alpha3) and beta
% alpha3 in [pi/6, pi/3]
% Enforce alpha1 + alpha2 + alpha3 = pi, and alpha1 >= alpha2 >= alpha3.

clear; clc;

%% ---- Angle grid (independent of beta) ----
nStep = 20;                         
alpha3_vec = linspace(pi/6, pi/3, nStep);       % α3 ∈ [π/6, π/3]
alpha2_vec = linspace(pi/6, 5*pi/12, nStep);   % α2 ∈ [π/6, 5π/12]

[alpha2_grid, alpha3_grid] = meshgrid(alpha2_vec, alpha3_vec);
alpha1_grid = pi - alpha2_grid - alpha3_grid;

tol = 1e-8;

mask = (alpha1_grid > -tol) & ...
       (alpha1_grid >= alpha2_grid - tol) & ...
       (alpha2_grid >= alpha3_grid - tol);

alpha1_list = alpha1_grid(mask);
alpha2_list = alpha2_grid(mask);
alpha3_list = alpha3_grid(mask);

nConfig = numel(alpha1_list);

%% ---- Beta sweep settings ----
beta_vec = linspace(0, pi/15, 20);
%beta_vec = pi/40;
nBeta = numel(beta_vec);

%% ---- Geometry & scale ----
edgeLen = 15;
scale   = 1.63;

%% ---- Storage for all (a1,a2,a3,eps_bist,eta_val,beta) ----
% OPT: Preallocate outputs to avoid dynamic growth
nTotal   = nBeta * nConfig;
a1_all   = repmat(alpha1_list, nBeta, 1);
a2_all   = repmat(alpha2_list, nBeta, 1);
a3_all   = repmat(alpha3_list, nBeta, 1);
beta_all = kron(beta_vec(:), ones(nConfig, 1));
eps_all  = NaN(nTotal, 1);
eta_all  = NaN(nTotal, 1);

% OPT: Precompute geometry that does not depend on beta
sin_a1 = sin(alpha1_list);
sin_a2 = sin(alpha2_list);
sin_a3 = sin(alpha3_list);

lam1_ratio = sin_a1 ./ sin_a3;
lam2_ratio = sin_a2 ./ sin_a3;
lam3_ratio = ones(nConfig, 1);

lambda1 = lam1_ratio * scale;
lambda2 = lam2_ratio * scale;
lambda3 = lam3_ratio * scale;

L1 = lambda1 * edgeLen;
L2 = lambda2 * edgeLen;
L3 = lambda3 * edgeLen;

q3_y = (L1.^2 - L2.^2 - L3.^2) ./ (2 * L3);
inside = L2.^2 - q3_y.^2;
valid  = inside > 0;
q3_x   = NaN(nConfig, 1);
q3_x(valid) = -sqrt(inside(valid));

%% ================== MAIN LOOPS ==================
% Geometry parameters independent of (a1,a2,a3)
l1_geom = edgeLen * 0.85;
l4_geom = edgeLen * 0.05;
t_geom  = edgeLen * 0.015;

% PARFOR over flattened (beta, config) index
parfor idx = 1:nTotal
    i  = mod(idx - 1, nConfig) + 1;
    ib = floor((idx - 1) / nConfig) + 1;

    % Skip invalid triangles early  
    if ~valid(i)
        continue;
    end

    beta = beta_vec(ib);

    % ----- Triangle coordinates -----
    q1 = [0, 0];
    q2 = [0, -L3(i)];
    q3 = [q3_x(i), q3_y(i)];

    [eps_bist, eta_val] = bistability_analysis( ...
        l1_geom, l4_geom, beta, t_geom, edgeLen, ...
        q1, q2, q3,0);

    if ~isnan(eps_bist) && eta_val >= 0.1
        eps_all(idx, 1) = eps_bist;
        eta_all(idx, 1) = eta_val;
    end
end

%% Build final table: anisotropy_study
anisotropy_study = table( ...
    a1_all, a2_all, a3_all, eps_all, eta_all, beta_all, ...
    'VariableNames', {'a1','a2','a3','eps_bist','eta_val','beta'});

save ternary_study_beta_test1

%% Plot ternary figure
% Keep rows with eta_val > threshold and finite eps_bist.
etaThreshold = 0.15;
isValid = isfinite(anisotropy_study.eta_val) & ...
          (anisotropy_study.eta_val > etaThreshold) & ...
          isfinite(anisotropy_study.eps_bist);
T_valid = anisotropy_study(isValid, :);

% Select beta slice for plotting.
T_plot = T_valid(abs(T_valid.beta - 0) < 1e-4, :);

% Remove isolated outliers in ternary coordinates (toolbox-free).
% outlier_strength: larger means more aggressive filtering.
outlier_strength = 3;
k_nn = 6;
n_pts = height(T_plot);
if n_pts > 3
    % Normalize ternary coordinates and map to 2D Cartesian.
    s = T_plot.a1 + T_plot.a2 + T_plot.a3;
    B = T_plot.a2 ./ s;
    C = T_plot.a3 ./ s;
    x = B + 0.5 * C;
    y = (sqrt(3)/2) * C;
    X = [x, y];

    % Pairwise distances and kNN mean distances.
    G = sum(X.^2, 2);
    D2 = max(G + G' - 2*(X*X'), 0);
    D = sqrt(D2);
    D(1:n_pts+1:end) = inf;

    k_eff = min(k_nn, n_pts - 1);
    if k_eff >= 1
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
        T_plot = T_plot(~is_outlier, :);
    end
end

plot_ternary( ...
    T_plot.a1, ...            % alpha1
    T_plot.a2, ...            % alpha2
    T_plot.a3, ...            % alpha3
    T_plot.eps_bist );        % index field (eps_bist)

%% Plot 3D ternary figure in term of beta
%plot_ternary_3D(anisotropy_study.a1, anisotropy_study.a2, anisotropy_study.a3, anisotropy_study.eps_bist, anisotropy_study.beta)

%% Test
% a1 = anisotropy_study(67,:).a1;
% a2 = anisotropy_study(67,:).a2;
% a3 = anisotropy_study(67,:).a3;
