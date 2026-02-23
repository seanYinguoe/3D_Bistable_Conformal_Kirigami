%% Sensitivity study over internal angles (alpha1, alpha2, alpha3) and beta
% alpha3 in [pi/6, pi/3]
% Enforce alpha1 + alpha2 + alpha3 = pi, and alpha1 >= alpha2 >= alpha3.

clear; clc;

%% ---- Angle grid (independent of beta) ----
interval = pi/36;

alpha3_vec = pi/6 : interval : pi/3;      % α3 ∈ [π/8, π/3]
alpha2_vec = pi/6 : interval : 5*pi/12;   % α2 ∈ [π/8, 7π/16]

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
        q1, q2, q3,0); % 1: plot energy curve; 0: not plotting
    eps_all(idx,1) = eps_bist;
    eta_all(idx,1) = eta_val;
end

%% Build final table: anisotropy_study
anisotropy_study = table( ...
    a1_all, a2_all, a3_all, eps_all, eta_all, beta_all, ...
    'VariableNames', {'a1','a2','a3','eps_bist','eta_val','beta'});

save ternary_study_beta_test2

%% Plot ternary figure
% % Keep rows with eta_val > threshold and finite eps_bist.
% etaThreshold = 0.05;
% isValid = isfinite(anisotropy_study.eta_val) & ...
%           (anisotropy_study.eta_val > etaThreshold) & ...
%           isfinite(anisotropy_study.eps_bist);
% T_valid = anisotropy_study(isValid, :);
% 
% % Select beta slice for plotting.
% T_filter = T_valid(abs(T_valid.beta - 0) < 1e-4, :);
% 
% plot_ternary(T_valid.a1, T_valid.a2, T_valid.a3, T_valid.eps_bist, 'scatter');
% 
% % Remove isolated outliers in ternary coordinates (toolbox-free).
% % outlier_strength: larger means more aggressive filtering.
% outlier_strength = 3;
% k_nn = 6;
% T_filter = filter_ternary_outliers_knn(T_filter, k_nn, outlier_strength);
% 
% plot_ternary(T_filter.a1, T_filter.a2, T_filter.a3, T_filter.eps_bist, 'scatter');
% 
% T_plot = interpolate_ternary_smooth(T_filter, 100, 1.3, etaThreshold);
% 
% plot_ternary(T_plot.a1, T_plot.a2, T_plot.a3, T_plot.eps_bist, 'scatter');
% 
% if ~isempty(T_plot)
%     eps_ref = T_filter.eps_bist(isfinite(T_filter.eps_bist));
%     clim_eps = [min(eps_ref) max(eps_ref)];
%     plot_ternary(T_plot.a1, T_plot.a2, T_plot.a3, T_plot.eps_bist, 'field', clim_eps);
% 
%     eta_ref = T_filter.eta_val(isfinite(T_filter.eta_val));
%     clim_eta = [min(eta_ref) max(eta_ref)];
%     plot_ternary(T_plot.a1, T_plot.a2, T_plot.a3, T_plot.eta_val, 'field', clim_eta);
% end

%% Plot 3D ternary figure in term of beta
%plot_ternary_3D(anisotropy_study.a1, anisotropy_study.a2, anisotropy_study.a3, anisotropy_study.eps_bist, anisotropy_study.beta)
