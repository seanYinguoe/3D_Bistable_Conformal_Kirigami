%% Sensitivity study over internal angles (alpha1, alpha2, alpha3) and beta
% alpha3 in [pi/6, pi/3]
% Enforce alpha1 + alpha2 + alpha3 = pi, and alpha1 >= alpha2 >= alpha3.

clear; clc;

%% ---- Angle grid (independent of beta) ----
interval = pi/36;

alpha3_vec = pi/6 : interval : pi/3;      % α3 ∈ [π/6, π/3]
alpha2_vec = pi/6 : interval : 5*pi/12;   % α2 ∈ [π/6, 5π/12]

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
nTotal   = nBeta * nConfig;
a1_all   = repmat(alpha1_list, nBeta, 1);
a2_all   = repmat(alpha2_list, nBeta, 1);
a3_all   = repmat(alpha3_list, nBeta, 1);
beta_all = kron(beta_vec(:), ones(nConfig, 1));
eps_all  = NaN(nTotal, 1);
eta_all  = NaN(nTotal, 1);

% Precompute geometry that does not depend on beta
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

% Batch by beta (checkpoint each beta slice)
outDir = 'beta_batches';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for ib = 1:nBeta
    beta = beta_vec(ib);
    batchFile = fullfile(outDir, sprintf('beta_%03d.mat', ib));

    if exist(batchFile, 'file')
        S = load(batchFile, 'eps_beta', 'eta_beta');
        eps_beta = S.eps_beta;
        eta_beta = S.eta_beta;
        fprintf('Loaded beta %d/%d from checkpoint.\n', ib, nBeta);
    else
        eps_beta = NaN(nConfig,1);
        eta_beta = NaN(nConfig,1);

        parfor i = 1:nConfig
            if ~valid(i)
                continue;
            end

            q1 = [0, 0];
            q2 = [0, -L3(i)];
            q3 = [q3_x(i), q3_y(i)];

            [eps_bist, eta_val] = bistability_analysis( ...
                l1_geom, l4_geom, beta, t_geom, edgeLen, ...
                q1, q2, q3, 0); % 1: plot energy curve; 0: not plotting

            eps_beta(i,1) = eps_bist;
            eta_beta(i,1) = eta_val;
        end

        save(batchFile, 'ib', 'beta', 'eps_beta', 'eta_beta');
        fprintf('Saved beta %d/%d checkpoint.\n', ib, nBeta);
    end

    idx0 = (ib - 1) * nConfig + 1;
    idx1 = ib * nConfig;
    eps_all(idx0:idx1,1) = eps_beta;
    eta_all(idx0:idx1,1) = eta_beta;
end

%% Build final table: anisotropy_study
anisotropy_study = table( ...
    a1_all, a2_all, a3_all, eps_all, eta_all, beta_all, ...
    'VariableNames', {'a1','a2','a3','eps_bist','eta_val','beta'});

% Save only final table
save('anisotropy_study.mat', 'anisotropy_study');

%% Clean islands first, then fill missing, then plot ternary figure (optional)
[anisotropy_clean, reportClean] = clean_isolated_outliers(anisotropy_study);
[anisotropy_filled, reportFill] = clean_and_fill_anisotropy(anisotropy_clean);

do_post_plot = true;   % set true to plot
etaThreshold = 0.03;
beta_plot = 0.209439510239320;
beta_tol = 1e-4;

if do_post_plot
    T0 = anisotropy_filled(abs(anisotropy_filled.beta - beta_plot) < beta_tol, :);
    isB = isfinite(T0.eps_bist) & isfinite(T0.eta_val) & (T0.eta_val > etaThreshold);
    T_bist = T0(isB,:);
    if ~isempty(T_bist)
        plot_ternary(T_bist.a1, T_bist.a2, T_bist.a3, T_bist.eps_bist, 'scatter');
        plot_ternary(T_bist.a1, T_bist.a2, T_bist.a3, T_bist.eta_val, 'scatter');
    end
end

%% Plot 3D ternary figure in terms of beta
% plot_ternary_3D(anisotropy_study.a1, anisotropy_study.a2, anisotropy_study.a3, anisotropy_study.eps_bist, anisotropy_study.beta)
