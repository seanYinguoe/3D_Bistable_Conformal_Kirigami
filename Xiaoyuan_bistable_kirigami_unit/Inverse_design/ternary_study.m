%% Sensitivity study over internal angles (alpha1, alpha2, alpha3) and beta
% alpha3 in [pi/6, pi/3]
% Enforce alpha1 + alpha2 + alpha3 = pi, and alpha1 >= alpha2 >= alpha3.

clear; clc;

%% ---- Output folder (avoid overwriting previous runs) ----
run_tag = ['ternary_refine_' datestr(now,'yyyymmdd_HHMMSS')];
out_root = fullfile('output', run_tag);
if ~exist(out_root, 'dir')
    mkdir(out_root);
end

%% ---- Angle grid (independent of beta) ----
interval = pi/144;
%interval = pi/36;
alpha3_vec = 10*pi/36 : interval : 12*pi/36;
alpha2_vec = 10*pi/36 : interval : 13*pi/36;
% alpha1 is always determined by angle sum:
% alpha1 = pi - alpha2 - alpha3

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
beta_vec = linspace(0, pi/20, 20);
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
l1_geom = edgeLen * 0.80;
l4_geom = edgeLen * 0.05;
t_geom  = edgeLen * 0.015;

% Parameter-specific save tag so repeated studies do not overwrite each other.
l1_tag = sprintf('l1_%03d', round(100 * l1_geom / edgeLen));
beta_tag = sprintf('beta_%03d_%03d', round(1000 * min(beta_vec)), round(1000 * max(beta_vec)));
study_tag = [l1_tag '_' beta_tag];

% Batch by beta (checkpoint each beta slice)
outDir = fullfile(out_root, 'beta_batches_refine');
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
anisotropy_study_refine = table( ...
    a1_all, a2_all, a3_all, eps_all, eta_all, beta_all, ...
    'VariableNames', {'a1','a2','a3','eps_bist','eta_val','beta'});

% Save only final table under a parameter-specific name.
study_var_name = ['anisotropy_study_refine_' study_tag];
S_study = struct();
S_study.(study_var_name) = anisotropy_study_refine;
save(fullfile(out_root, [study_var_name '.mat']), '-struct', 'S_study');

%% Clean islands first, then fill missing, then plot ternary figure (optional)
[anisotropy_clean, reportClean] = clean_isolated_outliers(anisotropy_study_refine);
[anisotropy_filled, reportFill] = clean_and_fill_anisotropy(anisotropy_clean);
mask = isfinite(anisotropy_filled.eta_val) & isfinite(anisotropy_filled.eps_bist);
anisotropy_filter = anisotropy_filled(mask,:);
filter_var_name = ['anisotropy_filter_' study_tag];
S_filter = struct();
S_filter.(filter_var_name) = anisotropy_filter;
save(fullfile(out_root, [filter_var_name '.mat']), '-struct', 'S_filter');

% Refine eps_bist toward endpoint-consistent deployment
opts_refine = struct();
opts_refine.edgeLen = edgeLen;
opts_refine.l1 = l1_geom;
opts_refine.l4 = l4_geom;
opts_refine.t = t_geom;
opts_refine.etaThreshold = 0.03;
opts_refine.alphaTarget = 0.99;
opts_refine.maxIter = 6;
opts_refine.epsFloor = 1e-4;
opts_refine.epsCeil = 2.0;
opts_refine.alphaTol = 1e-3;

anisotropy_filter_refine = refine_eps_bist_to_endpoint(anisotropy_filter, opts_refine);
filter_refine_var_name = ['anisotropy_filter_refine_' study_tag];
S_filter_refine = struct();
S_filter_refine.(filter_refine_var_name) = anisotropy_filter_refine;
S_filter_refine.opts_refine = opts_refine;
save(fullfile(out_root, [filter_refine_var_name '.mat']), '-struct', 'S_filter_refine');

% Save filled table
filled_var_name = ['anisotropy_filled_refine_' study_tag];
S_filled = struct();
S_filled.(filled_var_name) = anisotropy_filled;
S_filled.anisotropy_clean = anisotropy_clean;
S_filled.reportClean = reportClean;
S_filled.reportFill = reportFill;
save(fullfile(out_root, [filled_var_name '.mat']), '-struct', 'S_filled');

fprintf('Saved refined ternary outputs to: %s\n', out_root);

% anisotropy_filled = anisotropy_filled_pi12;
% 
% do_post_plot = true;   % set true to plot
% etaThreshold = 0.03;
% beta_plot = 0.220462642357178;
% beta_tol = 1e-4;
% 
% if do_post_plot
%     T0 = anisotropy_filled(abs(anisotropy_filled.beta - beta_plot) < beta_tol, :);
%     isB = isfinite(T0.eps_bist) & isfinite(T0.eta_val) & (T0.eta_val > etaThreshold);
%     T_bist = T0(isB,:);
%     if ~isempty(T_bist)
%         plot_ternary(T_bist.a1, T_bist.a2, T_bist.a3, T_bist.eps_bist, 'scatter');
%         plot_ternary(T_bist.a1, T_bist.a2, T_bist.a3, T_bist.eta_val, 'scatter');
%     end
% end

%% Plot 3D ternary figure in terms of beta (scatter version)
plot_ternary_3D(anisotropy_clean.a1, anisotropy_clean.a2, anisotropy_clean.a3, anisotropy_clean.eps_bist, anisotropy_clean.beta)

%% Plot bistable region contours stacked along beta axis
% opts.threshold : eps_bist cut-off that defines "bistable" (default 0)
% opts.nGrid     : interpolation grid resolution (default 100)
% opts.faceAlpha : filled-patch transparency (default 0.4)
opts_contour = struct('threshold', 0, 'nGrid', 100, 'faceAlpha', 0.4, 'cmap', 'parula');
plot_ternary_beta_contour(anisotropy_clean.a1, anisotropy_clean.a2, anisotropy_clean.a3, ...
    anisotropy_clean.eps_bist, anisotropy_clean.beta, opts_contour)
