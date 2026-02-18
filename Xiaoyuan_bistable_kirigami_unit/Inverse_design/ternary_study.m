%% Sensitivity study over internal angles (alpha1, alpha2, alpha3) and beta
% alpha3 in [pi/6, pi/3]
% Enforce alpha1 + alpha2 + alpha3 = pi, and alpha1 >= alpha2 >= alpha3.

clear; clc;

%% ---- Angle grid (independent of beta) ----
nStep = 10;                         
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
beta_vec = linspace(0, pi/15, 5);
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
        q1, q2, q3);

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

% %% Plot ternary figure
T0 = anisotropy_study(abs(anisotropy_study.beta - 0.209439510239320) < 1e-4, :);
% 
% mask = (T0.a1 > pi/2) | (T0.a3 < 5*pi/24 & T0.a2 > pi/3);
% T0.eps_bist(mask) = NaN;
% T0.eta_val(mask)  = NaN;
% T0 = anisotropy_study;
plot_ternary( ...
    T0.a1, ...            % alpha1
    T0.a2, ...            % alpha2
    T0.a3, ...            % alpha3
    T0.eps_bist );        % index field (e.g. eps_bist)

%% Plot 3D ternary figure in term of beta
%plot_ternary_3D(anisotropy_study.a1, anisotropy_study.a2, anisotropy_study.a3, anisotropy_study.eps_bist, anisotropy_study.beta)

%% Test
% a1 = anisotropy_study(67,:).a1;
% a2 = anisotropy_study(67,:).a2;
% a3 = anisotropy_study(67,:).a3;
