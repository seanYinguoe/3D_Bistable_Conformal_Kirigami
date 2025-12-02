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
beta_vec = linspace(0, pi/15, 50);
%beta_vec = pi/40;
nBeta = numel(beta_vec);

%% ---- Geometry & scale ----
edgeLen = 15;
scale   = 1.63;

%% ---- Storage for all (a1,a2,a3,eps_bist,eta_val,beta) ----
a1_all   = [];
a2_all   = [];
a3_all   = [];
eps_all  = [];
eta_all  = [];
beta_all = [];

%% ================== MAIN LOOPS ==================
for ib = 1:nBeta
    beta = beta_vec(ib);

    % Geometry parameters for this beta
    p_geom = struct();
    p_geom.l1   = edgeLen * 0.85;
    p_geom.l4   = edgeLen * 0.05;
    p_geom.t    = edgeLen * 0.015;
    p_geom.beta = beta;

    for i = 1:nConfig
        a1 = alpha1_list(i);
        a2 = alpha2_list(i);
        a3 = alpha3_list(i);

        % ----- Stretch ratios from angles -----
        lam1_ratio = sin(a1)/sin(a3);
        lam2_ratio = sin(a2)/sin(a3);
        lam3_ratio = 1.0;

        lambda1 = lam1_ratio * scale;
        lambda2 = lam2_ratio * scale;
        lambda3 = lam3_ratio * scale;

        % ----- Deformed edge lengths -----
        L1 = lambda1 * edgeLen;
        L2 = lambda2 * edgeLen;
        L3 = lambda3 * edgeLen;

        % ----- Triangle coordinates -----
        q1 = [0, 0];
        q2 = [0, -L3];
        q3_y = (L1^2 - L2^2 - L3^2) / (2 * L3);
        inside = L2^2 - q3_y^2;

        if inside <= 0
            eps_bist = NaN;
            eta_val  = NaN;
        else
            q3_x = -sqrt(inside);
            q3   = [q3_x, q3_y];

            % ----- Call bistability analysis -----
            try
                [eps_bist, eta_val] = bistability_analysis( ...
                    p_geom.l1, p_geom.l4, p_geom.beta, p_geom.t, edgeLen, ...
                    q1, q2, q3);
            catch ME   % ME = MException 异常对象
                fprintf('alpha3 = %.6f, alpha2 = %.6f\n,, alpha1 = %.6f\n', a3, a2, a1);
                disp('Inputs to bistability_analysis:');  % 显示传入的参数
                % 比如：
                % disp(p_geom);
                % disp(stretch_vec);

                fprintf('Error message: %s\n', ME.message);

                keyboard;  % 进入调试模式，停在这里
            end

            if isnan(eps_bist) || eta_val < 0.1
                eps_bist = NaN;
                eta_val  = NaN;
            end
        end

        % store results
        a1_all(end+1,1)   = a1;
        a2_all(end+1,1)   = a2;
        a3_all(end+1,1)   = a3;
        eps_all(end+1,1)  = eps_bist;
        eta_all(end+1,1)  = eta_val;
        beta_all(end+1,1) = beta;
    end
end

%% Build final table: anisotropy_study
anisotropy_study = table( ...
    a1_all, a2_all, a3_all, eps_all, eta_all, beta_all, ...
    'VariableNames', {'a1','a2','a3','eps_bist','eta_val','beta'});

save ternary_study_beta

%plot_ternary(a1_all, a2_all, a3_all, eps_all)
%plot_ternary(alpha1_list, alpha2_list, alpha3_list, epsb_list);
% T0 = anisotropy_study(abs(anisotropy_study.beta - 0.005416539057913) < 1e-4, :);
% 
% if T0.a3 < 1.33 && T0.a2 < 1.33
%     T0.eps_bist = NaN;
%     T0.eta_val = NaN;
% end
% 
% plot_ternary( ...
%     T0.a1, ...            % alpha1
%     T0.a2, ...            % alpha2
%     T0.a3, ...            % alpha3
%     T0.eps_bist );        % index field (e.g. eps_bist)



% a1 = anisotropy_study(67,:).a1;
% a2 = anisotropy_study(67,:).a2;
% a3 = anisotropy_study(67,:).a3;

