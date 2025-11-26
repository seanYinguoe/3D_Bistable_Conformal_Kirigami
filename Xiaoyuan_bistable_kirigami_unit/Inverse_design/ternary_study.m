%% Sensitivity study over internal angles (alpha1, alpha2, alpha3)
% alpha3 in [pi/6, pi/3]
% Enforce alpha1 + alpha2 + alpha3 = pi, and alpha1 >= alpha2 >= alpha3.

clear; clc;

nStep = 10;                         
alpha3_vec = linspace(pi/6, pi/3, nStep);   % α3 ∈ [π/6, π/3]
alpha2_vec = linspace(pi/6, 5*pi/12, nStep);   % α2 ∈ [π/6, π/2]

[alpha2_grid, alpha3_grid] = meshgrid(alpha2_vec, alpha3_vec);

alpha1_grid = pi - alpha2_grid - alpha3_grid;

mask = (alpha1_grid > 0) & ...
       (alpha1_grid >= alpha2_grid) & ...
       (alpha2_grid >= alpha3_grid);

alpha1_list = alpha1_grid(mask);
alpha2_list = alpha2_grid(mask);
alpha3_list = alpha3_grid(mask);

% Storage arrays
eta_list    = [];
epsb_list   = [];

% Geometry parameters
edgeLen = 15;
p_geom = struct();
p_geom.l1   = edgeLen * 0.85;
p_geom.l4   = edgeLen * 0.05;
p_geom.t    = edgeLen * 0.015;
p_geom.beta = 0;

scale = 1.6;                             % overall stretch scale

%% Sweep alpha3 and alpha2
for i = 1:numel(alpha1_list)
    a1 = alpha1_list(i);
    a2 = alpha2_list(i);
    a3 = alpha3_list(i);
    % ----- Stretch ratios from angles -----
    lam1_ratio = sin(a1)/sin(a3);
    lam2_ratio = sin(a2)/sin(a3);
    lam3_ratio = 1.0;

    lambda1 = lam1_ratio * scale;
    lambda2 = lam2_ratio * scale;
    lambda3 = lam3_ratio * scale;    % = scale

    % ----- Corresponding deformed edge lengths -----
    L1 = lambda1 * edgeLen;   % edge 1 (between q2 and q3)
    L2 = lambda2 * edgeLen;   % edge 2 (between q1 and q3)
    L3 = lambda3 * edgeLen;   % edge 3 (between q2 and q1)

    % Build triangle coordinates
    q1 = [0, 0];
    q2 = [0, -L3];
    q3_y = (L1^2 - L2^2 - L3^2) / (2 * L3);
    inside = L2^2 - q3_y^2;
    if inside <= 0
        continue;
    end
    q3_x = -sqrt(inside);
    q3   = [q3_x, q3_y];

    % ----- Call bistability analysis -----
    [strain_bist, bistability] = bistability_analysis( ...
        p_geom.l1, p_geom.l4, p_geom.beta, p_geom.t, edgeLen, ...
        q1, q2, q3);

    % Store results (treat eta < 0.1 as monostable if you like)
    if ~isnan(strain_bist) && bistability >= 0.1
        eta_list(end+1,1)  = bistability;
        epsb_list(end+1,1) = strain_bist;
    else
        eta_list(end+1,1)  = NaN;
        epsb_list(end+1,1) = NaN;
    end
end

%% Ternary plot (uses your custom function)
plot_ternary(alpha1_list, alpha2_list, alpha3_list, epsb_list);

