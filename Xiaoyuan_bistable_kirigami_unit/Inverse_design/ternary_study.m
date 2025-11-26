%% Sensitivity study over internal angles (alpha1, alpha2, alpha3)
% alpha3 in [pi/6, pi/3]
% For each alpha3, alpha2 in [alpha3, (pi - alpha3)/2]
% Enforce alpha1 + alpha2 + alpha3 = pi, and alpha1 >= alpha2 >= alpha3.

clear; clc;

nStep = 100;                            % number of samples per direction
alpha3_vec = linspace(pi/6, pi/3, nStep);

% Storage arrays
alpha1_list = [];
alpha2_list = [];
alpha3_list = [];
eta_list    = [];
epsb_list   = [];

% Geometry parameters
edgeLen = 15;
p_geom = struct();
p_geom.l1   = edgeLen * 0.85;
p_geom.l4   = edgeLen * 0.05;
p_geom.t    = edgeLen * 0.015;
p_geom.beta = 0;

scale = 1.7;                             % overall stretch scale

%% Sweep alpha3 and alpha2
for i3 = 1:numel(alpha3_vec)
    
    a3 = alpha3_vec(i3);
    
    % For a given a3:
    % alpha2 ∈ [a3, (pi - a3)/2]
    a2_min = a3;
    a2_max = (pi - a3) / 2;
    
    % Use linspace so endpoints are always included
    a2_vec = linspace(a2_min, a2_max, nStep);
    
    for a2 = a2_vec
        
        a1 = pi - a2 - a3;      % triangle sum
        
        % Enforce ordering alpha1 >= alpha2 >= alpha3
        if ~(a1 >= a2 && a2 >= a3 && a1 > 0)
            continue;
        end
        
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
            % degenerate / numerical issue, skip this configuration
            continue;
        end
        q3_x = -sqrt(inside);
        q3   = [q3_x, q3_y];

        % ----- Call bistability analysis -----
        [strain_bist, bistability] = bistability_analysis( ...
            p_geom.l1, p_geom.l4, p_geom.beta, p_geom.t, edgeLen, ...
            q1, q2, q3);

        % Store angles
        alpha1_list(end+1,1) = a1;
        alpha2_list(end+1,1) = a2;
        alpha3_list(end+1,1) = a3;
        
        % Store results (treat eta < 0.1 as monostable if you like)
        if ~isnan(strain_bist) && bistability >= 0.1
            eta_list(end+1,1)  = bistability;
            epsb_list(end+1,1) = strain_bist;
        else
            eta_list(end+1,1)  = NaN;
            epsb_list(end+1,1) = NaN;
        end
    end
end

plot_ternary(alpha1_list, alpha2_list, alpha3_list, epsb_list, eta_list)