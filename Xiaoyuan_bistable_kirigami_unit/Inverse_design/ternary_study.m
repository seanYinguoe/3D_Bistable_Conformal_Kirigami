%% Sensitivity study over internal angles (alpha1, alpha2, alpha3)
clear; clc;

dA = pi/24; 
alpha3_vec = pi/6 : dA : pi/3;

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

%% Sweep alpha3 and alpha2
for i3 = 1:numel(alpha3_vec)
    
    a3 = alpha3_vec(i3);
    
    % For a given a3:
    % alpha2 ∈ [a3, (pi - a3)/2]
    a2_min = a3;
    a2_max = (pi - a3) / 2;
    
    a2_vec = a2_min : dA : a2_max;
    
    for a2 = a2_vec
        
        a1 = pi - a2 - a3;      % triangle sum
        
        % enforce ordering
        if ~(a1 >= a2 && a2 >= a3 && a1 > 0)
            continue;
        end
        
        % Compute stretch ratios
        lam1_ratio = sin(a1)/sin(a3);
        lam2_ratio = sin(a2)/sin(a3);
        lam3_ratio = 1.0;

        scale = 1.7;

        lambda1 = lam1_ratio * scale;
        lambda2 = lam2_ratio * scale;
        lambda3 = lam3_ratio * scale; 

        % Corresponding edge lengths
        L1 = lambda1 * edgeLen;   % edge 1 (between q2 and q3)
        L2 = lambda2 * edgeLen;   % edge 2 (between q1 and q3)
        L3 = lambda3 * edgeLen;   % edge 3 (between q2 and q1)

        q1 = [0, 0];
        q2 = [0, -L3];
        q3_y = (L1^2 - L2^2 - L3^2) / (2 * L3);
        q3_x = -sqrt(max(L2^2 - q3_y^2, 0));
        q3 = [q3_x, q3_y];

        % Call bistability analysis
        [strain_bist, bistability] = bistability_analysis( ...
            p_geom.l1, p_geom.l4, p_geom.beta, p_geom.t, edgeLen, ...
            q1, q2, q3);

        % Store angles
        alpha1_list(end+1,1) = a1;
        alpha2_list(end+1,1) = a2;
        alpha3_list(end+1,1) = a3;
        
        % Store results (use NaN for monostable)
        if ~isnan(strain_bist) && bistability > 0
            eta_list(end+1,1)  = bistability;
            epsb_list(end+1,1) = strain_bist;
        else
            eta_list(end+1,1)  = NaN;
            epsb_list(end+1,1) = NaN;
        end
    end
end

plot_ternary(alpha1_list, alpha2_list, alpha3_list, epsb_list, eta_list)