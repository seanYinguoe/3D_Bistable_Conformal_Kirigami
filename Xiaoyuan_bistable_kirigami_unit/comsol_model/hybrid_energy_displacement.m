function [strain, U_vec] = hybrid_energy_displacement(l1,l4,t,beta,edgeLen)
q1 = [0, 0, 0] * 1.56;
q3 = [0, -edgeLen, 0] * 1.56;
q2 = [-sqrt(3)/2 * edgeLen, -0.5 * edgeLen, 0] * 1.56;

%% Parameters
n_step = 1000; % Number of interpolation steps
strain = zeros(1,n_step);        % Physical strain
U_vec = nan(1, n_step);          % Total energy
energy_b_vec = nan(1, n_step);   % Bending energy
energy_s_vec = nan(1, n_step);   % Stretch energy

%% Original (undeformed) boundary nodes
p1_orig = [0, 0, 0];
p3_orig = [0, -edgeLen, 0];
p2_orig = [-sqrt(3)/2 * edgeLen, -0.5 * edgeLen, 0];

%% Deformed boundary nodes (from q1, q2, q3)
edge1_def = norm(q1 - q2);
edge2_def = norm(q2 - q3);
edge3_def = norm(q1 - q3);
p1_def = [0, 0, 0];
p3_def = [0, -edge3_def, 0];
p2_y = (edge2_def^2 - edge1_def^2 - edge3_def^2) / (2 * edge3_def);
p2_x = -sqrt(edge1_def^2 - p2_y^2);
p2_def = [p2_x, p2_y, 0];

%% Loop through interpolation steps
for i = 1:n_step
    alpha = (i - 1) / (n_step - 1); % interpolation fraction 0 → 1

    % Interpolate node positions
    p1 = (1 - alpha) * p1_orig + alpha * p1_def;
    p2 = (1 - alpha) * p2_orig + alpha * p2_def;
    p3 = (1 - alpha) * p3_orig + alpha * p3_def;

    % Compute physical strain
    edge1_i = norm(p1 - p2);
    edge2_i = norm(p2 - p3);
    edge3_i = norm(p1 - p3);
    strain(i) = (((edge1_i + edge2_i + edge3_i) / 3) - edgeLen)/edgeLen;

    % Compute energy for current configuration
    try
        [~, energy_b, energy_s] = deform_triangle(p1, p2, p3, edgeLen, l1, l4, beta, t, 0);
        U_vec(i) = energy_b + energy_s;  
        energy_b_vec(i) = energy_b;
        energy_s_vec(i) = energy_s;
    catch ME
        warning('Step %.3f failed: %s', alpha, ME.message);
        U_vec(i) = NaN;
    end
end
