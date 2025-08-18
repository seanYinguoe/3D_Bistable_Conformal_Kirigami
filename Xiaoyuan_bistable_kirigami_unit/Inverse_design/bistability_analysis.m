function [strain_bist, bistability] = bistability_analysis(l1, l4, beta, t, edgeLen, q1, q2, q3)
% Compute bistable displacement and energy barrier of single unit using deform_triangle energy model
%
% Inputs:
%   l1, l4, beta : Geometric parameters (check reference graph)
%   q1, q2, q3 : Target deformed boundary vertices
%   t          : Filament width
%   edgeLen    : Original edge length of outer triangle
%
% Outputs:
%   strain_bist  : Physical displacement (mean edge change) at bistable state
%   bistability  : Energy barrier between local max and local min

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
        [~, energy_b, energy_s] = deform_triangle(p1, p2, p3, edgeLen, l1, l4, beta,t,0);
        U_vec(i) = energy_b + energy_s;
        energy_b_vec(i) = energy_b;
        energy_s_vec(i) = energy_s;
    catch ME
        warning('Step %.3f failed: %s', alpha, ME.message);
        U_vec(i) = NaN;
    end
end

%% Compute derivative to find local maxima and minima
dU = gradient(U_vec, strain);
critical_points = struct('strain', {}, 'U', {}, 'type', {}, 'index', {});

for i = 2:length(dU)-1
    if dU(i-1) > 0 && dU(i) < 0
        critical_points(end+1) = struct('strain', strain(i), 'U', U_vec(i), 'type', 'Local Maximum', 'index', i);
    elseif dU(i-1) < 0 && dU(i) > 0
        critical_points(end+1) = struct('strain', strain(i), 'U', U_vec(i), 'type', 'Local Minimum', 'index', i);
    end
end

%% Keep only first local max and first local min after it
first_max_idx = find(strcmp({critical_points.type}, 'Local Maximum'), 1);
first_min_idx = find(strcmp({critical_points.type}, 'Local Minimum') & [critical_points.index] > critical_points(first_max_idx).index, 1);

strain_bist = NaN;
bistability = 0;

if ~isempty(first_max_idx) && ~isempty(first_min_idx)
    max1 = critical_points(first_max_idx);
    min1 = critical_points(first_min_idx);
    strain_bist = min1.strain;
    bistability = (max1.U - min1.U)/max1.U;
end

%% Plot energy vs physical displacement
figure;
set(gca, 'FontSize', 18);  % applies to both x and y tick labels
plot(strain, U_vec, 'b-', 'LineWidth', 2); hold on;
plot(strain, energy_b_vec, '--r', 'LineWidth', 1.5);
plot(strain, energy_s_vec, '--g', 'LineWidth', 1.5);
xlabel('Mean Strain','FontSize', 20);
ylabel('Energy','FontSize', 20);
title('Energy-Strain Curve','FontSize', 20);
grid on;

% Mark only max1 and min1
if exist('max1', 'var') && exist('min1', 'var')
    plot(max1.strain, max1.U, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
    text(max1.strain, max1.U, ' localmax', 'VerticalAlignment', 'bottom', 'FontSize', 12);

    plot(min1.strain, min1.U, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
    text(min1.strain, min1.U, ' localmin', 'VerticalAlignment', 'bottom', 'FontSize', 12);
end

legend({'Total Energy','Bending Energy','Stretch Energy'}, 'Location', 'best');
end
