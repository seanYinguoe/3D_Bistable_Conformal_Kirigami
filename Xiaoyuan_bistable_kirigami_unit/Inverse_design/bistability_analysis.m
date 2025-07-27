function [disp_bist, bistability] = bistability_analysis(l1, l2, l3, t, edgeLen, q1, q2, q3)
% Compute bistable displacement and energy barrier of single unit using deform_triangle energy model
%
% Inputs:
%   l1, l2, l3 : Geometric parameters (flank, filament, inner triangle lengths)
%   q1, q2, q3 : Target deformed boundary vertices
%   t          : Filament width
%   edgeLen    : Original edge length of outer triangle
%
% Outputs:
%   disp_bist  : Physical displacement (mean edge change) at bistable state
%   bistability: Energy barrier between local max and local min

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
        [~, energy_b, energy_s] = deform_triangle(p1, p2, p3, edgeLen, l1, l2, l3, t, 0);
        U_vec(i) = energy_b + energy_s;
        energy_b_vec(i) = energy_b;
        energy_s_vec(i) = energy_s;
    catch ME
        warning('Step %.3f failed: %s', alpha, ME.message);
        U_vec(i) = NaN;
    end
end

%% Compute derivative to find local maxima and minima
% dU = gradient(U_vec, displacement);
% critical_points = struct('delta', {}, 'U', {}, 'type', {});
% 
% for i = 2:length(dU)-1
%     if dU(i-1) > 0 && dU(i) < 0
%         critical_points(end+1) = struct('delta', displacement(i), 'U', U_vec(i), 'type', 'Local Maximum');
%     elseif dU(i-1) < 0 && dU(i) > 0
%         critical_points(end+1) = struct('delta', displacement(i), 'U', U_vec(i), 'type', 'Local Minimum');
%     end
% end

%% Remove duplicates
% [~, unique_idx] = unique([critical_points.delta]);
% critical_points = critical_points(unique_idx);

% %% Classify bistability
% max_points = critical_points(strcmp({critical_points.type}, 'Local Maximum'));
% min_points = critical_points(strcmp({critical_points.type}, 'Local Minimum'));
% min_points = min_points([min_points.delta] > 0.0001); % Ignore near-zero min
% 
% if length(max_points) == 1 && ~isempty(min_points)
%     disp_bist = min_points(1).delta;
%     bistability = max_points(1).U - min_points(1).U;
%     fprintf('Bistable system detected:\n');
%     fprintf('  Local Maximum: δ = %.4f, U = %.4f\n', max_points(1).delta, max_points(1).U);
%     fprintf('  Local Minimum: δ = %.4f, U = %.4f\n', disp_bist, min_points(1).U);
%     fprintf('  Energy Barrier: %.4f\n', bistability);
% else
%     disp('This is not a bistable unit.');
%     disp_bist = NaN;
%     bistability = 0;
% end
% 
%% Plot energy vs physical displacement
figure;
plot(strain, U_vec, 'b-', 'LineWidth', 2); hold on;
plot(strain, energy_b_vec, '--r', 'LineWidth', 1.5);
plot(strain, energy_s_vec, '--g', 'LineWidth', 1.5);
legend('Total Energy', 'Bending Energy', 'Stretch Energy', 'Location', 'best');
xlabel('\delta (Mean Edge Length Change)');
ylabel('Energy');
title('Energy-Displacement Curve');
grid on;
% 
% % Mark critical points on plot
% for i = 1:length(critical_points)
%     plot(critical_points(i).delta, critical_points(i).U, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
%     text(critical_points(i).delta, critical_points(i).U, sprintf(' %s', critical_points(i).type), ...
%         'VerticalAlignment', 'bottom', 'FontSize', 8);
% end

end
