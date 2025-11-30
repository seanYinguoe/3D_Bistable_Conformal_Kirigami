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


%% Get the E_total energy and alpha(alpha is the interpolation)
nD = 200;
N = 8;
[E_total,alpha] = deform_triangle_anisotropic(q1,q2,q3,edgeLen,l1,l4,beta,t,nD,N);

%% Compute derivative to find local maxima and minima
dU = gradient(E_total, alpha);
critical_points = struct('strain', {}, 'U', {}, 'type', {}, 'index', {});
strain = alpha*(norm(q1-q2)/edgeLen-1);

for i = 2:length(dU)-1
    if dU(i-1) > 0 && dU(i) < 0
        critical_points(end+1) = struct('strain', strain(i), 'U', E_total(i), 'type', 'Local Maximum', 'index', i);
    elseif dU(i-1) < 0 && dU(i) > 0
        critical_points(end+1) = struct('strain', strain(i), 'U', E_total(i), 'type', 'Local Minimum', 'index', i);
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
plot(strain, E_total, 'b-', 'LineWidth', 2); hold on;
xlabel('Strain','FontSize', 20);
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

legend({'Total Energy'}, 'Location', 'best');
end
