function [disp_bist, bistability] = bistability_analysis(l2, l3, delta_range)
%% Parameters
k = 1;
rest_alpha1 = pi/3;
rest_alpha2 = 2*pi/3;

%% Define system functions
f1 = @(a1, a2) l2*(-sqrt(3) + sin(2*pi/3 - a1) + sin(a1)) + ...
    l3*(-sqrt(3)/2 + sin(a1 + a2 - 2*pi/3));
f2 = @(a1, a2) l2*(-1 - cos(2*pi/3 - a1) + cos(a1)) + ...
    l3*(-0.5 + cos(a1 + a2 - 2*pi/3)) + l2;
U = @(a1, a2) 0.5*k*(a1 - rest_alpha1).^2 + 0.5*k*(a2 - rest_alpha2).^2;

%% Solve system for U(delta)
U_vec = nan(size(delta_range));
a1_guess = rest_alpha1;
a2_guess = rest_alpha2;
options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-10);

for i = 1:length(delta_range)
    fun = @(x) [f1(x(1), x(2));
        f2(x(1), x(2)) + delta_range(i)];
    [sol, ~, exitflag] = fsolve(fun, [a1_guess; a2_guess], options);
    if exitflag > 0
        U_vec(i) = U(sol(1), sol(2));
        a1_guess = sol(1);
        a2_guess = sol(2);
    end
end

%% Identify critical points using first derivative
dU_ddelta = gradient(U_vec, delta_range);
critical_points = struct('delta', {}, 'U', {}, 'type', {});

% Find sign changes indicating critical points
for i = 2:length(dU_ddelta)-1
    % Local maximum detection (positive to negative)
    if dU_ddelta(i-1) > 0 && dU_ddelta(i) <= 0 && dU_ddelta(i+1) < 0
        [~, idx] = min(abs(dU_ddelta(max(1,i-10):min(i+10,length(delta_range)))));
        idx = idx + max(1,i-10) - 1;
        critical_points(end+1) = struct('delta', delta_range(idx), ...
            'U', U_vec(idx), ...
            'type', 'Local Maximum');

        % Local minimum detection (negative to positive)
    elseif dU_ddelta(i-1) < 0 && dU_ddelta(i) >= 0 && dU_ddelta(i+1) > 0
        [~, idx] = min(abs(dU_ddelta(max(1,i-10):min(i+10,length(delta_range)))));
        idx = idx + max(1,i-10) - 1;
        critical_points(end+1) = struct('delta', delta_range(idx), ...
            'U', U_vec(idx), ...
            'type', 'Local Minimum');
    end
end

%% Remove duplicates
[~, unique_idx] = unique([critical_points.delta]);
critical_points = critical_points(unique_idx);

%% Classify bistability
max_points = critical_points(strcmp('Local Maximum', {critical_points.type}));
min_points = critical_points(strcmp('Local Minimum', {critical_points.type}));

% Remove initial state (delta=0) from min_points
min_points = min_points([min_points.delta] > 0.1);

% Check bistability conditions
if length(max_points) == 1 && ~isempty(min_points)
    % Valid bistable system
    disp_bist = min_points(1).delta;  % Use first post-zero minimum
    bistability = max_points(1).U - min_points(1).U;
    fprintf('Bistable system detected:\n');
    fprintf('  Local Maximum: δ = %.4f, U = %.4f\n', max_points(1).delta, max_points(1).U);
    fprintf('  Bistable State: δ = %.4f, U = %.4f\n', disp_bist, min_points(1).U);
    fprintf('  Energy Barrier: %.4f\n', bistability);
else
    % Not bistable
    disp_bist = [];
    bistability = 0;
    disp('This is not a bistable unit.');
end
end
