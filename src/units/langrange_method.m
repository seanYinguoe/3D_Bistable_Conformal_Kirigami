
%% Parameters
k = 1;
l2 = 1.5;
l3 = 8.25;
rest_alpha1 = pi/3;
rest_alpha2 = 2*pi/3;
delta_min = 0;
delta_max = 9;
num_points = 1000; % High resolution for accurate sign detection

%% Define functions
f1 = @(a1, a2) l2*(-sqrt(3) + sin(2*pi/3 - a1) + sin(a1)) + ...
                l3*(-sqrt(3)/2 + sin(a1 + a2 - 2*pi/3));
f2 = @(a1, a2) l2*(-1 - cos(2*pi/3 - a1) + cos(a1)) + ...
                l3*(-0.5 + cos(a1 + a2 - 2*pi/3)) + l2;
U = @(a1, a2) 0.5*k*(a1 - rest_alpha1).^2 + 0.5*k*(a2 - rest_alpha2).^2;

%% Solve system for U(delta)
delta_vec = linspace(delta_min, delta_max, num_points);
U_vec = nan(size(delta_vec));
a1_guess = rest_alpha1;
a2_guess = rest_alpha2;

options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-10);

for i = 1:length(delta_vec)
    fun = @(x) [f1(x(1), x(2)); 
                 f2(x(1), x(2)) + delta_vec(i)];
    [sol, ~, exitflag] = fsolve(fun, [a1_guess; a2_guess], options);
    if exitflag > 0
        U_vec(i) = U(sol(1), sol(2));
        a1_guess = sol(1);
        a2_guess = sol(2);
    end
end

%% Compute dU/ddelta and classify critical points
dU_ddelta = gradient(U_vec, delta_vec);
critical_points = [];

% Detect sign changes in dU/ddelta
for i = 2:length(dU_ddelta)-1
    % Sign change: positive to negative (local maximum)
    if dU_ddelta(i-1) > 0 && dU_ddelta(i) <= 0 && dU_ddelta(i+1) < 0
        [~, idx] = min(abs(dU_ddelta(max(1,i-10):min(i+10,num_points))));
        idx = idx + max(1,i-10) - 1;
        critical_points = [critical_points; 
                          struct('delta', delta_vec(idx), ...
                                 'U', U_vec(idx), ...
                                 'type', 'Local Maximum')];
    
    % Sign change: negative to positive (local minimum)
    elseif dU_ddelta(i-1) < 0 && dU_ddelta(i) >= 0 && dU_ddelta(i+1) > 0
        [~, idx] = min(abs(dU_ddelta(max(1,i-10):min(i+10,num_points))));
        idx = idx + max(1,i-10) - 1;
        critical_points = [critical_points; 
                          struct('delta', delta_vec(idx), ...
                                 'U', U_vec(idx), ...
                                 'type', 'Local Minimum')];
    end
end

%% Remove duplicates and sort
[~, unique_idx] = unique([critical_points.delta], 'stable');
critical_points = critical_points(unique_idx);

%% Display results
disp('Critical Points (First Derivative Test):');
for i = 1:length(critical_points)
    fprintf('  %s at δ = %.4f, U = %.4f\n', ...
            critical_points(i).type, ...
            critical_points(i).delta, ...
            critical_points(i).U);
end

%% Plot results
figure;
plot(delta_vec, U_vec, 'b-', 'LineWidth', 1.5);
hold on;
for i = 1:length(critical_points)
    plot(critical_points(i).delta, critical_points(i).U, 'ro', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'r');
    text(critical_points(i).delta, critical_points(i).U, ...
         sprintf('  %s\n  δ=%.3f, U=%.3f', critical_points(i).type, ...
                 critical_points(i).delta, critical_points(i).U), ...
         'VerticalAlignment', 'bottom', 'FontSize', 9);
end
xlabel('\delta');
ylabel('U(\delta)');
title('Energy-Displacement Curve with Critical Points');
grid on;
legend('U(\delta)', 'Critical Points', 'Location', 'best');
