function [opt_l2, opt_l3, max_bistability] = unit_design(desired_disp, delta_range, lb, ub, edgeLen)
% Initial geometry guess



% Parameters
tol = 0.1;  % Tolerance for displacement match
options = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'UseCompletePoll', true, ...
    'MaxIterations', 30);

% Define objective function
objective = @(x) inverse_objective(x, desired_disp, tol, delta_range);

% Initial guess (centered within bounds)
x0 = [mean(lb(1), ub(1)), mean(lb(2), ub(2))];

% Run optimization
[x_opt, ~] = patternsearch(objective, x0, [], [], [], [], lb, ub, [], options);

% Extract results
opt_l2 = x_opt(1);
opt_l3 = x_opt(2);
[~, max_bistability] = inverse_objective(x_opt, desired_disp, tol, delta_range);
end

function [obj_value, disp_bist] = inverse_objective(x, desired_disp, tol, delta_range)
l2 = x(1);
l3 = x(2);

% Run forward analysis
[disp_bist, bistability] = bistability_analysis(l2, l3, delta_range);

% Penalty for non-bistable or displacement mismatch
if isempty(disp_bist)
    obj_value = -1e10;  % Large penalty for non-bistable
else
    displacement_error = abs(disp_bist - desired_disp);
    if displacement_error <= tol
        obj_value = bistability;  % Maximize bistability
    else
        obj_value = -1e10 * displacement_error;  % Penalize displacement mismatch
    end
end
end
