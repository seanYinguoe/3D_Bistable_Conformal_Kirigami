function [opt_l2, opt_l3, max_bistability] = unit_design(desired_disp, edgeLen)
    % Parameters
    tol = 0.1;  % Tolerance for displacement match
    options = optimoptions('patternsearch', ...
        'Display', 'iter', ...
        'UseCompletePoll', true, ...
        'MaxIterations', 30);
    
    % Set bounds based on geometric constraint
    lb = [0, 0];
    ub = [0.85*edgeLen/3, 0.85*edgeLen];  % Max l2 = 0.85*edgeLen/3, max l3 = 0.85*edgeLen
    
    % Define equality constraint: l3 + 3*l2 = 0.85*edgeLen
    Aeq = [3, 1];          % 3*l2 + 1*l3
    beq = 0.85 * edgeLen;  % Target value
    
    % Define objective function
    objective = @(x) inverse_objective(x, desired_disp, tol, edgeLen);
    
    % Feasible initial guess
    l2_0 = 0.1 * edgeLen;
    l3_0 = 0.85 * edgeLen - 3*l2_0;
    x0 = [l2_0, l3_0];
    
    % Run optimization with equality constraint
    [x_opt, ~] = patternsearch(objective, x0, [], [], Aeq, beq, lb, ub, [], options);
    
    % Extract results
    opt_l2 = x_opt(1);
    opt_l3 = x_opt(2);
    [~, max_bistability] = inverse_objective(x_opt, desired_disp, tol, edgeLen);
end

function [obj_value, disp_bist] = inverse_objective(x, desired_disp, tol, edgeLen)
    l2 = x(1);
    l3 = x(2);
    
    % Set delta range proportional to edge length
    delta_range = linspace(0, edgeLen, 1000);
    
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
