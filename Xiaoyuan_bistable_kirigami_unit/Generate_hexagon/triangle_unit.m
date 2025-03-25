% Create single bistable unit:
% Input parameters:
% l1: Length of flanks
% l2: thickness of flanks
% l3: length of filaments
% l4: length of inner triangle
% t: thickness of filaments
% theta: deploying angle
% Output:
% Deployed unit

% Using fmincon to solve the constrained optimisation problem
% There are some implicit relationship among alpha_1, alpha_2 and delta
% Define the desired delta

function [triangle,alpha_1_optimal,alpha_2_optimal] = triangle_unit(prev_alpha_1, prev_alpha_2, delta,l1,l2,l3,t)
% Adjust this value to control the motion range

% Define the constraint function
constraint_func = @(params) constraint_function(params,l1,l2,l3,t);

% Define the objective function
objective = @(params) (params(3) - delta)^2;


% Use previous values if available, otherwise use default initial guess
if nargin < 2
    x0 = [pi/3, 2*pi/3, delta];
else
    x0 = [prev_alpha_1, prev_alpha_2, delta];
end

% Lower and upper bounds for parameters
lb = [0,0,0];
ub = [pi,pi,1];

% Options for fmincon
options = optimoptions('fmincon', 'Display', 'final-detailed', 'Algorithm', 'sqp', ...
    'ConstraintTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12, ...
    'StepTolerance', 1e-10,...
    'MaxFunctionEvaluations', 10000);
  

% Solve the constrained optimization problem
[optimal_params, fval, ~, ~] = fmincon(@(params) objective(params), x0, [], [], [], [], lb, ub, constraint_func, options);
display(fval);
[~, ceq] = constraint_function(optimal_params,l1,l2,l3,t);
display(ceq)

% Extract the optimal parameters
alpha_1_optimal = optimal_params(1);
alpha_2_optimal = optimal_params(2);

% Generate the final triangle
triangle = triangle_expression(alpha_1_optimal, alpha_2_optimal,l1,l2,l3,t);

    function [c, ceq] = constraint_function(params,l1,l2,l3,t)
        alpha_1 = params(1);
        alpha_2 = params(2);
        delta = params(3);
        triangle = triangle_expression(alpha_1,alpha_2,l1,l2,l3,t);
        flank_3 = triangle(27:30,:);


        % Define inequality constraints
        %     ceq = [
        %         (1/2)*(-sqrt(3)*(2*l2 + l3) + sqrt(3)*l2*cos(alpha_1) + 3*l2*sin(alpha_1) - l3*(sqrt(3)*cos(alpha_1 + alpha_2) + sin(alpha_1 + alpha_2)));
        %         (1/2)*(-l3 + 3*l2*cos(alpha_1) - l3*cos(alpha_1 + alpha_2) - sqrt(3)*l2*sin(alpha_1) + sqrt(3)*l3*sin(alpha_1 + alpha_2)) - (-delta)
        %     ]; some problem with this equations from the paper should be fixed
        l4 = l1 - 2*l2 - l3;
        ceq = [flank_3(3,1);
            flank_3(3,2) + (l2+l4+delta)
            ];
        % Define inequality constraints
        c = [];
    end

end