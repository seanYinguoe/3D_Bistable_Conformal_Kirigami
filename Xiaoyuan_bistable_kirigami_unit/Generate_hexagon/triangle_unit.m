function [triangle,alpha_1_optimal,alpha_2_optimal] = triangle_unit(prev_alpha_1, prev_alpha_2, beta, delta,l1,l2,l3,t)
% Create single bistable unit:
% Input parameters:
% l1: Length of flanks
% l2: length of filaments
% l3: length of inner triangle
% t: thickness of filaments
% theta: deploying angle
% beta: tilting angle
% Output:
% Deployed unit
% There are some implicit relationship among alpha_1, alpha_2 and delta
% Define the desired delta

% Use previous values if available, otherwise use default initial guess
if nargin < 2
    x0 = [pi/3, 2*pi/3];
else
    x0 = [prev_alpha_1, prev_alpha_2];
end
  

% Solve the system
options = optimoptions('fsolve', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt');
[alpha, fval] = fsolve(@(x) equations(x, l1, l2, l3, beta, delta), x0, options);

display(fval);

% Extract the optimal parameters
alpha_1_optimal = alpha(1);
alpha_2_optimal = alpha(2);

% Generate the final triangle
triangle = triangle_expression(alpha_1_optimal, alpha_2_optimal,l1,l2,l3,t,beta);

function F = equations(x, l1, l2, l3, beta, delta)
    
    alpha_1 = x(1);
    alpha_2 = x(2);
    l4 = l1 - l3 - (1+cos(pi/3)+sin(pi/3)*cos(pi/3-beta)/sin(pi/3-beta))* l2;
    % l2_ = l2;
    % F = [
    %     l2 * (-2*sin(pi/3) + sin(2*pi/3-alpha_1) + l2_/l2*sin(alpha_1))+...
    %     l3 * (-sin(pi/3) + sin(alpha_1 + alpha_2 - 2*pi/3));
    %     l2 * (-2*cos(pi/3) - cos(2*pi/3-alpha_1) + l2_/l2*cos(alpha_1))+...
    %     l3 * (-cos(pi/3) + cos(alpha_1 + alpha_2 - 2*pi/3)) + l2 + delta;
    %     ];
    F = [
        (l4 - l1)*cos(pi/6 + beta) ...
        + l2*cos(beta + pi/6 - alpha_1) ...
        - l3*cos(pi/6 + beta - alpha_1 - alpha_2) ...
        + l2*cos(pi/6 + beta - alpha_1 + pi/3);                          % = 0

        (l4 - l1)*sin(pi/6 + beta) ...
        + l2*sin(beta +pi/6 - alpha_1) ...
        - l3*sin(pi/6 + beta - alpha_1 - alpha_2) ...
        + l2*sin(pi/6 + beta - alpha_1 + pi/3) ...
        - ( -delta - (l2/sin(pi/3 - beta))*sin(pi/3) );                   % = 0
        ];
end


end