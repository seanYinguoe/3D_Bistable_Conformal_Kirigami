function [triangle,alpha_1_optimal,alpha_2_optimal] = triangle_unit(prev_alpha_1, prev_alpha_2, delta,l1,l2,l3,t)
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
[alpha, fval] = fsolve(@(x) equations(x, l2, l3, t, delta), x0, options);

display(fval);

% Extract the optimal parameters
alpha_1_optimal = alpha(1);
alpha_2_optimal = alpha(2);

% Generate the final triangle
triangle = triangle_expression(alpha_1_optimal, alpha_2_optimal,l1,l2,l3,t);

function F = equations(x, l2, l3, t, delta)
    
    alpha_1 = x(1);
    alpha_2 = x(2);
%     F = [
%         (1/2) * (-sqrt(3)*(2*l2 + l3) + sqrt(3)*l2*cos(alpha_1) + ...
%         3*l2*sin(alpha_1) - l3*(sqrt(3)*cos(alpha_1 + alpha_2) + sin(alpha_1 + alpha_2)));
%         (1/2) * (-l3 + 3*l2*cos(alpha_1) - l3*cos(alpha_1 + alpha_2) - ...
%         sqrt(3)*l2*sin(alpha_1) + sqrt(3)*l3*sin(alpha_1 + alpha_2)) + delta
%         ];
    l2_ = l2 + t*sin(alpha_2 + alpha_1)/sin(alpha_1);
    %l2_ = l2;
    F = [
        l2 * (-2*sin(pi/3) + sin(2*pi/3-alpha_1) + l2_/l2*sin(pi/3))+...
        l3 * (-sin(pi/3) + sin(alpha_1 + alpha_2 - 2*pi/3));
        l2 * (-2*cos(pi/3) - cos(2*pi/3-alpha_1) + l2_/l2*cos(pi/3))+...
        l3 * (-cos(pi/3) + cos(alpha_1 + alpha_2 - 2*pi/3)) + l2 + delta;
        ];
end


end