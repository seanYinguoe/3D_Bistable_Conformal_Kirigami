function bistability_analysis(l2, l3)
syms x1 x2 delta real
k = 1; % stiffness constant

% Potential energy as function of x1 and x2
U = 0.5 * k * (x1 - pi/3)^2 + 0.5 * k * (x2 - 2*pi/3)^2;

% Constraint equations
f1 = l2*(-sqrt(3) + sin(2*pi/3 - x1) + sin(x1)) + ...
    l3*(-sqrt(3)/2 + sin(x1 + x2 - 2*pi/3));

f2 = l2*(-1 - cos(2*pi/3 - x1) + cos(x1)) + ...
    l3*(-1/2 + cos(x1 + x2 - 2*pi/3)) + (delta + l2);

% Implicitly solve constraints: x1 = x1(delta), x2 = x2(delta)
[solx1, solx2] = solve([f1 == 0, f2 == 0], [x1, x2], 'Real', true, 'IgnoreAnalyticConstraints', true);

% Loop through each solution
for i = 1:length(solx1)
    x1_sol = simplify(solx1(i));
    x2_sol = simplify(solx2(i));

    % Substitute into energy function
    U_delta = simplify(subs(U, [x1, x2], [x1_sol, x2_sol]));
    dU_ddelta = simplify(diff(U_delta, delta));

    fprintf("Solution %d:\n", i);
    fprintf("U(delta) = "); disp(U_delta);
    fprintf("dU/d(delta) = "); disp(dU_ddelta);
    fprintf("\n");
end
end
