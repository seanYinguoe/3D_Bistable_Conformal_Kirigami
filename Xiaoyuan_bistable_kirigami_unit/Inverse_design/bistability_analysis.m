syms x1 x2 d lambda1 lambda2
k = 1; l2 = 1.5; l3 = 8.25;

f = 0.5 * k * (x1 - pi/3)^2 + 0.5 * k * (x2 - 2*pi/3)^2;
g1 = l2*(-sqrt(3) + sin(2*pi/3 - x1) + sin(x1)) + ...
     l3*(-sqrt(3)/2 + sin(x1 + x2 - 2*pi/3));
g2 = l2*(-1 - cos(2*pi/3 - x1) + cos(x1)) + ...
     l3*(-1/2 + cos(x1 + x2 - 2*pi/3)) + (d + l2);
L = f - lambda1 * g1 - lambda2 * g2;

eqns = [
    diff(L, x1) == 0;
    diff(L, x2) == 0;
    diff(L, d) == 0;
    g1 == 0;
    g2 == 0
];

S = solve(eqns, [x1, x2, d, lambda1, lambda2], 'Real', true, 'IgnoreAnalyticConstraints', true);
disp(S)
