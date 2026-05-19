% 1. Define your angles and strain
i = 30;
a1 = data.a1(i); 
a2 = data.a2(i); 
a3 = data.a3(i);   % radians
beta = data.beta(i);
eps = data.eps_bist(i);   % strain (eps_bist from ternary study)

% 2. Compute lambdas
lam3 = 1 + eps;
lam1 = lam3 * sin(a1)/sin(a3);
lam2 = lam3 * sin(a2)/sin(a3);
edgeLen = 15;
l1 = 0.80*edgeLen;
l4 = 0.05*edgeLen;
t = 0.015*edgeLen;

% 3. Convert to geometry
[q1, q2, q3] = scale_facs_to_q([lam1, lam2, lam3], edgeLen);

% 4. Plot energy profile
bistability_analysis(l1, l4, beta, t, edgeLen, q1, q2, q3, 1);
%  OR
%deform_triangle_anisotropic(q1, q2, q3, edgeLen, l1, l4, beta, t, nD, N, true);
