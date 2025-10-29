function [E_one, pack] = energy_lig(x0, N, p)
%% ENERGY: Single ligament energy minimization (standalone)
% Inputs:
%   th : initial theta guess
%   N  : number of segments
%   B0 : Left node
%   p  : struct with fields {t, E, b, beta, G, r_vertex, l2, l3, alphaL}

wrap = @(a) atan2(sin(a), cos(a));

%% Material and stiffness
t_eff = p.t * sqrt(3)/2;
EI = p.E * (p.b * t_eff^3) / 12;
EA = p.E * (p.b * t_eff);

%% Clustering
a0vec = href_end_cluster(p.l2, N, 0.9, 2.0);
%a0vec = p.l2/N * ones(1,N);
K = N - 1;
hhinge = 0.5*(a0vec(1:end-1) + a0vec(2:end));
kb_vec = EI ./ hhinge;
ks_vec = EA ./ a0vec;

%% Bounds
lb = [ -pi*ones(K,1);   -0.6*a0vec;   -2*pi ];
ub = [  pi*ones(K,1);    0.6*a0vec;    2*pi ];

%% Objective function
    function [f, g] = obj_fun(x)
        phi = x(1:K);
        e   = x(K+1:K+N);
        f = 0.5*( phi.'*(kb_vec.*phi) + e.'*(ks_vec.*e) );
        if nargout > 1
            g = [kb_vec.*phi; ks_vec.*e; 0];
        end
    end

%% Constraint function
    function [c, ceq, gc, gceq] = cons_fun(x)
        phi = x(1:K);
        e   = x(K+1:K+N);
        theta = x(end);

        % Geometry
        psi = zeros(N,1);
        psi(1) = p.alphaL;
        psi(2:end) = p.alphaL + cumsum(phi);

        u  = [cos(psi), sin(psi)];
        up = [-sin(psi), cos(psi)];

        % Position closure
        res_pos = sum(((a0vec + e).*u), 1).' - ( p.G(:) + p.r_vertex*[cos(theta); sin(theta)] - p.B_left(:));
        res_ang = sum(phi) - theta - pi + p.beta;

        ceq = [res_pos; res_ang];
        c = []; gc = [];

        % Jacobian (transpose for fmincon)
        Jpos_phi = zeros(2, K);
        for i = 1:K
            idx = (i+1):N;
            Jpos_phi(:, i) = sum( (a0vec(idx)+e(idx)) .* up(idx,:), 1 ).';
        end

        uthp = [-sin(theta); cos(theta)];

        Jpos_e  = (u.');
        Jpos_th = (-p.r_vertex * uthp);

        Jang_phi = ones(1, K);
        Jang_e   = zeros(1, N);
        Jang_th  = -1;

        J = [Jpos_phi, Jpos_e, Jpos_th;
            Jang_phi, Jang_e, Jang_th];

        gceq = J.';  % fmincon wants nvars×neq
    end

%% Solve optimization
opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true, ...
    'Display','off', ...
    'MaxIterations',300, ...
    'OptimalityTolerance',1e-12, ...
    'ConstraintTolerance',1e-12, ...
    'StepTolerance',1e-12);

[x_opt, fval] = fmincon(@obj_fun, x0, [], [], [], [], lb, ub, @cons_fun, opts);

%% Output
phi   = x_opt(1:K);
e     = x_opt(K+1:K+N);
theta = x_opt(end);
r = p.r_vertex;
B = p.G + r * [cos(theta), sin(theta)];
C = p.G + r * [cos(theta + 2*pi/3), sin(theta + 2*pi/3)];

E_one = fval;
pack = struct('phi',phi,'e',e,'theta',theta,...
    'B',B,'C',C,'G',p.G);

    function href = href_end_cluster(L0, N, r, p)
        % End-clustered segment lengths that sum to L0
        if nargin < 4, p = 2.0; end
        r = max(1e-6, min(0.9999, r));
        s = ((1:N)' - 0.5)/N;
        d = min(s, 1 - s);
        w = ((1 - r) + d).^p;
        href = (w / sum(w)) * L0;
    end

%% Plot the deformed ligament
ell = a0vec + e;                         % current lengths 
psi = zeros(N,1);                        % absolute segment angles
psi(1) = p.alphaL;
psi(2:end) = p.alphaL + cumsum(phi);

% % Deformed chain
% XY = zeros(N+1,2);
% XY(1,:) = p.B_left(:).';
% for i = 1:N
%     XY(i+1,:) = XY(i,:) + ell(i)*[cos(psi(i)), sin(psi(i))];
% end
% % undefomred chain
% XY0 = zeros(N+1,2);
% XY0(1,:) = p.B_left(:).';
% for i = 1:N
%     XY0(i+1,:) = XY0(i,:) + a0vec(i)*[cos(p.alphaL), sin(p.alphaL)];
% end
% figure; hold on; box on; axis equal;
% plot(XY(:,1),  XY(:,2),  '-o', 'LineWidth', 1.6, 'MarkerSize', 4, ...
%     'DisplayName', 'deformed');
% plot(XY0(:,1),  XY0(:,2),  '--', 'LineWidth', 1.0,'DisplayName', 'undeformed');
% hold off;
end
