function [E_one, pack] = energy_lig_anisotropic(x0, N, p)
%% ENERGY: Single ligament energy minimization (standalone)
% Inputs:
%   th : initial theta guess
%   N  : number of segments
%   B0 : Left node
%   p = struct('t',t,'E',Emod,'b',b,'beta',beta,'B_flank',B_flank,...
%        'A_flank',A_flank,'C_flank',C_flank,'r_vertex',r_vertex,'l2',l2,...
%        'l3',l3,'alphaL_B',alphaL_B,'alphaL_A',alphaL_A,'alphaL_C',alphaL_C,...
%        'xG0',xG0,'yG0',yG0);

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
xG0 = p.xG0;
yG0 = p.yG0;
G0 = [xG0, yG0];

% Cache constant index/layout blocks (no math change)
i_phiB = 1:K;
i_eB   = K + (1:N);
i_phiA = (K+N) + (1:K);
i_eA   = (K+N+K) + (1:N);
i_phiC = (K+N+K+N) + (1:K);
i_eC   = (K+N+K+N+K) + (1:N);
i_th   = 3*K + 3*N + 1;
i_xG   = 3*K + 3*N + 2;
i_yG   = 3*K + 3*N + 3;
r_posB = 1:2;   r_angB = 3;
r_posA = 4:5;   r_angA = 6;
r_posC = 7:8;   r_angC = 9;
nvars  = 3*K + 3*N + 3;
negI2  = -eye(2);
onesK  = ones(1,K);
zerosN = zeros(1,N);

%% Bounds
lb = [ -pi*ones(K,1);   -0.6*a0vec;...
    -pi*ones(K,1);   -0.6*a0vec;...
    -pi*ones(K,1);   -0.6*a0vec;...
    -2*pi; xG0 - 0.2*p.l3; yG0 - 0.2*p.l3];
ub = [  pi*ones(K,1);    0.6*a0vec;...
    pi*ones(K,1);    0.6*a0vec;...
    pi*ones(K,1);    0.6*a0vec;...
    2*pi; xG0 + 0.2*p.l3; yG0 + 0.2*p.l3];

%% Objective function
    function [f, g] = obj_fun(x)
        phiB = x(1:K);
        eB   = x(K+1:K+N);
        phiA = x(K+N+1:2*K+N);
        eA   = x(2*K+N+1:2*K+2*N);
        phiC = x(2*K+2*N+1:3*K+2*N);
        eC   = x(3*K+2*N+1:3*K+3*N);
        f = 0.5*( phiB.'*(kb_vec.*phiB) + eB.'*(ks_vec.*eB) ) + ...
            0.5*( phiA.'*(kb_vec.*phiA) + eA.'*(ks_vec.*eA) ) + ...
            0.5*( phiC.'*(kb_vec.*phiC) + eC.'*(ks_vec.*eC) );
        if nargout > 1
            g = [kb_vec.*phiB; ks_vec.*eB;...
                kb_vec.*phiA; ks_vec.*eA;...
                kb_vec.*phiC; ks_vec.*eC;...
                0; 0; 0];
        end
    end

%% Constraint function
    function [c, ceq, gc, gceq] = cons_fun(x)
        phiB = x(1:K);
        eB   = x(K+1:K+N);
        phiA = x(K+N+1:2*K+N);
        eA   = x(2*K+N+1:2*K+2*N);
        phiC = x(2*K+2*N+1:3*K+2*N);
        eC   = x(3*K+2*N+1:3*K+3*N);
        theta = x(end-2);
        xG = x(end-1);
        yG = x(end);
        G = [xG; yG];

        % ---------- Geometry----------
        % B
        psiB = zeros(N,1); psiB(1) = p.alphaL_B; psiB(2:end) = p.alphaL_B + cumsum(phiB);
        uB  = [cos(psiB), sin(psiB)];
        upB = [-sin(psiB), cos(psiB)];
        % A
        psiA = zeros(N,1); psiA(1) = p.alphaL_A; psiA(2:end) = p.alphaL_A + cumsum(phiA);
        uA  = [cos(psiA), sin(psiA)];
        upA = [-sin(psiA), cos(psiA)];
        % C
        psiC = zeros(N,1); psiC(1) = p.alphaL_C; psiC(2:end) = p.alphaL_C + cumsum(phiC);
        uC  = [cos(psiC), sin(psiC)];
        upC = [-sin(psiC), cos(psiC)];

        % ---------- Constraints----------
        res_posB = sum(((a0vec + eB).*uB), 1).' - ( G(:) + p.r_vertex*[cos(theta);               sin(theta)]               - p.B_flank(:));
        res_angB = sum(phiB) + p.alphaL_B - (theta + 5*pi/6);

        res_posA = sum(((a0vec + eA).*uA), 1).' - ( G(:) + p.r_vertex*[cos(theta - 2*pi/3);      sin(theta - 2*pi/3)]      - p.A_flank(:));
        res_angA = sum(phiA) + p.alphaL_A - (theta + pi/6);

        res_posC = sum(((a0vec + eC).*uC), 1).' - ( G(:) + p.r_vertex*[cos(theta + 2*pi/3);      sin(theta + 2*pi/3)]      - p.C_flank(:));
        res_angC = sum(phiC) + p.alphaL_C - (theta + 3*pi/2);

        ceq = [res_posB; res_angB; ...
            res_posA; res_angA; ...
            res_posC; res_angC];
        c = []; gc = [];

        % d(pos)/d(e) = u^T
        JposB_e = uB.';   % 2xN
        JposA_e = uA.';   % 2xN
        JposC_e = uC.';   % 2xN

        % d(pos)/d(phi): sum_{j>i} (a0_j + e_j) * up_j
        JposB_phi = zeros(2,K);
        JposA_phi = zeros(2,K);
        JposC_phi = zeros(2,K);
        for i = 1:K
            idx = (i+1):N;
            JposB_phi(:,i) = sum( (a0vec(idx) + eB(idx)) .* upB(idx,:), 1 ).';
            JposA_phi(:,i) = sum( (a0vec(idx) + eA(idx)) .* upA(idx,:), 1 ).';
            JposC_phi(:,i) = sum( (a0vec(idx) + eC(idx)) .* upC(idx,:), 1 ).';
        end

        % d(pos)/d(theta): -r * R'(theta±shift) * [1;0]
        uthpB = [-sin(theta);             cos(theta)];
        uthpA = [-sin(theta - 2*pi/3);    cos(theta - 2*pi/3)];
        uthpC = [-sin(theta + 2*pi/3);    cos(theta + 2*pi/3)];
        JposB_th = -p.r_vertex * uthpB;   % 2x1
        JposA_th = -p.r_vertex * uthpA;   % 2x1
        JposC_th = -p.r_vertex * uthpC;   % 2x1

        % d(ang)/d(·)
        Jang_phiB = onesK;  Jang_eB = zerosN;  Jang_thB = -1;
        Jang_phiA = onesK;  Jang_eA = zerosN;  Jang_thA = -1;
        Jang_phiC = onesK;  Jang_eC = zerosN;  Jang_thC = -1;

        % Initialize J and place blocks
        J = zeros(9, nvars);

        % ---- B rows ----
        J(r_posB, i_phiB) = JposB_phi;
        J(r_posB, i_eB)   = JposB_e;
        J(r_posB, i_th)   = JposB_th;
        J(r_posB, [i_xG i_yG]) = negI2;                 % d res_posB / d(xG,yG) = -I2
        J(r_angB, i_phiB) = Jang_phiB;
        J(r_angB, i_eB)   = Jang_eB;
        J(r_angB, i_th)   = Jang_thB;

        % ---- A rows ----
        J(r_posA, i_phiA) = JposA_phi;
        J(r_posA, i_eA)   = JposA_e;
        J(r_posA, i_th)   = JposA_th;
        J(r_posA, [i_xG i_yG]) = negI2;                 % d res_posA / d(xG,yG) = -I2
        J(r_angA, i_phiA) = Jang_phiA;
        J(r_angA, i_eA)   = Jang_eA;
        J(r_angA, i_th)   = Jang_thA;

        % ---- C rows ----
        J(r_posC, i_phiC) = JposC_phi;
        J(r_posC, i_eC)   = JposC_e;
        J(r_posC, i_th)   = JposC_th;
        J(r_posC, [i_xG i_yG]) = negI2;                 % d res_posC / d(xG,yG) = -I2
        J(r_angC, i_phiC) = Jang_phiC;
        J(r_angC, i_eC)   = Jang_eC;
        J(r_angC, i_th)   = Jang_thC;

        % fmincon needs nvars x neq（转置）
        gceq = J.';
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
phiB  = x_opt(1:K);
eB    = x_opt(K+1:K+N);
phiA  = x_opt(K+N+1:2*K+N);
eA    = x_opt(2*K+N+1:2*K+2*N);
phiC  = x_opt(2*K+2*N+1:3*K+2*N);
eC    = x_opt(3*K+2*N+1:3*K+3*N);
theta = x_opt(end-2);
xG    = x_opt(end-1);
yG    = x_opt(end);
G     = [xG, yG];                 % optimized centroid (质心)

r = p.r_vertex;
B = G + r * [cos(theta),               sin(theta)];
C = G + r * [cos(theta + 2*pi/3),      sin(theta + 2*pi/3)];
A = G + r * [cos(theta - 2*pi/3),      sin(theta - 2*pi/3)];

% Deformed chain
XYB = zeros(N+1,2);
XYB(1,:) = p.B_flank(:).';
ellB = a0vec + eB;                         % current lengths 
psiB = zeros(N,1);                        % absolute segment angles
psiB(1) = p.alphaL_B;
psiB(2:end) = p.alphaL_B + cumsum(phiB);
for j = 1:N
    XYB(j+1,:) = XYB(j,:) + ellB(j)*[cos(psiB(j)), sin(psiB(j))];
end

XYA = zeros(N+1,2);
XYA(1,:) = p.A_flank(:).';
ellA = a0vec + eA;                         % current lengths 
psiA = zeros(N,1);                        % absolute segment angles
psiA(1) = p.alphaL_A;
psiA(2:end) = p.alphaL_A + cumsum(phiA);
for j = 1:N
    XYA(j+1,:) = XYA(j,:) + ellA(j)*[cos(psiA(j)), sin(psiA(j))];
end

XYC = zeros(N+1,2);
XYC(1,:) = p.C_flank(:).';
ellC = a0vec + eC;                         % current lengths 
psiC = zeros(N,1);                        % absolute segment angles
psiC(1) = p.alphaL_C;
psiC(2:end) = p.alphaL_C + cumsum(phiC);
for j = 1:N
    XYC(j+1,:) = XYC(j,:) + ellC(j)*[cos(psiC(j)), sin(psiC(j))];
end

E_one = fval;
pack = struct( ...
    'phiB',  phiB, ...
    'eB',    eB, ...
    'phiA',  phiA, ...
    'eA',    eA, ...
    'phiC',  phiC, ...
    'eC',    eC, ...
    'theta', theta, ...
    'B',     B, ...
    'A',     A, ...
    'C',     C,...
    'XYB',   XYB,...
    'XYA',   XYA,...
    'XYC',   XYC,...
    'xG',    xG,...
    'yG',    yG);

    function href = href_end_cluster(L0, N, r, p)
        % End-clustered segment lengths that sum to L0
        if nargin < 4, p = 2.0; end
        r = max(1e-6, min(0.9999, r));
        s = ((1:N)' - 0.5)/N;
        d = min(s, 1 - s);
        w = ((1 - r) + d).^p;
        href = (w / sum(w)) * L0;
    end
end
