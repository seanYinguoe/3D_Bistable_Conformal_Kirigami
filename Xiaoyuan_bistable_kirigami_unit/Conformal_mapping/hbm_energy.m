function [Etotal, XY, springs, err, info] = hbm_energy(L0, A1,B1, vecA1,vecB1, N, E, b, t, varargin)
% Newton–KKT (Augmented Lagrangian) solver for Hencky bar-chain (bending + axial)
% Robust: box bounds on kappa/ell, fraction-to-boundary step, backtracking, AL multiplier update.
%
% Inputs:
%   L0          : undeformed length (reference) of the centerline
%   A1,B1       : end-point positions in the deformed configuration (2x1 each)
%   vecA1,vecB1 : end directions (normals if 'VectorsAreNormals'==true; otherwise tangents)
%   N           : number of links (>=2), K = N-1 hinges
%   E,b,t       : material & section (EI = E*b*t^3/12, EA = E*b*t)
%
% Name-Value options:
%   'VectorsAreNormals' (true) : if true, convert normals to tangents along A1->B1
%   'MaxIter' (1000)           : max Newton iterations
%   'TolC' (1e-10)             : constraint norm tolerance for convergence
%   'TolStep' (1e-12)          : step norm tolerance
%   'Damping' (1.0)            : base damping (0< Damping <= 1)
%   'KappaMax' (pi)            : bound |kappa_i| <= KappaMax
%   'Stretch' (0.4)            : bound ell_i in [a0*(1-Stretch), a0*(1+Stretch)]
%   'AngScale' (10)            : scale for angle residual (sin(Δ) term)
%   'Penalty' (10)             : AL penalty rho (will adapt)
%
% Outputs:
%   Etotal      : total energy (Eb + Ex)
%   XY          : (N+1)x2 node coordinates of the deformed centerline
%   springs     : struct with plotting helpers
%   err         : final scaled constraint norm
%   info        : diagnostics

% ---------- options ----------
p = inputParser;
addParameter(p,'VectorsAreNormals',true,@islogical);
addParameter(p,'MaxIter',1000,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TolC',1e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TolStep',1e-12,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Damping',1.0,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
addParameter(p,'KappaMax',pi,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Stretch',0.4,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'AngScale',10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Penalty',10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(p,varargin{:});
opt = p.Results;

% ---------- small helpers ----------
wrap  = @(th) atan2(sin(th),cos(th));      % wrap to [-pi,pi]
v2ang = @(v) atan2(v(2),v(1));             % angle of 2D vector
nrm1  = @(v) v./max(norm(v),eps);          % normalize

% ---------- convert normals -> tangents if required ----------
if opt.VectorsAreNormals
    vecA1 = normal_to_tangent(vecA1, A1, B1);
    vecB1 = normal_to_tangent(vecB1, A1, B1);
else
    vecA1 = nrm1(vecA1); vecB1 = nrm1(vecB1);
end
thetaA = v2ang(vecA1);  thetaB = v2ang(vecB1);

% ---------- constants ----------
if N < 2, error('N must be >= 2.'); end
a0  = L0/N;
EI  = E*b*t^3/12;
EA  = E*b*t;

K  = N-1;                      % number of hinges
kr = (EI/a0);                  % rotational "spring" per hinge (EI/a0)
ka = (EA/a0);                  % axial "spring" per segment (EA/a0)

% Hessian of the (pure) energy is constant diagonal
H  = diag([kr*ones(K,1); ka*ones(N,1)]);

% ---------- initial guess ----------
kappa = (wrap(thetaB-thetaA)/K)*ones(K,1);
ell   = a0*ones(N,1);
z     = [kappa; ell];

% ---------- box bounds (上下界) ----------
stretch = opt.Stretch;
e0 = a0*ones(N,1);
lb = [ -opt.KappaMax*ones(K,1);  e0*(1 - stretch) ];
ub = [  opt.KappaMax*ones(K,1);  e0*(1 + stretch) ];
epsb = 1e-12;
z = min(max(z, lb+epsb), ub-epsb);   % ensure feasible start

% ---------- AL multipliers & penalty ----------
m      = 3;                            % constraints: gx, gy, gtheta(sin)
lambda = zeros(m,1);
rho    = opt.Penalty;

% ---------- iteration ----------
step_norm = NaN;  c_norm = NaN;
for it = 1:opt.MaxIter

    % unpack
    kappa = z(1:K);
    ell   = z(K+1:end);

    % geometry: angles, tangents, chain coordinates
    phi = thetaA + [0; cumsum(kappa(:))];      % N x 1
    txy = [cos(phi), sin(phi)];                % N x 2
    XY  = forward_chain(A1, phi, ell);

    % constraints g(z) = [gx; gy; gtheta] (使用 sin(Δ) 更平滑)
    pos_err = (A1(:) + sum(ell(:).*txy,1).' - B1(:));   % 2x1
    Delta   = phi(end) - thetaB;                        % scalar
    s_ang   = opt.AngScale;                             % angle scale
    g       = [ pos_err/a0;  s_ang * sin(Delta) ];      % 3x1

    % gradient of energy (pure E): gradE = [kr*kappa; ka*(ell-a0)]
    gradE = [kr*kappa; ka*(ell - a0)];

    % Jacobian A = dg/dz
    % dPos/dell = t^T
    dPos_dell = txy.';                    % 2 x N
    % dPos/dkappa
    dPos_dk = zeros(2,K);
    for k = 1:K
        idx = (k+1):N;
        if ~isempty(idx)
            dPos_dk(:,k) = [ -sum(ell(idx).*sin(phi(idx)));  sum(ell(idx).*cos(phi(idx))) ];
        end
    end
    % angle row: d/dkappa sin(Delta) = cos(Delta)*1, d/dell = 0
    cD       = cos(Delta);
    dAng_dk   = s_ang * cD * ones(1,K);
    dAng_dell = zeros(1,N);

    A = [ dPos_dk/a0,  dPos_dell/a0;      % 2 x (K+N)
          dAng_dk,     dAng_dell ];       % 1 x (K+N)  -> total 3 rows

    % ----- Augmented Lagrangian KKT system -----
    % Lagrangian gradient:
    gradL = gradE + A.'*lambda;
    % Stabilized left block (positive-definite if rho>0):
    KKT  = [ H + rho*(A.'*A),  A.' ;
             A,                zeros(m,m) ];
    rhs  = -[ gradL + rho*A.'*g ;
              g ];

    sol  = KKT \ rhs;                      % prefer sparse if large
    dz   = sol(1:K+N);
    dl   = sol(K+N+1:end);

    % ----- fraction-to-the-boundary step (stay inside lb/ub) -----
    alpha_bd = 1.0;
    pos = dz > 0;
    if any(pos)
        alpha_bd = min(alpha_bd, min( (ub(pos) - z(pos))./dz(pos) ));
    end
    neg = dz < 0;
    if any(neg)
        alpha_bd = min(alpha_bd, min( (lb(neg) - z(neg))./dz(neg) ));
    end
    if ~isfinite(alpha_bd), alpha_bd = 1.0; end
    alpha_bd = max(0, 0.99*alpha_bd);

    % ----- backtracking on AL merit  φ(z) = E(z) + (rho/2)||g(z)||^2 -----
    alpha   = min(opt.Damping, alpha_bd);
    [E0, g0] = eval_state(z);               % current merit pieces
    phi0 = E0 + 0.5*rho*(g0.'*g0);
    c1   = 1e-4;

    while true
        z_try = z + alpha*dz;
        z_try = min(max(z_try, lb+epsb), ub-epsb);         % keep inside bounds
        [E1, g1] = eval_state(z_try);
        phi1 = E1 + 0.5*rho*(g1.'*g1);

        % simple decrease check (Armijo-like)
        % use quadratic model along dz: dz'*(H+rho A'A)*dz
        suff_dec = phi1 <= phi0 - c1*alpha*(dz.'*(H + rho*(A.'*A))*dz);
        if suff_dec || alpha < 1e-8, break; end
        alpha = 0.5*alpha;
    end

    % accept
    z = z_try;

    % multiplier update (Augmented Lagrangian)
    % Option A (robust): lambda = lambda + rho * g_new
    lambda = lambda + rho * g1;
    % Option B (quasi-Newton multipliers): lambda = lambda + dl;

    % diagnostics at accepted step
    step_norm = norm(alpha*dz);
    c_norm    = norm(g1);

    % adapt rho if constraints stagnate
    if norm(g1) > 0.7*norm(g0)
        rho = min(1e8, 5*rho);
    end

    % stopping
    if step_norm < opt.TolStep && c_norm < opt.TolC
        break;
    end
end

% ---------- final geometry and energy ----------
kappa = z(1:K); ell = z(K+1:end);
phi = thetaA + [0; cumsum(kappa(:))];
txy = [cos(phi), sin(phi)];
XY  = forward_chain(A1, phi, ell);
Eb = 0.5*kr * sum(kappa.^2);
Ex = 0.5*ka * sum((ell - a0).^2);
Etotal = Eb + Ex;

springs.axial      = [XY(1:end-1,1) XY(1:end-1,2) XY(2:end,1) XY(2:end,2)];
springs.rotational = XY(2:end-1,:);
% final constraint (scaled) using sin residual
pos_err = (A1(:) + sum(ell(:).*txy,1).' - B1(:));
Delta   = phi(end) - thetaB;
g_final = [ pos_err/a0;  opt.AngScale * sin(Delta) ];
err = norm(g_final);

% info
info.iters    = it;
info.kr       = kr;
info.ka       = ka;
info.a0       = a0;
info.kappa    = kappa;
info.ell      = ell;
info.lambda   = lambda;
info.rho      = rho;
info.lb       = lb; 
info.ub       = ub;
info.constr   = g_final;

fprintf('KKT-AL+Bnds: it=%d, ||g||=%.3e, step=%.3e, E=%.6g (Eb=%.6g, Ex=%.6g), rho=%g\n',...
    it, err, step_norm, Etotal, Eb, Ex, rho);

% ---------- nested: evaluate E(z) & g(z) ----------
    function [E, gvec] = eval_state(zvec)
        k = zvec(1:K); L = zvec(K+1:end);
        ph = thetaA + [0; cumsum(k(:))];
        tx = [cos(ph), sin(ph)];
        pos_loc = (A1(:) + sum(L(:).*tx,1).' - B1(:));
        dlt     = ph(end) - thetaB;
        gvec    = [ pos_loc/a0;  opt.AngScale * sin(dlt) ];
        Eb_ = 0.5*kr*sum(k.^2);
        Ex_ = 0.5*ka*sum((L - a0).^2);
        E   = Eb_ + Ex_;
    end
end

% ---------- helpers ----------
function XY = forward_chain(A, phi, ell)
N  = numel(phi);
XY = zeros(N+1,2);
XY(1,:) = A;
for i=1:N
    XY(i+1,:) = XY(i,:) + ell(i)*[cos(phi(i)), sin(phi(i))];
end
end

function v_tan = normal_to_tangent(v_normal, A, B)
R90  = [0 -1; 1 0];                  % +90° CCW
v1   = (R90 * v_normal(:)).';
v2   = -v1;
pref = (B - A);
if dot(v1, pref) >= dot(v2, pref), v_tan = v1; else, v_tan = v2; end
n = norm(v_tan); if n>0, v_tan = v_tan / n; end
end
