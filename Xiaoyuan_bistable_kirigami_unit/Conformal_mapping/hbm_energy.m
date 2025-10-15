function [Etotal, XY, springs, err, info] = hbm_energy(L0, A1,B1, vecA1,vecB1, N, E, b, t, varargin)
% HENCKY BAR-CHAIN (bending + axial) with TRUE end clustering.
% - Shorter segments near both clamps (end clustering) -> denser hinges.
% - Per-hinge (bending) and per-segment (axial) stiffness use local lengths.
% - Augmented Lagrangian (AL) Newton-KKT solver with box bounds.
%
% Inputs:
%   L0          : undeformed centerline length
%   A1,B1       : end-point positions in the deformed configuration (2x1 each)
%   vecA1,vecB1 : end directions (normals if 'VectorsAreNormals'==true; else tangents)
%   N           : number of links (>=2). Hinges K = N-1
%   E,b,t       : material & section (EI = E*b*t^3/12, EA = E*b*t)
%
% Name-Value options:
%   'VectorsAreNormals' (true) : if true, convert given normals to tangents
%   'MaxIter' (1000)           : max Newton iterations
%   'TolC' (1e-10)             : constraint norm tolerance
%   'TolStep' (1e-12)          : step norm tolerance
%   'Damping' (0.6)            : base step damping (0<Damping<=1)
%   'KappaMax' (2*pi)          : bound |kappa_i| <= KappaMax
%   'Stretch' (0.6)            : per-segment stretch bound: ell_i in [h_i*(1-..), h_i*(1+..)]
%   'AngScale' (20)            : scaling for angle residual sin(Δ)
%   'Penalty' (1e11)           : initial AL penalty ρ
%   'EndCluster' (true)        : use end clustering (shorter segments near ends)
%   'ClusterRatio' (0.85)      : geometric ratio r in (0,1); smaller -> stronger end clustering
%
% Outputs:
%   Etotal      : total energy (Eb + Ex)
%   XY          : (N+1) x 2 node coordinates
%   springs     : struct with plotting helpers
%   err         : final (unscaled) constraint norm
%   info        : diagnostics and internal arrays

% ---------- options ----------
p = inputParser;
addParameter(p,'VectorsAreNormals',true);
addParameter(p,'MaxIter',1000);
addParameter(p,'TolC',1e-10);
addParameter(p,'TolStep',1e-12);
addParameter(p,'Damping',0.6);
addParameter(p,'KappaMax',2*pi);
addParameter(p,'Stretch',0.6);
addParameter(p,'AngScale',20);
addParameter(p,'Penalty',1e11);
addParameter(p,'EndCluster',true);
addParameter(p,'ClusterRatio',0.85);
parse(p,varargin{:});
opt = p.Results;

% ---------- helpers ----------
wrap  = @(th) atan2(sin(th),cos(th));  % wrap angle to [-pi,pi]
v2ang = @(v) atan2(v(2),v(1));         % angle of a 2D vector
nrm1  = @(v) v./max(norm(v),eps);      % safe normalization

% ---------- convert normals -> tangents if requested ----------
if opt.VectorsAreNormals
    vecA1 = normal_to_tangent(vecA1, A1, B1);
    vecB1 = normal_to_tangent(vecB1, A1, B1);
else
    vecA1 = nrm1(vecA1); vecB1 = nrm1(vecB1);
end
thetaA = v2ang(vecA1);  thetaB = v2ang(vecB1);

% ---------- constants ----------
if N < 2, error('N must be >= 2.'); end
K  = N-1;                 % number of hinges
EI = E*b*t^3/12;          % bending rigidity
EA = E*b*t;               % axial rigidity

% ---------- reference segmentation (TRUE end clustering) ----------
if opt.EndCluster
    a0vec = href_end_cluster(L0, N, opt.ClusterRatio); % N x 1, ends are shortest
else
    a0vec = (L0/N)*ones(N,1);
end
% representative hinge length = average of adjacent segments
hhinge = 0.5*(a0vec(1:end-1) + a0vec(2:end));          % K x 1

% ---------- local stiffness (consistent discrete energy) ----------
kr_vec = EI ./ hhinge;    % bending torsional stiffness per hinge:   Kb_i = EI / h̄_i
ka_vec = EA ./ a0vec;     % axial stiffness per segment:              Ka_i = EA / h_i

% Hessian (pure energy), diagonal
H = diag([kr_vec; ka_vec]);  % size (K+N) x (K+N)

% ---------- initial guess ----------
kappa = (wrap(thetaB-thetaA)/K)*ones(K,1);   % uniform curvature
ell   = a0vec;                               % start at reference lengths
z     = [kappa; ell];

% ---------- box bounds ----------
stretch = opt.Stretch;
lb = [ -opt.KappaMax*ones(K,1);  a0vec.*(1 - stretch) ];
ub = [  opt.KappaMax*ones(K,1);  a0vec.*(1 + stretch) ];
epsb = 1e-12;
z = min(max(z, lb+epsb), ub-epsb);           % feasible start

% ---------- AL multipliers & penalty ----------
m      = 3;                  % constraints: gx, gy, gtheta(sin)
lambda = zeros(m,1);
rho    = opt.Penalty;

ascale = L0;                 % position constraint scaling

% ---------- AL Newton-KKT iteration ----------
step_norm = NaN;  c_norm = NaN;
for it = 1:opt.MaxIter

    % unpack
    kappa = z(1:K);
    ell   = z(K+1:end);

    % geometry: angles, tangents, chain coordinates
    phi = thetaA + [0; cumsum(kappa(:))];  % N x 1 angles for segments
    txy = [cos(phi), sin(phi)];            % N x 2 tangents
    XY  = forward_chain(A1, phi, ell);     % (N+1) x 2 nodes

    % constraints g(z) = [gx; gy; s_ang*sin(Delta)]
    pos_err = (A1(:) + sum(ell(:).*txy,1).' - B1(:));   % 2x1
    Delta   = phi(end) - thetaB;                        % scalar
    s_ang   = opt.AngScale;
    g       = [ pos_err/ascale;  s_ang * sin(Delta) ];  % 3x1

    % energy gradient ∇E = [Kb.*kappa; Ka.*(ell - a0)]
    gradE = [kr_vec .* kappa;  ka_vec .* (ell - a0vec)];

    % Jacobian A = dg/dz
    % dPos/dell = t^T
    dPos_dell = txy.';                    % 2 x N
    % dPos/dkappa (accumulated downstream effect)
    dPos_dk = zeros(2,K);
    for k = 1:K
        idx = (k+1):N;
        if ~isempty(idx)
            dPos_dk(:,k) = [ -sum(ell(idx).*sin(phi(idx)));  sum(ell(idx).*cos(phi(idx))) ];
        end
    end
    % angle row
    cD        = cos(Delta);
    dAng_dk   = s_ang * cD * ones(1,K);
    dAng_dell = zeros(1,N);

    A = [ dPos_dk/ascale,  dPos_dell/ascale;    % 2 x (K+N)
          dAng_dk,         dAng_dell ];         % 1 x (K+N)

    % Augmented Lagrangian KKT system
    gradL = gradE + A.'*lambda;
    KKT  = [ H + rho*(A.'*A),  A.' ;
             A,                zeros(m,m) ];
    rhs  = -[ gradL + rho*A.'*g ; g ];

    sol  = KKT \ rhs;                 % for large N, switch to sparse
    dz   = sol(1:K+N);
    dl   = sol(K+N+1:end); %#ok<NASGU>

    % fraction-to-the-boundary step
    alpha_bd = 1.0;
    pos = dz > 0;
    if any(pos), alpha_bd = min(alpha_bd, min( (ub(pos) - z(pos))./dz(pos) )); end
    neg = dz < 0;
    if any(neg), alpha_bd = min(alpha_bd, min( (lb(neg) - z(neg))./dz(neg) )); end
    if ~isfinite(alpha_bd), alpha_bd = 1.0; end
    alpha_bd = max(0, 0.99*alpha_bd);

    % backtracking on AL merit φ(z) = E(z) + (ρ/2)||g(z)||^2
    alpha   = min(opt.Damping, alpha_bd);
    [E0, g0] = eval_state(z);
    phi0 = E0 + 0.5*rho*(g0.'*g0);
    c1   = 1e-4;

    while true
        z_try = z + alpha*dz;
        z_try = min(max(z_try, lb+epsb), ub-epsb);
        [E1, g1] = eval_state(z_try);
        phi1 = E1 + 0.5*rho*(g1.'*g1);

        suff_dec = phi1 <= phi0 - c1*alpha*(dz.'*(H + rho*(A.'*A))*dz);
        if suff_dec || alpha < 1e-8, break; end
        alpha = 0.5*alpha;
    end

    % accept
    z = z_try;

    % multiplier update
    lambda = lambda + rho * g1;

    % diagnostics
    step_norm = norm(alpha*dz);
    c_norm    = norm(g1);

    % adapt penalty if constraints stagnate
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

Eb = 0.5*sum( kr_vec .* (kappa.^2) );
Ex = 0.5*sum( ka_vec .* ((ell - a0vec).^2) );
Etotal = Eb + Ex;

springs.axial      = [XY(1:end-1,1) XY(1:end-1,2) XY(2:end,1) XY(2:end,2)];
springs.rotational = XY(2:end-1,:);

% final (unscaled) constraint
pos_err = (A1(:) + sum(ell(:).*txy,1).' - B1(:));
Delta   = phi(end) - thetaB;
g_final = [ pos_err;  opt.AngScale * sin(Delta) ];
err = norm(g_final);

% info
info.iters    = it;
info.a0vec    = a0vec;
info.hhinge   = hhinge;
info.kr_vec   = kr_vec;
info.ka_vec   = ka_vec;
info.kappa    = kappa;
info.ell      = ell;
info.lambda   = lambda;
info.rho      = rho;
info.lb       = lb;
info.ub       = ub;
info.constr   = g_final;

fprintf('HBC(AL): it=%d, ||g||=%.3e, step=%.3e, E=%.6g (Eb=%.6g, Ex=%.6g), rho=%g\n',...
    it, err, step_norm, Etotal, Eb, Ex, rho);

% ---------- nested: evaluate E(z) & g(z) ----------
    function [E, gvec] = eval_state(zvec)
        k = zvec(1:K); L = zvec(K+1:end);
        ph = thetaA + [0; cumsum(k(:))];
        tx = [cos(ph), sin(ph)];
        pos_loc = (A1(:) + sum(L(:).*tx,1).' - B1(:));
        dlt     = ph(end) - thetaB;
        gvec    = [ pos_loc/ascale;  opt.AngScale * sin(dlt) ];
        Eb_ = 0.5*sum( kr_vec .* (k.^2) );
        Ex_ = 0.5*sum( ka_vec .* ((L - a0vec).^2) );
        E   = Eb_ + Ex_;
    end
end

% ---------- helpers ----------
function XY = forward_chain(A, phi, ell)
% Accumulate chain coordinates from A using segment angles and lengths.
N  = numel(phi);
XY = zeros(N+1,2);
XY(1,:) = A(:).';
for i=1:N
    XY(i+1,:) = XY(i,:) + ell(i)*[cos(phi(i)), sin(phi(i))];
end
end

function v_tan = normal_to_tangent(v_normal, A, B)
% Convert a normal vector into a tangent consistent with A->B direction.
R90  = [0 -1; 1 0];                 % +90° CCW rotation
v1   = (R90 * v_normal(:)).';
v2   = -v1;
pref = (B(:).' - A(:).');
if dot(v1, pref) >= dot(v2, pref), v_tan = v1; else, v_tan = v2; end
n = norm(v_tan); if n>0, v_tan = v_tan / n; end
end

function href = href_end_cluster(L0, N, r, p)
% TRUE end-clustered segment lengths that sum to L0.
% r in (0,1): clustering strength (smaller r -> stronger end clustering)
% p >= 1    : exponent shaping how fast segments grow toward center
%
% Method:
%   Use a symmetric weight w(s) that is SMALL near the ends and LARGER
%   near the midspan. Segment length h_i ∝ w(s_i) so ends get shorter
%   segments => denser torsion springs near clamps.

if nargin < 4, p = 2.0; end
r = max(1e-6, min(0.9999, r));   % keep r in (0,1)
eps = (1 - r);                   % small epsilon; smaller r => smaller eps

% segment centers in (0,1)
s = ((1:N)' - 0.5)/N;
d = min(s, 1 - s);               % distance to nearest end

% KEY: weight small at ends, large in middle
w = (eps + d).^p;                % <-- NOT inverse power

% normalize to total length L0
href = (w / sum(w)) * L0;
end
