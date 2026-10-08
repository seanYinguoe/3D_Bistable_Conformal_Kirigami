function [Etotal, XY, springs, err, info] = hbm_energy(L0, A1,B1, vecA1,vecB1, N, E, b, t, varargin)
% Hencky bar-chain (bending + axial) with end clustering and EXACT right-end angle.
% - Augmented Lagrangian with bounds.
% - Position constraints are tightened via:
%     * stronger penalty schedule,
%     * Second-Order Correction (SOC),
%     * optional feasibility projection on lengths only (minimal-norm).
% - Consistent discrete energy: bending Kb_i = EI / hhinge_i, axial Ka_i = EA / h_i.
%
% Inputs:
%   L0          : undeformed centerline length
%   A1,B1       : end positions in the deformed configuration (2x1 each)
%   vecA1,vecB1 : end directions (normals if 'VectorsAreNormals'==true; else tangents)
%   N           : number of links (>=2). number of hinges K = N-1
%   E,b,t       : material & section (EI = E*b*t^3/12, EA = E*b*t)
%
% Name-Value options:
%   'VectorsAreNormals' (true)
%   'MaxIter' (1000)
%   'TolC' (1e-14)          : target norm for end-position constraint (unscaled)
%   'TolStep' (1e-12)
%   'Damping' (0.9)         : base step damping
%   'KappaMax' (2*pi)
%   'Stretch' (0.6)
%   'Penalty' (1e10)        : initial AL penalty
%   'EndCluster' (true)
%   'ClusterRatio' (0.85)   : smaller -> stronger end clustering (both ends denser)
%   'DoSOC' (true)          : apply Second-Order Correction each accepted step
%   'DoLenProjection' (true): apply minimal-norm feasibility projection on lengths


%% Define parameters
% options
p = inputParser;
addParameter(p,'VectorsAreNormals',true);
addParameter(p,'MaxIter',1000);
addParameter(p,'TolC',1e-10);
addParameter(p,'TolStep',1e-12);
addParameter(p,'Damping',0.9);
addParameter(p,'KappaMax',2*pi);
addParameter(p,'Stretch',0.6);
addParameter(p,'Penalty',1e8);
addParameter(p,'EndCluster',true);
addParameter(p,'ClusterRatio',0.85); % 0.7 - 0.9
addParameter(p,'DoSOC',true);
addParameter(p,'DoLenProjection',false);
addParameter(p,'z0',[]);
parse(p,varargin{:});
opt = p.Results;

% helpers
wrap  = @(th) atan2(sin(th),cos(th));
v2ang = @(v) atan2(v(2),v(1));
nrm1  = @(v) v./max(norm(v),eps);

% convert normals -> tangents if requested
if opt.VectorsAreNormals
    vecA1 = normal_to_tangent(vecA1, A1, B1);
    vecB1 = normal_to_tangent(vecB1, A1, B1);
else
    vecA1 = nrm1(vecA1); vecB1 = nrm1(vecB1);
end
thetaA = v2ang(vecA1);  thetaB = v2ang(vecB1);

% constants
if N < 2, error('N must be >= 2.'); end
K  = N-1;                 % total hinges
Kv = K-1;                 % free hinges after eliminating the last one
EI = E*b*t^3/12;
EA = E*b*t;

% reference segmentation
if opt.EndCluster % end clustering
    a0vec = href_end_cluster(L0, N, opt.ClusterRatio); % N x 1; ends shortest
else % distrubuted evenly
    a0vec = (L0/N)*ones(N,1);
end
hhinge = 0.5*(a0vec(1:end-1) + a0vec(2:end));          % K x 1

% local stiffness
kr_vec = EI ./ hhinge;          % K x 1 (bending)
ka_vec = EA ./ a0vec;           % N x 1 (axial)
krK    = kr_vec(end);           % last hinge stiffness (eliminated one)
Dk     = diag(kr_vec(1:end-1)); % Kv x Kv

% initial guess
DeltaTot = wrap(thetaB - thetaA);     % total rotation needed
if ~isempty(opt.z0) && numel(opt.z0) == (Kv + N)
    z = opt.z0(:);                    % warm-start: [kfree; ell]
else
    kfree = (DeltaTot / K) * ones(Kv,1);   % distribute over free hinges
    ell   = a0vec;                          % start at reference lengths
    z     = [kfree; ell];
end

% upper and lower bounds
stretch = opt.Stretch;
lb = [ -opt.KappaMax*ones(Kv,1);  a0vec.*(1 - stretch) ];
ub = [  opt.KappaMax*ones(Kv,1);  a0vec.*(1 + stretch) ];
epsb = 1e-12;
z = min(max(z, lb+epsb), ub-epsb);

% AL multipliers & penalty
m      = 2;                          % gx, gy (angle is exact now)
lambda = zeros(m,1);
rho    = opt.Penalty;
ascale = max(L0, norm(B1 - A1));

%% optimisation
step_norm = NaN;  c_norm = NaN;
for it = 1:opt.MaxIter
    % unpack
    kfree = z(1:Kv);
    ell   = z(Kv+1:end);
    kK    = DeltaTot - sum(kfree);   % exact last hinge
    kappa = [kfree; kK];             % K x 1

    % geometry
    phi = thetaA + [0; cumsum(kappa)];     % N x 1
    txy = [cos(phi), sin(phi)];
    XY  = forward_chain(A1, phi, ell);

    % constraints(angle has already been considered)
    pos_err = (A1(:) + sum(ell(:).*txy,1).' - B1(:));   % 2x1
    g = pos_err/ascale;                                  % scaled 2x1

    % energy and gradient (with elimination coupling)
    onesv = ones(Kv,1);
    gradEb_k = Dk*kfree - krK*kK*onesv;
    Hk = Dk + krK*(onesv*onesv.');
    % axial part
    gradEx_l = ka_vec .* (ell - a0vec);
    Hl = diag(ka_vec);
    % total gradient and Hessian (block diagonal)
    gradE = [gradEb_k; gradEx_l];
    H     = blkdiag(Hk, Hl);

    % Jacobian A = dg/dz
    % First build dPos/dkappa for all K hinges (like standard HBC)
    dPos_dk_full = zeros(2,K);
    for k = 1:K
        idx = (k+1):N;
        if ~isempty(idx)
            dPos_dk_full(:,k) = [ -sum(ell(idx).*sin(phi(idx)));  sum(ell(idx).*cos(phi(idx))) ];
        end
    end
    % Because kK = DeltaTot - sum(kfree), dkK/dkfree = -1.
    % Reduced Jacobian wrt kfree: dPos/dkfree = dPos/dk(1:Kv) + dPos/dkK * (-1)
    dPos_dk = dPos_dk_full(:,1:Kv) - dPos_dk_full(:,end);
    % dPos/dell = t^T
    dPos_dell = txy.';                    % 2 x N

    A = [ dPos_dk/ascale,  dPos_dell/ascale ];   % 2 x (Kv+N)

    % Augmented Lagrangian KKT system
    gradL = gradE + A.'*lambda;
    KKT = [ H + rho*(A.'*A),  A.' ;
        A,                zeros(m,m) ];
    rhs = -[ gradL + rho*A.'*g ; g ];

    sol = KKT \ rhs;  % check the matrix singularity
    dz  = sol(1:Kv+N);

    % fraction-to-the-boundary step
    alpha_bd = 1.0;
    posmask = dz > 0;
    if any(posmask), alpha_bd = min(alpha_bd, min( (ub(posmask) - z(posmask))./dz(posmask) )); end
    negmask = dz < 0;
    if any(negmask), alpha_bd = min(alpha_bd, min( (lb(negmask) - z(negmask))./dz(negmask) )); end
    if ~isfinite(alpha_bd), alpha_bd = 1.0; end
    alpha_bd = max(0, 0.99*alpha_bd);

    % backtracking on AL merit
    alpha   = min(opt.Damping, alpha_bd);
    [E0, g0] = eval_state(z);
    phi0 = E0 + 0.5*rho*(g0.'*g0);
    c1   = 1e-4;

    while true
        z_try = z + alpha*dz;
        z_try = min(max(z_try, lb+epsb), ub-epsb);
        [E1, g1] = eval_state(z_try);
        phi1 = E1 + 0.5*rho*(g1.'*g1);
        suff_dec = (phi1 <= phi0 - c1*alpha*(dz.'*(H + rho*(A.'*A))*dz)) || (alpha < 1e-10);
        if suff_dec, break; end
        alpha = 0.5*alpha;
    end

    % accept
    z = z_try; g = g1; % update to accepted state

    % ---- stronger penalty growth if constraints stall ----
    if norm(g1) > 0.25*norm(g0)
        rho = min(1e12, 10*rho);
    end

    % === Second-Order Correction (SOC) to drive constraints tighter ===
    if opt.DoSOC
        % Recompute Jacobian A and residual g at the accepted point
        kfree = z(1:Kv); ell = z(Kv+1:end);
        kK    = DeltaTot - sum(kfree);
        kappa = [kfree; kK];
        phi   = thetaA + [0; cumsum(kappa)];
        txy   = [cos(phi), sin(phi)];
        pos_err = (A1(:) + sum(ell(:).*txy,1).' - B1(:));
        g_soc   = pos_err/ascale;

        dPos_dk_full = zeros(2,K);
        for k = 1:K
            idx = (k+1):N;
            if ~isempty(idx)
                dPos_dk_full(:,k) = [ -sum(ell(idx).*sin(phi(idx)));  sum(ell(idx).*cos(phi(idx))) ];
            end
        end
        dPos_dk   = dPos_dk_full(:,1:Kv) - dPos_dk_full(:,end);
        dPos_dell = txy.';
        A_soc     = [ dPos_dk/ascale,  dPos_dell/ascale ];  % 2 x (Kv+N)

        % Minimal-norm correction: solve min ||d|| s.t. A_soc d = -g_soc
        G   = (A_soc*A_soc.');     % 2x2
        rhs_soc = -g_soc;
        dlam = G \ rhs_soc;        % multipliers of LS problem
        d_soc = A_soc.' * dlam;    % minimal-norm step

        % Bound-safe SOC step
        alpha_bd2 = 1.0;
        pos2 = d_soc > 0;
        if any(pos2), alpha_bd2 = min(alpha_bd2, min( (ub(pos2) - z(pos2))./d_soc(pos2) )); end
        neg2 = d_soc < 0;
        if any(neg2), alpha_bd2 = min(alpha_bd2, min( (lb(neg2) - z(neg2))./d_soc(neg2) )); end
        if ~isfinite(alpha_bd2), alpha_bd2 = 1.0; end
        alpha_bd2 = max(0, 0.99*alpha_bd2);

        z = z + 0.9*alpha_bd2 * d_soc;
    end

    % Optional feasibility projection on lengths only (minimal-norm)
    if opt.DoLenProjection
        kfree = z(1:Kv); ell = z(Kv+1:end);
        kK    = DeltaTot - sum(kfree);
        phi   = thetaA + [0; cumsum([kfree; kK])];
        txy   = [cos(phi), sin(phi)];
        pos_err = (A1(:) + sum(ell(:).*txy,1).' - B1(:));   % unscaled

        if norm(pos_err) > opt.TolC
            J = txy.';   % 2 x N = d pos / d ell

            % Minimal-norm correction: solve min ||dell|| s.t. J*dell = -pos_err
            % Stabilize the 2x2 normal equations with a tiny Tikhonov term.
            JJt = J*J.';                  % 2x2
            reg = 1e-15 * trace(JJt);
            dlam = (JJt + reg*eye(2)) \ (-pos_err);

            dell = J.' * dlam;            % N x 1
            % Make sure everything is column vectors
            ell  = ell(:);  dell = dell(:);

            % Bound-safe step length (fraction-to-the-boundary) on lengths only
            ubL = ub(Kv+1:end);  ubL = ubL(:);
            lbL = lb(Kv+1:end);  lbL = lbL(:);

            alpha_bd3 = 1.0;

            idxP = dell > 0;                  % moving up -> check upper bound
            if any(idxP)
                num = ubL(idxP) - ell(idxP);
                den = dell(idxP);
                % avoid division by ~0
                den = max(den, 1e-16);
                alpha_bd3 = min(alpha_bd3, min(num ./ den));
            end

            idxN = dell < 0;                  % moving down -> check lower bound
            if any(idxN)
                num = lbL(idxN) - ell(idxN);
                den = dell(idxN);
                den = min(den, -1e-16);
                alpha_bd3 = min(alpha_bd3, min(num ./ den));
            end

            if ~isfinite(alpha_bd3), alpha_bd3 = 1.0; end
            alpha_bd3 = max(0, 0.99*alpha_bd3);

            % apply a conservative projection
            ell = ell + 0.95*alpha_bd3 * dell;

            % write back
            z(Kv+1:end) = ell;
        end
    end

    % diagnostics for stopping
    step_norm = norm(alpha*dz);
    kfree = z(1:Kv); ell = z(Kv+1:end);
    kK    = DeltaTot - sum(kfree);
    kappa = [kfree; kK];
    phi = thetaA + [0; cumsum(kappa)];
    txy = [cos(phi), sin(phi)];
    pos_err = (A1(:) + sum(ell(:).*txy,1).' - B1(:));
    c_norm  = norm(pos_err);   % unscaled for stopping

    if step_norm < opt.TolStep && c_norm < opt.TolC
        break;
    end
end

% ---------- final geometry & energies ----------
kfree = z(1:Kv);
ell   = z(Kv+1:end);
kK    = DeltaTot - sum(kfree);
kappa = [kfree; kK];

phi = thetaA + [0; cumsum(kappa)];
txy = [cos(phi), sin(phi)];
XY  = forward_chain(A1, phi, ell);

Eb = 0.5*( kfree.'*Dk*kfree + krK*kK^2 );
Ex = 0.5*sum( ka_vec .* ((ell - a0vec).^2) );
Etotal = Eb + Ex;

springs.axial      = [XY(1:end-1,1) XY(1:end-1,2) XY(2:end,1) XY(2:end,2)];
springs.rotational = XY(2:end-1,:);

% report (both unscaled and scaled)
pos_err = (A1(:) + sum(ell(:).*txy,1).' - B1(:));
err = norm(pos_err);
g_scaled = pos_err/ascale;

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
info.constr_unscaled = pos_err;
info.constr_scaled   = g_scaled;

fprintf('HBC exact-angle + AL: it=%d, ||g_pos||=%.3e (scaled=%.3e), step=%.3e, E=%.6g (Eb=%.6g, Ex=%.6g), rho=%g\n',...
    it, err, norm(g_scaled), step_norm, Etotal, Eb, Ex, rho);

% ---------- nested: evaluate E(z) & g(z) ----------
    function [E, gvec] = eval_state(zvec)
        kf = zvec(1:Kv); L = zvec(Kv+1:end);
        kKloc = DeltaTot - sum(kf);
        kvec  = [kf; kKloc];
        ph = thetaA + [0; cumsum(kvec)];
        tx = [cos(ph), sin(ph)];
        pos_loc = (A1(:) + sum(L(:).*tx,1).' - B1(:));
        % energy
        Eb_ = 0.5*( kf.'*Dk*kf + krK*kKloc^2 );
        Ex_ = 0.5*sum( ka_vec .* ((L - a0vec).^2) );
        E   = Eb_ + Ex_;
        gvec = pos_loc/ascale;  % scaled position constraints
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
d = min(s, 1 - s);               % distance to nearest end(left or right)

% KEY: weight small at ends, large in middle
w = (eps + d).^p;

% normalize to total length L0
href = (w / sum(w)) * L0;
end
