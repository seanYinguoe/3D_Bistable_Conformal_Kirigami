function [Etotal, X, Rlist, springs, err, info] = hbm_energy_3d(L0, A1, B1, RA, RB, N, E, b, t, varargin)
% 3D Hencky bar-chain (bending + axial) with AL-KKT and exact end pose via constraints.
% Inputs:
%   A1,B1 : 3x1 vectors
%   RA,RB : 3x3 rotations; end-frame at A and target end-frame at B 
%   N     : number of links; K = N-1 hinges
%   E,b,t : EI = E*b*t^3/12, EA = E*b*t
%
% Options (Name-Value):
%   'MaxIter'(1000), 'TolC'(1e-10), 'TolStep'(1e-12), 'Damping'(0.9),
%   'KappaMax'(2*pi), 'Stretch'(0.6), 'Penalty'(1e7),
%   'EndCluster'(true), 'ClusterRatio'(0.85), 'DoSOC'(true), 'DoLenProjection'(false),
%   'z0'([]), 'FDStep'(1e-8)
%
% Outputs:
%   X      : (N+1) x 3 node coordinates
%   Rlist  : 1..N rotation matrices along chain
%   springs: geometry for plotting
%   err    : ||position error|| + ||orientation error||
%   info   : diagnostics

% ---- options ----
p = inputParser;
addParameter(p,'MaxIter',1000);
addParameter(p,'TolC',1e-10);
addParameter(p,'TolStep',1e-12);
addParameter(p,'Damping',0.9);
addParameter(p,'KappaMax',2*pi);
addParameter(p,'Stretch',0.6);
addParameter(p,'Penalty',1e7);
addParameter(p,'EndCluster',true);
addParameter(p,'ClusterRatio',0.85);
addParameter(p,'DoSOC',true);
addParameter(p,'DoLenProjection',false);
addParameter(p,'FDStep',1e-8);
addParameter(p,'z0',[]);
parse(p,varargin{:});
opt = p.Results;

% ---- constants ----
if N < 2, error('N must be >= 2.'); end
K  = N-1;
EI = E*b*t^3/12;
EA = E*b*t;

% reference segmentation with end clustering
if opt.EndCluster
    a0 = href_end_cluster_1D(L0, N, opt.ClusterRatio);
else
    a0 = (L0/N)*ones(N,1);
end
hhinge = 0.5*(a0(1:end-1)+a0(2:end));           % Kx1
kr     = (EI ./ hhinge);                         % bending coeff per hinge (scalar here)
ka     = (EA ./ a0);                             % axial per segment

% decision vector: [omega(3*K); ell(N)]
if isempty(opt.z0)
    omega = zeros(3*K,1);                        % small-curvature start
    ell   = a0;                                   % start at ref lengths
    z     = [omega; ell];
else
    z = opt.z0(:);
    assert(numel(z)==(3*K+N),'z0 size mismatch');
end

% bounds
lb = [-opt.KappaMax*ones(3*K,1); a0.*(1-opt.Stretch)];
ub = [ opt.KappaMax*ones(3*K,1); a0.*(1+opt.Stretch)];
epsb = 1e-12; z = min(max(z, lb+epsb), ub-epsb);

% AL setup
m = 6;                          % 3 pos + 3 rot
lambda = zeros(m,1);
rho    = opt.Penalty;
ascale = max(L0, norm(B1 - A1));

% helpers
hat  = @(v)[  0   -v(3)  v(2); v(3)   0   -v(1); -v(2) v(1)  0];
vee  = @(M)[M(3,2); M(1,3); M(2,1)];
wrapR = @(R) project_SO3(R);    % small projection to SO(3) for numeric stability

% main loop
step_norm = NaN;  c_norm = NaN;
for it = 1:opt.MaxIter
    [E, g, X, Rlist, gradE, H, A] = eval_all(z);
    gradL = gradE + A.'*lambda;
    KKT   = [H + rho*(A.'*A),  A.';  A,  zeros(m,m)];
    rhs   = -[ gradL + rho*A.'*g ; g ];
    % Newton step
    sol   = KKT \ rhs;
    dz    = sol(1:numel(z));

    % fraction-to-boundary
    alpha_bd = 1.0;
    posmask = dz > 0; if any(posmask), alpha_bd = min(alpha_bd, min( (ub(posmask)-z(posmask))./dz(posmask) )); end
    negmask = dz < 0; if any(negmask), alpha_bd = min(alpha_bd, min( (lb(negmask)-z(negmask))./dz(negmask) )); end
    if ~isfinite(alpha_bd), alpha_bd = 1.0; end
    alpha_bd = max(0, 0.99*alpha_bd);

    % backtracking on AL merit
    alpha = min(opt.Damping, alpha_bd);
    phi0  = E + 0.5*rho*(g.'*g);
    c1    = 1e-4;  Mmat = (H + rho*(A.'*A));
    while true
        z_try = min(max(z + alpha*dz, lb+epsb), ub-epsb);
        [E1, g1] = eval_Eg(z_try);
        phi1 = E1 + 0.5*rho*(g1.'*g1);
        suff = (phi1 <= phi0 - c1*alpha*(dz.'*Mmat*dz)) || (alpha < 1e-12);
        if suff, break; end
        alpha = 0.5*alpha;
    end

    % accept
    z = z_try; g = g1;

    % penalty growth if stalled
    if norm(g1) > 0.25*norm(eval_Eg(z - alpha*dz))  % simple check
        rho = min(1e12, 10*rho);
    end

    % optional SOC on constraints
    if opt.DoSOC
        [~, g_soc, ~, ~, ~, ~, A_soc] = eval_all(z);
        % minimal-norm d: min ||d|| s.t. A d = -g
        Gsoc = (A_soc*A_soc.');
        dlam = Gsoc \ (-g_soc);
        dmin = A_soc.' * dlam;
        % bound-safe
        alpha2 = 1.0;
        pos2 = dmin > 0; if any(pos2), alpha2 = min(alpha2, min((ub(pos2)-z(pos2))./dmin(pos2))); end
        neg2 = dmin < 0; if any(neg2), alpha2 = min(alpha2, min((lb(neg2)-z(neg2))./dmin(neg2))); end
        if ~isfinite(alpha2), alpha2 = 1.0; end
        z = z + 0.9*0.99*alpha2*dmin;
    end

    % optional length-only feasibility projection
    if opt.DoLenProjection
        [~, g_now, X_now, R_now] = eval_Eg(z);
        Jlen = dpos_dell(X_now, R_now);           % 3xN
        pos_err = g_now(1:3)*ascale;              % unscale
        JJt = Jlen*Jlen.'; reg = 1e-15*trace(JJt);
        dlam = (JJt + reg*eye(3)) \ (-pos_err);
        dell = Jlen.' * dlam;                     % N x 1
        % bound-safe on ell only
        ell = z(3*K+1:end); ubL = ub(3*K+1:end); lbL = lb(3*K+1:end);
        alpha3 = 1.0;
        p = dell>0; if any(p), alpha3 = min(alpha3, min((ubL(p)-ell(p))./dell(p))); end
        n = dell<0; if any(n), alpha3 = min(alpha3, min((lbL(n)-ell(n))./dell(n))); end
        if ~isfinite(alpha3), alpha3 = 1.0; end
        ell = ell + 0.95*0.99*alpha3*dell;
        z(3*K+1:end) = ell;
    end

    % stopping
    step_norm = norm(alpha*dz);
    [~, gfin] = eval_Eg(z);
    c_norm = norm(gfin);
    if step_norm < opt.TolStep && c_norm < opt.TolC, break; end
end

% final state
[Etotal, gfin, X, Rlist] = eval_Eg(z);
err = norm(gfin);
omega = reshape(z(1:3*K),3,[]);
ell   = z(3*K+1:end);
Eb = 0.5*sum( sum( (sqrt(kr(:)).'.*omega).^2 ,1) );
Ex = 0.5*sum( ka(:).*(ell(:)-a0(:)).^2 );

springs.axial = [X(1:end-1,:), X(2:end,:)];  % for plotting
info = struct('iters',it,'rho',rho,'kappa',omega,'ell',ell,'a0',a0,'Eb',Eb,'Ex',Ex);

fprintf('HBC-3D AL: it=%d, ||g||=%.3e, step=%.3e, E=%.6g (Eb=%.6g, Ex=%.6g), rho=%g\n',...
    it, err, step_norm, Etotal, Eb, Ex, rho);

% ---------- nested: energy, constraints, jacobians ----------
    function [E, g, X, Rlst, gradE, H, A] = eval_all(zv)
        [E, g, X, Rlst] = eval_Eg(zv);

        % energy gradient/Hessian (bending + axial, block diagonal, simple)
        om = reshape(zv(1:3*K),3,[]);
        el = zv(3*K+1:end);

        % bending
        gradEb_omega = zeros(3*K,1);
        Hbo = zeros(3*K,3*K);
        for i=1:K
            Di = kr(i)*eye(3);
            oi = om(:,i);
            gradEb_omega(3*(i-1)+(1:3)) = Di*oi;
            Hbo(3*(i-1)+(1:3), 3*(i-1)+(1:3)) = Di;
        end
        % axial
        gradEx_l = ka(:).*(el(:)-a0(:));
        Hl = diag(ka(:));

        gradE = [gradEb_omega; gradEx_l];
        H = blkdiag(Hbo, Hl);

        % constraint Jacobian via finite difference (central)
        if nargout >= 7
            A = jacobian_fd(@eval_g_only, zv, opt.FDStep);
        end
    end

    function [E, g, X, Rlst] = eval_Eg(zv)
        % unpack
        om = reshape(zv(1:3*K),3,[]);
        el = zv(3*K+1:end);

        % forward chain (positions & rotations)
        X = zeros(N+1,3); Rlst = cell(N,1);
        X(1,:) = A1(:).';
        R = RA;  % start frame at A
        for i=1:N
            if i<=K
                R = wrapR(R * expm(hat(om(:,i))));
            end
            Rlst{i} = R;
            X(i+1,:) = X(i,:) + (R*[1;0;0]*el(i)).';
        end

        % energy
        Eb = 0.0;
        for i=1:K
            Eb = Eb + 0.5*om(:,i).'* (kr(i)*eye(3)) * om(:,i);
        end
        Ex = 0.5*sum( ka(:).*(el(:)-a0(:)).^2 );
        E  = Eb + Ex;

        % constraints: pos + rot
        pos_err = (X(end,:).' - B1(:));
        R_end   = Rlst{end};
        R_err   = R_end.' * RB;              % want R_end == RB  => log(R_endᵀ RB) = 0
        rot_err = vee( logm_clamped(R_err) );
        g = [pos_err; rot_err] / ascale;
    end

    function g = eval_g_only(zv)
        [~, g] = eval_Eg(zv);
    end
end

% ---- helpers ----
function a0 = href_end_cluster_1D(L0,N,r)
    if nargin<3, r=0.85; end
    r = max(1e-6, min(0.9999,r));
    s = ((1:N)' - 0.5)/N;
    d = min(s,1-s);
    w = (1-r + d).^2;
    a0 = (w/sum(w))*L0;
end

function A = jacobian_fd(fun, z, h)
    % central FD Jacobian of g(z) (m x n)
    g0 = fun(z);
    m  = numel(g0); n = numel(z);
    A  = zeros(m,n);
    for j=1:n
        zj = z; zj(j) = zj(j)+h;
        zjm = z; zjm(j) = zjm(j)-h;
        gp = fun(zj); gm = fun(zjm);
        A(:,j) = (gp - gm)/(2*h);
    end
end

function R = project_SO3(R)
    % project to nearest SO(3) (polar decomposition)
    [U,~,V] = svd(R);
    R = U*diag([1,1,det(U*V')])*V';
end

function L = logm_clamped(R)
    % numerically stable log for rotation (skew-symmetric)
    th = acos( max(-1,min(1,(trace(R)-1)/2)) );
    if th < 1e-9
        L = 0.5*(R - R');  % small-angle
    else
        L = th/(2*sin(th))*(R - R');
    end
end

function J = dpos_dell(X, Rlist)
    % d( end position ) / d ell  (3xN)
    N = size(X,1)-1;
    J = zeros(3,N);
    for i=1:N
        e1i = Rlist{i}(:,1);
        J(:,i) = e1i;
        for k=i+1:N
            % position of end depends only through subsequent rotations (ignored here)
            % but for d pos / d ell, only direct contribution is e1_i (OK)
        end
    end
end
