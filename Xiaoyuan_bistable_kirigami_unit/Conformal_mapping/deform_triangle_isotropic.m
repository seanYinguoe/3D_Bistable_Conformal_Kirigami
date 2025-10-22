function [Delta_out, E_total, Triangles] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta, t)
% DEFORM_TRIANGLE_ISOTROPIC (warm start, N=20, triangle_unit-based clamps)
% Inputs:
%   delta     : scalar or vector of signed edge-length changes (+ extend, - shorten)
%   edgeLen   : original equilateral outer edge length
%   l1, l4    : flank geometry parameters (your notation)
%   beta      : flank tilt angle (rad)
%   t         : filament thickness (used for visual offset of outer clamp normal only)
%
% Outputs:
%   Delta_out : the input delta as a column vector (displacements)
%   E_total   : total energy (3 ligaments) at each step
%   Triangles : 1 x nSteps cell; each is a 45x2 snapshot (your indexing)
%
% Notes:
%   - Left clamp direction:  vecA0 = (d1_ - d1), then rotate by +pi/6.
%   - Right clamp direction: vecB0 = (A   - B ), then rotate by +pi/6.
%   - These are passed as TANGENTS => 'VectorsAreNormals', True.
%   - d1 and d1_ are obtained from triangle_unit(...) (kept exactly as in your workflow).
%   - Warm start: theta uses centroid-transport; inner chain uses info.(kappa,ell) if available.

%% ----------------- helpers -----------------
rot2  = @(th) [cos(th), -sin(th); sin(th), cos(th)];
unit  = @(v) (v(:).'/max(norm(v),eps));
wrapA = @(th) atan2(sin(th), cos(th));   % angle wrap to (-pi,pi]

%% ----------------- geometry precompute (your formulas) -----------------
Aang = pi/3 - beta;
% ligament reference length (constant)
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) * sin(Aang);
% auxiliaries
l6 = (2/sqrt(3)) * sin(pi/3 - beta) * (l1 - 0.5*l4) - cos(pi/3 - beta) * l4;
l5 = ((sqrt(3)/2) * l4 + sin(beta) * l6) / sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 - l2 * ((sqrt(3)/2) * (cos(pi/3 - beta) / sin(pi/3 - beta)));
R  = (sqrt(3)/3) * l3;          % inner rigid-triangle radius

%% ----------------- discretization & material -----------------
N = 20;               % per your request
E = 4.3e11; b = 1.0;  % EI = E*b*t^3/12, EA = E*b*t (used inside hbm_energy)

%% ----------------- inputs -> vector -----------------
Delta      = delta(:);
nSteps     = numel(Delta);
Delta_out  = Delta;
E_total    = nan(nSteps,1);
Triangles  = cell(1,nSteps);

%% ----------------- initial theta (centroid -> base-center) -----------------
edge0    = edgeLen + Delta(1);
centroid = [-sqrt(3)/6*edge0, -0.5*edge0];
prev_alpha_1 = pi/3; prev_alpha_2 = 2*pi/3;
[tri0,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, Delta(1), beta, edgeLen, l1, l4, t);
d1 = tri0(20,:);
d1_ = tri0(19,:);
B0    = tri0(44,:);
theta = atan2(B0(2)-centroid(2), B0(1)-centroid(1));   % initial warm-start angle

% Optional inner-state warm start (requires hbm_energy to support 'z0')
z0 = [];

% tiny theta-solver settings
theta_tol = 1e-10;
alpha0    = 0.5;
maxIt     = 25;

% history for centroid-transport warm start
have_prev = false;
B_prev = []; C_prev = []; theta_prev2 = [];

%% ----------------- sweep over Delta -----------------
for k = 1:nSteps
    d     = Delta(k);
    edgek = edgeLen + d;                                 % signed convention (+ extend, - shorten)
    Ck    = [-sqrt(3)/6*edgek, -0.5*edgek];             % centroid from your orientation

    % -------- theta warm start with centroid transport --------
    if ~have_prev
        % keep initial theta as above for first step
    else
        v = B_prev - Ck;
        if norm(v) > 0
            theta_transport = atan2(v(2), v(1));
        else
            theta_transport = theta; % degenerate fallback
        end
        theta_extrap = theta;
        if ~isempty(theta_prev2)
            dth = wrapA(theta - theta_prev2);
            theta_extrap = wrapA(theta + dth);
        end
        theta = wrapA(0.8*theta_transport + 0.2*theta_extrap);
    end

    % -------- tiny scalar loop on theta: enforce R · t_perp = 0 --------
    theta_it = theta;
    zstart   = z0;
    for it = 1:maxIt
        [E_one, pack, ztmp] = eval_one_lig(theta_it, d1, d1_, Ck, R, ...
                                           l2, N, E, b, t, zstart);
        % reaction at B from AL multipliers
        Rvec   = reaction_at_B1(pack.info);
        t_perp = [-sin(theta_it), cos(theta_it)];
        gtheta = - R * (Rvec * t_perp.');

        if abs(gtheta) < theta_tol
            z0 = ztmp;   % accept inner-state warm start for next step
            break;
        end

        % backtracking line search along -grad(theta)
        alpha = alpha0;
        while true
            theta_try = theta_it - alpha * gtheta;
            [E_try, ~] = eval_one_lig(theta_try, d1, d1_, Ck, R, ...
                                      l2, N, E, b, t, zstart);
            if E_try <= E_one || alpha < 1e-8
                theta_it = theta_try; break;
            end
            alpha = 0.5 * alpha;
        end

        zstart = ztmp;   % refresh inner-state seed
    end

    % -------- finalize at converged theta_it --------
    [E_one, pack, z0] = eval_one_lig(theta_it, d1, d1_, Ck, R, ...
                                     l2, N, E, b, t, z0);
    A_new = pack.A; B_new = pack.B;
    C_new = Ck + (B_new - Ck) * rot2(-2*pi/3);  % third inner vertex

    % assemble snapshot with your indexing
    tri(43,:) = A_new; tri(44,:) = B_new; tri(45,:) = C_new;
    % filaments (your original pattern; thickness offset along chord)
    tri(32,:) = tri(20,:);  tri(33,:) = tri(44,:);
    tri(34,:) = tri(33,:) + t/(norm(A_new-B_new)+eps) * (A_new-B_new);
    tri(37,:) = tri(45,:);  tri(36,:) = tri(24,:);
    tri(38,:) = tri(37,:) + t/(norm(B_new-C_new)+eps) * (B_new-C_new);
    tri(40,:) = tri(28,:);  tri(41,:) = tri(43,:);
    tri(42,:) = tri(41,:) + t/(norm(C_new-A_new)+eps) * (C_new-A_new);

    Triangles{k} = tri;
    E_total(k)   = 3 * E_one;

    % update warm-start history
    theta_prev2 = have_prev * theta + (~have_prev)*[];
    theta       = theta_it;
    B_prev      = B_new;   C_prev = Ck; %#ok<NASGU>
    have_prev   = true;
end
end

%% ====================== subfunctions ======================

function [E_one, pack, z0_next] = eval_one_lig(theta, d1, d1_, Ck, R, ...
                                               l2, N, E, b, t, z0_in)
% Build inner triangle (B at angle theta), set clamp directions, call HBC once.

rot2  = @(th) [cos(th), -sin(th); sin(th), cos(th)];
unit  = @(v) (v(:).'/max(norm(v),eps));

% inner triangle vertices
B = Ck + R * [cos(theta), sin(theta)];
A = Ck + (B - Ck) * rot2( 2*pi/3);
% C3 is not needed for energy here

% clamp directions AS TANGENTS (per your correction), then rotate +pi/6
Rtilt = rot2(pi/6);
vecA  = unit( (Rtilt * (d1_ - d1).').' );   % left:  d1_-d1, +pi/6
vecB  = unit( (Rtilt * (A   - B ).').' );   % right: A-B,     +pi/6

% call Hencky bar-chain with tangents (VectorsAreNormals=false)
args = {'VectorsAreNormals', true, 'EndCluster', true, 'ClusterRatio', 0.9, ...
        'TolC', 1e-12, 'TolStep', 1e-12, 'DoSOC', true, 'DoLenProjection', true};
if ~isempty(z0_in)
    args = [args, {'z0', z0_in}];
end
[E_lig, ~, ~, ~, info] = hbm_energy(l2, d1, B, vecA, vecB, N, E, b, t, args{:});

E_one     = E_lig;
pack.A    = A;
pack.B    = B;
pack.info = info;

% build next warm start state if available
if isfield(info,'kappa') && isfield(info,'ell') && numel(info.kappa) >= 2
    z0_next = [info.kappa(1:end-1); info.ell(:)];
else
    z0_next = [];
end
end

function R = reaction_at_B1(info)
% Map augmented Lagrangian multipliers back to physical reaction at B.
lam = info.lambda(:);
if numel(lam) >= 2 && isfield(info,'constr_unscaled') && isfield(info,'constr_scaled')
    g_un = info.constr_unscaled(:);
    g_sc = info.constr_scaled(:);
    if norm(g_sc) > 0
        ascale_est = norm(g_un)/norm(g_sc);
    else
        ascale_est = 1.0;
    end
    R = (lam(1:2).') / ascale_est;   % [Rx, Ry]
else
    R = [NaN, NaN];
end
end
