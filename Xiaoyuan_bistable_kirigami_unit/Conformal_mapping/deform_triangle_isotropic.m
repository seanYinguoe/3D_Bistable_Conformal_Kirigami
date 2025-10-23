function [triangle_new,E_total] = deform_triangle_isotropic(delta,edgeLen,l1,l4,beta,t)
% Deform an isotropic triangle and solve theta with fast Newton/Armijo.
%
% Inputs:
%   delta     : signed edge-length change for the outer equilateral triangle
%   edgeLen   : original outer edge length
%   l1, l4    : flank geometry parameters
%   beta      : flank tilt angle (rad)
%   t         : filament thickness (用于可视化厚度/偏置)
%
% Outputs:
%   triangle_new : deformed unit key coordinates (same layout you used)
%   E_total      : total energy (3 ligaments)
%
% Notes:
% - Left clamp tangent: vecA = (d1 - d1_), then +pi/6 (if you want a tilt, add it yourself)
% - Right clamp tangent: vecB = (C  - B )
% - Passed as TANGENTS => 'VectorsAreNormals', false
% - Theta is solved by Newton/secant + Armijo using HBM reaction as gradient;
%   if info is missing, we fallback to finite-difference (FD) gradient.
%
% 英文为主，括号里附简短中文注释

%% ---------- helpers ----------
rotation = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)]; % 2D rotation
nrm1     = @(v) v / max(norm(v), eps);                              % normalize
wrapA    = @(th) atan2(sin(th), cos(th));                           % wrap to (-pi,pi]

%% ---------- geometry precompute ----------
% your original formulas
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(pi/3 - beta); % ligament reference length L0
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) - cos(pi/3 - beta) .* l4;
l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) );
R  = sqrt(3)/3 * l3;                       % inner rigid triangle radius (理论半径)
N  = 15;  E = 4.3e11;  b = 1.0;

%% ---------- initial triangle (your way) ----------
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;
[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1, l4, t);
triangle_new   = triangle;

% flank nodes (follow your current indices)
d1  = triangle(20,:);
d1_ = triangle(21,:);

% outer centroid in your orientation
edge     = edgeLen + delta;
centroid = [-sqrt(3)/6*edge, -1/2*edge];

% initial B and theta0 from the triangle you built
B_orig = triangle(44,:);
theta0 = atan2(B_orig(2)-centroid(2), B_orig(1)-centroid(1));

%% ---------- energy wrapper: one-ligament energy @ theta ----------
% Returns: E_oneLig (scalar), pack (B,C,info), XYdef, springs
    function [E_oneLig, pack, XYdef, springs, info] = energy(th)
        % place B on circle about outer centroid (与你原来一致)
        B = centroid + R*[cos(th), sin(th)];
        % choose the "other" inner vertex C by -120° rotation
        C = centroid + (B - centroid) * rotation(-2*pi/3);

        % tangents at clamps (注意：传切向，VectorsAreNormals=false)
        vecA = (d1 - d1_);        % left clamp tangent
        vecB = (C  - B );         % right clamp tangent
        vecA = nrm1(vecA);
        vecB = nrm1(vecB);

        % call HBM (EI = E*b*t^3/12, EA = E*b*t; here we kept your t*sqrt(3)/2)
        % IMPORTANT: pass tangents -> 'VectorsAreNormals', false
        try
            [E_oneLig, XYdef, springs, ~, info] = hbm_energy( ...
                l2, d1, B, vecA, vecB, N, E, b, t*sqrt(3)/2, ...
                'VectorsAreNormals', false);
            if ~isfinite(E_oneLig), E_oneLig = 1e50; end
        catch
            % if HBM fails (near-singular, etc.), return big penalty so line-search backs off
            E_oneLig = 1e50; XYdef = []; springs = []; info = struct();
        end

        if nargout > 1
            pack.B = B; pack.C = C; pack.info = info;
        end
    end

%% ---------- fast 1D theta solve: Newton/secant + Armijo ----------
theta_it = theta0;                 % warm start
maxIt    = 20;
tol_g    = 1e-10;                  % |dE/dθ| tolerance
damp     = 0.8;                    % damping for Newton step
last_good = struct('theta', theta0, 'E', inf);  % track best-so-far

for it = 1:maxIt
    [E_now, pack, XYdef, springs, info] = energy(theta_it);
    if E_now < last_good.E
        last_good.theta = theta_it; last_good.E = E_now;
    end

    % gradient g(θ) = -R * (R_B · t_perp)  (用反力做导数；缺失则退化为差分)
    Rvec   = reaction_at_B1(pack.info);     % [Rx,Ry] or [NaN,NaN]
    t_perp = [-sin(theta_it), cos(theta_it)];
    gtheta = - R * (Rvec * t_perp.');       % scalar

    % fallback: finite-difference if gradient unusable（若梯度异常，改用差分）
    if ~isfinite(gtheta) || abs(gtheta) < 1e-20
        h = 1e-5;
        E_p = energy(theta_it + h); 
        gtheta = (E_p - E_now) / h;
    end

    if abs(gtheta) < tol_g
        break;   % converged
    end

    % numeric curvature H(θ) ≈ d²E/dθ² （数值二阶导）
    h = 5e-4;
    E_f = energy(theta_it + h);
    g_f = (E_f - E_now) / h;
    Htheta = (g_f - gtheta)/h;

    % step proposal: damped Newton if curvature OK, else gradient step
    if isfinite(Htheta) && Htheta > 0
        dtheta = - damp * gtheta / Htheta;
    else
        dtheta = - damp * gtheta;
    end

    % Armijo backtracking（回溯线搜索，单调下降）
    alpha = 1.0;
    while true
        theta_try = wrapA(theta_it + alpha * dtheta);
        E_try = energy(theta_try);
        if E_try <= E_now || alpha < 1e-8
            theta_it = theta_try; 
            break;
        end
        alpha = 0.5 * alpha;
    end
end

% accept best found (robust against non-ideal gradients)
if last_good.E < inf
    theta_opt = last_good.theta;
else
    theta_opt = theta_it;
end

% final energy at optimal theta
[E_oneLig, pack, XYdef, ~] = energy(theta_opt);
E_total = 3 * E_oneLig;

%% ---------- rebuild final geometry (your layout) ----------
B_new = pack.B;
C_new = pack.C;
A_new = centroid + (B_new - centroid) * rotation(2*pi/3); 

triangle_new(43,:) = A_new;
triangle_new(44,:) = B_new;
triangle_new(45,:) = C_new;

% filaments (keep your pattern)
triangle_new(32,:) = triangle_new(20,:);
triangle_new(33,:) = triangle_new(44,:);
triangle_new(34,:) = triangle_new(33,:) + t/(norm(A_new-B_new)+eps)*(A_new-B_new);

triangle_new(37,:) = triangle_new(45,:);
triangle_new(36,:) = triangle_new(24,:);
triangle_new(38,:) = triangle_new(37,:) + t/(norm(B_new-C_new)+eps)*(B_new-C_new);

triangle_new(40,:) = triangle_new(28,:);
triangle_new(41,:) = triangle_new(43,:);
triangle_new(42,:) = triangle_new(41,:) + t/(norm(C_new-A_new)+eps)*(C_new-A_new);

% quick plot (kept)
figure(); hold on; axis equal; box on; grid on;
plot_triangle(triangle_new);
if ~isempty(XYdef)
    plot(XYdef(:,1), XYdef(:,2), '-o', 'Color', [1 0 0], 'DisplayName', 'Deformed');
end
title('Isotropic unit — fast \theta solve (Newton/Armijo)');
legend show;

end

%% ===== helper: map AL multipliers -> reaction at B (scale-aware) =====
function R = reaction_at_B1(info)
% Return [Rx,Ry] at the right clamp using AL/Lagrange multipliers in info.
% If 'constr_unscaled'/'constr_scaled' exist, auto-correct the scaling.
R = [NaN, NaN];
if ~isstruct(info) || ~isfield(info,'lambda'), return; end
lam = info.lambda(:);
if numel(lam) < 2, return; end
ascale = 1.0;
if isfield(info,'constr_unscaled') && isfield(info,'constr_scaled')
    gu = info.constr_unscaled(:); gs = info.constr_scaled(:);
    if norm(gs) > 0, ascale = norm(gu)/norm(gs); end
end
R = (lam(1:2).') / ascale;   % [Rx,Ry]
end
