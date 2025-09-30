function [Etotal, XY, springs, error] = hbm_energy(L0, A1,B1, vecA1,vecB1,  N, E, b, t, varargin)
% Hencky bar-chain with bending + axial springs
%   Etotal  : total elastic energy (bending + axial)
%   XY      : (N+1)x2 node coordinates of the deformed centerline
%   springs : struct with coordinates to plot springs
%             - springs.axial       : N x 4  -> [x_i y_i x_{i+1} y_{i+1}]
%             - springs.rotational  : (N-1) x 2 -> joint points [x y]
%
% Inputs:
%   A0,B0        : undeformed end points (for reference length a0 = |B0-A0|/N)
%   vecA0,vecB0  : undeformed end vectors (only used if you set 'VectorsAreNormals', true)
%   A1,B1        : deformed end points (constraints)
%   vecA1,vecB1  : deformed end vectors (clamped directions at ends)
%   N            : number of links (>=2)
%   E,b,t        : material/section (EI = E*b*t^3/12, EA = E*b*t)

p = inputParser;
addParameter(p,'VectorsAreNormals',true,@islogical);
parse(p,varargin{:});
vectorsAreNormals = p.Results.VectorsAreNormals;

% helpers
wrap  = @(th) atan2(sin(th),cos(th)); % wrap angles to [-pi,pi]
v2ang = @(v) atan2(v(2),v(1)); % compute angle of 2D vector
nrm1  = @(v) v./max(norm(v),eps); % normalize vector

% convert normals
if vectorsAreNormals
    vecA1 = normal_to_tangent(vecA1, A1, B1);
    vecB1 = normal_to_tangent(vecB1, A1, B1);
else
    vecA1 = nrm1(vecA1); vecB1 = nrm1(vecB1);
end
thetaA1 = v2ang(vecA1);  thetaB1 = v2ang(vecB1);

% reference length per segment from undeformed geometry
if N < 2, error('N must be >= 2.'); end
a0 = L0 / N;

EI = E * b * t^3 / 12;
EA = E * b * t;

% solve the deformed configuration
[phi, XY, kappa, ell, Eb, Ex] = ...
    solve_chain_axial(A1,B1,thetaA1,thetaB1,N,a0,EI,EA,wrap);

% outputs
Etotal = Eb + Ex;

% springs geometry for plotting
springs.axial = [XY(1:end-1,1) XY(1:end-1,2) XY(2:end,1) XY(2:end,2)]; % N x 4
springs.rotational = XY(2:end-1,:);                                    % (N-1) x 2

% compute and print equality-constraint residuals
dxdy   = XY(end,:).' - B1(:);              % [dx; dy]
error = norm(dxdy);
dtheta = wrap(phi(end) - thetaB1);         % angle residual

fprintf('ceq residuals: dx=%.3e, dy=%.3e, dtheta=%.3e rad (||pos||=%.3e)\n', ...
        dxdy(1), dxdy(2), dtheta, norm(dxdy));

end

%% minimize Ebend + Eax with end position & end angle constraints
% solve: z = [kappa(1..N-1); ell(1..N)] curvature;elongation
% axial energy per segment: (1/2)*(EA/a0)*(ell_i - a0)^2  (k=EA/a0)
% bending energy per hinge : (1/2)*(EI/a0)*(kappa_i)^2    (k=EI/a0)
function [phi, XY, kappa_opt, ell_opt, Eb, Ex, exitflag, output] = ...
    solve_chain_axial(A, B, thetaA, thetaB, N, a0, EI, EA, wrap)

K  = N-1;
% initial guess
k0 = (wrap(thetaB - thetaA)/K) * ones(K,1) + 1e-6*randn(K,1); % initial curvature
e0 = a0*ones(N,1); % initial elongation
z0 = [k0; e0]; % initial optimised variables
stretch = 0.4;

% bounds on curvatures and lengths
lb = [-ones(K,1)*pi;  e0*(1-stretch)];
ub = [ones(K,1)*pi;  e0*(1+stretch)];

% objective and constraints
obj  = @(z) obj_total(z, K, EI, EA, a0);
nonl = @(z) cons_end(z, A, B, thetaA, thetaB, N);

opts = optimoptions('fmincon', 'Algorithm','sqp', ...
    'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',true, ...
    'Display','off', ...
    'MaxIterations',2000, ...
    'ConstraintTolerance',1e-8, ...
    'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-10,...
    'HessianApproximation','lbfgs');       % robust with analytic J;

[z, ~, exitflag, output] = fmincon(obj, z0, [],[],[],[], lb, ub, nonl, opts);

kappa_opt = z(1:K);
ell_opt   = z(K+1:end);

% geometry
phi = thetaA + [0; cumsum(kappa_opt(:))];  % N x 1 (segment angles)
XY  = forward_chain(A, phi, ell_opt);

% energies (report)
Eb = 0.5*(EI/a0) * sum(kappa_opt.^2);
Ex = 0.5*(EA/a0) * sum((ell_opt - a0).^2);
end

%% Define the objective function(energy function)
function [f, g] = obj_total(z, K, EI, EA, a0)
kappa = z(1:K);
ell   = z(K+1:end);
f     = 0.5*(EI/a0)*sum(kappa.^2) + 0.5*(EA/a0)*sum((ell - a0).^2);
if nargout>1
    % gradients for end-spring terms
    g = [ (EI/a0)*kappa ; (EA/a0)*(ell - a0) ];
end
end

%% Define the constraints function(boundary conditions)
function [c, ceq, gc, gceq] = cons_end(z, A, B, thetaA, thetaB, N)
K     = N-1;
kappa = z(1:K);
ell   = z(K+1:end);

phi = thetaA + [0; cumsum(kappa(:))];   % N x 1
t   = [cos(phi), sin(phi)];             % N x 2

XYend = A + sum(ell(:).*t, 1);
pos_err = XYend.' - B(:);
ang_err = atan2(sin(phi(end)-thetaB), cos(phi(end)-thetaB));
a0 = max(norm(B - A)/N, 1e-12);
% scaled constraints 
%ceq = [pos_err*1e2; ang_err*1e2];
ceq = [pos_err / a0; ang_err * 1e2];
c   = [];

if nargout>2
    dPos_dell = t.';                    % 2 x N
    dPos_dk   = zeros(2,K);
    for k = 1:K
        idx = (k+1):N;
        if ~isempty(idx)
            dPos_dk(:,k) = [ -sum(ell(idx).*sin(phi(idx)));  sum(ell(idx).*cos(phi(idx))) ];
        end
    end
    dAng_dk   = ones(1,K);
    dAng_dell = zeros(1,N);

    G    = [dPos_dk, dPos_dell; dAng_dk, dAng_dell]; % 3 x (K+N)
    gceq = G.';                                      % (K+N) x 3
    gc   = [];
end
end

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
