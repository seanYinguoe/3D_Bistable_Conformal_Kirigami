function [Etotal, XY, springs, err, info] = hbm_energy(L0, A1,B1, vecA1,vecB1, N, E, b, t, varargin)
% Newton-KKT solver for Hencky bar-chain (bending + axial)
% Replaces fmincon with stiffness + constraints linear solve

p = inputParser;
addParameter(p,'VectorsAreNormals',true,@islogical);
addParameter(p,'MaxIter',1000,@(x)isnumeric(x)&&x>0);
addParameter(p,'TolC',1e-10,@(x)x>0);
addParameter(p,'TolStep',1e-12,@(x)x>0);
addParameter(p,'Damping',1.0,@(x)x>0);    % line-search damping
parse(p,varargin{:});
opt = p.Results;

wrap  = @(th) atan2(sin(th),cos(th));
v2ang = @(v) atan2(v(2),v(1));
nrm1  = @(v) v./max(norm(v),eps);


% convert normals -> tangents
if opt.VectorsAreNormals
    vecA1 = normal_to_tangent(vecA1, A1, B1);
    vecB1 = normal_to_tangent(vecB1, A1, B1);
else
    vecA1 = nrm1(vecA1); vecB1 = nrm1(vecB1);
end
thetaA = v2ang(vecA1);  thetaB = v2ang(vecB1);

% constants
a0  = L0/N;
EI  = E*b*t^3/12;
EA  = E*b*t;

K = N-1;                % number of hinges
kr = (EI/a0);           % rotational spring stiffness per hinge
ka = (EA/a0);           % axial spring stiffness per segment

% Hessian of energy (constant diagonal)
H  = diag([kr*ones(K,1); ka*ones(N,1)]);

% initial guess (same as your fmincon version)
kappa = (wrap(thetaB-thetaA)/K)*ones(K,1);
ell   = a0*ones(N,1);
z     = [kappa; ell];

% Newton iterations
for it = 1:opt.MaxIter
    % unpack
    kappa = z(1:K);
    ell   = z(K+1:end);

    % angles and tangents
    phi = thetaA + [0; cumsum(kappa(:))];    % N x 1
    txy = [cos(phi), sin(phi)];              % N x 2

    % forward kinematics
    XY = forward_chain(A1, phi, ell);

    % constraints g(z) = 0
    pos_err = (A1(:) + sum(ell(:).*txy,1).' - B1(:));     % 2x1
    ang_err = wrap(phi(end) - thetaB);                     % 1x1
    g = [pos_err/a0; 1e2*ang_err];                         % scaling 与你一致

    % gradients of E
    gradE = [kr*kappa; ka*(ell - a0)];

    % Jacobian A = dg/dz  (与您 cons_end 里推导一致)
    % dPos/dell = t^T
    dPos_dell = txy.';             % 2 x N
    % dPos/dkappa
    dPos_dk = zeros(2,K);
    for k = 1:K
        idx = (k+1):N;
        if ~isempty(idx)
            dPos_dk(:,k) = [ -sum(ell(idx).*sin(phi(idx)));  sum(ell(idx).*cos(phi(idx))) ];
        end
    end
    dAng_dk   = ones(1,K);
    dAng_dell = zeros(1,N);

    A = [ dPos_dk/a0,  dPos_dell/a0;    % 2 x (K+N)
          1e2*dAng_dk, 1e2*dAng_dell ]; % 1 x (K+N)

    % KKT system
    KKT = [H, A'; A, zeros(3,3)];
    rhs = [-gradE; -g];

    % solve (prefer sparse backslash or LDL^T)
    sol = KKT \ rhs;
    dz  = sol(1:K+N);
    %lambda = sol(K+N+1:end); %#ok<NASGU>

    % line search (simple damping)
    z_new = z + opt.Damping*dz;

    % convergence checks
    step_norm = norm(dz);
    c_norm    = norm(g);

    z = z_new;

    if step_norm < opt.TolStep && c_norm < opt.TolC
        break;
    end
end

% final geometry & energy
kappa = z(1:K); ell = z(K+1:end);
phi = thetaA + [0; cumsum(kappa(:))];
XY  = forward_chain(A1, phi, ell);
Eb = 0.5*kr * sum(kappa.^2);
Ex = 0.5*ka * sum((ell - a0).^2);
Etotal = Eb + Ex;

springs.axial      = [XY(1:end-1,1) XY(1:end-1,2) XY(2:end,1) XY(2:end,2)];
springs.rotational = XY(2:end-1,:);
err = norm([ (A1(:) + sum(ell(:).*[cos(phi),sin(phi)],1).' - B1(:))/a0 ; 1e2*wrap(phi(end)-thetaB) ]);

info.iters   = it;
info.kr      = kr;
info.ka      = ka;
info.a0      = a0;
info.kappa   = kappa;
info.ell     = ell;
info.constr  = g;

fprintf('KKT-Newton: it=%d, ||g||=%.3e, step=%.3e, E=%.6g (Eb=%.6g, Ex=%.6g)\n',...
    it, norm(g), norm(dz), Etotal, Eb, Ex);

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