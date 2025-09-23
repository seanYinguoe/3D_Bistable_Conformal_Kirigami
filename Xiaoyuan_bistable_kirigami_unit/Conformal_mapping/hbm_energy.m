function sol = hbm_energy(A0,B0, vecA0,vecB0,  A1,B1, vecA1,vecB1,  N, E, b, t, varargin)
% Geometrically exact Hencky bar-chain WITH axial compliance (longitudinal springs).
%
% Required inputs:
%   A0,B0         : [1x2] undeformed centerline end points
%   vecA0,vecB0   : [1x2] undeformed end vectors (normals by default)
%   A1,B1         : [1x2] deformed   centerline end points
%   vecA1,vecB1   : [1x2] deformed   end vectors (normals by default)
%   N             : number of links (segments), N >= 2
%   E,b,t         : Young's modulus, out-of-plane width, thickness
%
% Name-Value options:
%   'VectorsAreNormals' : true/false (default: true).
%                         If true, rotate by +90° and choose sign consistent with A->B.
%   'StretchBounds'     : scalar s in (0,1] or [lo hi] (default: []).
%                         If provided, enforce a0*(1-s) <= ell_i <= a0*(1+s), or lo/hi form.
%                         Helps stability with very small N.
%
% Output struct:
%   .XY0,.phi0,.kappa0,.ell0,.Ebend0,.Eax0,.E0
%   .XY1,.phi1,.kappa1,.ell1,.Ebend1,.Eax1,.E1
%   .thetaA0,.thetaB0,.thetaA1,.thetaB1
%   .EI,.EA,.a0,.L0,.N
%   .exit0/.out0, .exit1/.out1

    % ---------------- options ----------------
    p = inputParser;
    addParameter(p,'VectorsAreNormals',true,@islogical);
    addParameter(p,'StretchBounds',[],@(x)isnumeric(x) && (isscalar(x) || (isvector(x)&&numel(x)==2)));
    parse(p,varargin{:});
    vectorsAreNormals = p.Results.VectorsAreNormals;
    stretchBounds     = p.Results.StretchBounds;

    % -------------- helpers --------------
    wrap  = @(th) atan2(sin(th),cos(th));        % wrap to [-pi,pi]
    v2ang = @(v) atan2(v(2),v(1));
    nrm   = @(v) v./max(norm(v),eps);

    % -------------- sanitize end data --------------
    % Convert provided vectors → tangents if they are normals
    if vectorsAreNormals
        vecA0 = normal_to_tangent(vecA0, A0, B0);
        vecB0 = normal_to_tangent(vecB0, A0, B0);
        vecA1 = normal_to_tangent(vecA1, A1, B1);
        vecB1 = normal_to_tangent(vecB1, A1, B1);
    else
        vecA0 = nrm(vecA0); vecB0 = nrm(vecB0);
        vecA1 = nrm(vecA1); vecB1 = nrm(vecB1);
    end

    thetaA0 = v2ang(vecA0);  thetaB0 = v2ang(vecB0);
    thetaA1 = v2ang(vecA1);  thetaB1 = v2ang(vecB1);

    % -------------- reference geometry & stiffness --------------
    if N < 2, error('N must be >= 2.'); end
    L0 = norm(B0 - A0);                 % reference length from undeformed ends
    if L0 <= 0, error('A0 and B0 must be distinct.'); end
    a0 = L0 / N;                         % rest length per segment

    EI = E * b * t^3 / 12;               % bending rigidity
    EA = E * b * t;                       % axial rigidity

    % (informative) axial stretch will handle ‖B1−A1‖ > L0, so no error here

    % -------------- solve: undeformed & deformed --------------
    [phi0, XY0, kappa0, ell0, Eb0, Ex0, exit0, out0] = ...
        solve_chain_axial(A0,B0,thetaA0,thetaB0,N,a0,EI,EA,wrap,stretchBounds);

    [phi1, XY1, kappa1, ell1, Eb1, Ex1, exit1, out1] = ...
        solve_chain_axial(A1,B1,thetaA1,thetaB1,N,a0,EI,EA,wrap,stretchBounds);

    % -------------- pack output --------------
    sol.XY0=XY0; sol.phi0=phi0; sol.kappa0=kappa0; sol.ell0=ell0;
    sol.Ebend0=Eb0; sol.Eax0=Ex0; sol.E0=Eb0+Ex0;

    sol.XY1=XY1; sol.phi1=phi1; sol.kappa1=kappa1; sol.ell1=ell1;
    sol.Ebend1=Eb1; sol.Eax1=Ex1; sol.E1=Eb1+Ex1;

    sol.A0=A0; sol.B0=B0; sol.A1=A1; sol.B1=B1;
    sol.thetaA0=thetaA0; sol.thetaB0=thetaB0; sol.thetaA1=thetaA1; sol.thetaB1=thetaB1;

    sol.EI=EI; sol.EA=EA; sol.a0=a0; sol.L0=L0; sol.N=N;
    sol.exit0=exit0; sol.out0=out0; sol.exit1=exit1; sol.out1=out1;
end

% ========================================================================
% Core solver with axial compliance: minimize Ebend + Eax subject to
% end position & end angle constraints.
% Unknowns: z = [kappa(1..N-1); ell(1..N)]
% ========================================================================
function [phi, XY, kappa_opt, ell_opt, Eb, Ex, exitflag, output] = ...
    solve_chain_axial(A, B, thetaA, thetaB, N, a0, EI, EA, wrap, stretchBounds)

    K  = N-1;
    % Initial guess
    k0 = (wrap(thetaB - thetaA)/K) * ones(K,1) + 1e-6*randn(K,1);
    e0 = a0*ones(N,1);
    z0 = [k0; e0];

    % Bounds
    lb = [-inf(K,1);  1e-9*ones(N,1)];   % ell_i > 0
    ub = [];
    if ~isempty(stretchBounds)
        if isscalar(stretchBounds)
            s = stretchBounds;
            if ~(s>0 && s<=1), error('StretchBounds scalar must be in (0,1].'); end
            lo = a0*(1 - s);
            hi = a0*(1 + s);
        else
            lo = a0*(1 + stretchBounds(1));
            hi = a0*(1 + stretchBounds(2));
        end
        lb = [-inf(K,1);  lo*ones(N,1)];
        ub = [ inf(K,1);  hi*ones(N,1)];
    end

    % Objective (with gradient)
    obj   = @(z) obj_bend_axial(z, K, N, EI, EA, a0);

    % Constraints (with analytic Jacobian) — NOTE: returns gceq as (nVars x nCeq)
    nonl  = @(z) cons_end_axial(z, A, B, thetaA, thetaB, N, wrap);

    opts = optimoptions('fmincon', ...
        'Algorithm','sqp', ...
        'SpecifyObjectiveGradient', true, ...
        'SpecifyConstraintGradient', true, ...
        'Display','off', 'MaxIterations', 1000, 'MaxFunctionEvaluations', 2e6);

    [z, ~, exitflag, output] = fmincon(obj, z0, [],[],[],[], lb, ub, nonl, opts);

    kappa_opt = z(1:K);
    ell_opt   = z(K+1:end);

    % Reconstruct geometry
    phi = thetaA + [0; cumsum(kappa_opt(:))];           % N x 1
    XY  = forward_chain_varlen(A, phi, ell_opt);

    % Energies
    Eb = 0.5*(EI/a0) * sum(kappa_opt.^2);
    Ex = 0.5*EA       * sum((ell_opt - a0).^2);
end

function [f, g] = obj_bend_axial(z, K, N, EI, EA, a0)
    kappa = z(1:K);
    ell   = z(K+1:end);
    f_b   = 0.5*(EI/a0) * sum(kappa.^2);
    f_x   = 0.5*EA       * sum((ell - a0).^2);
    f     = f_b + f_x;
    if nargout>1
        gk = (EI/a0) * kappa;
        ge = EA * (ell - a0);
        g  = [gk; ge];
    end
end

function [c, ceq, gc, gceq] = cons_end_axial(z, A, B, thetaA, thetaB, N, wrap)
    % Variables
    K     = N-1;
    kappa = z(1:K);
    ell   = z(K+1:end);               % N x 1

    % Angles & tangents
    phi = thetaA + [0; cumsum(kappa(:))];   % N x 1
    t   = [cos(phi), sin(phi)];             % N x 2

    % End position
    XY  = A + sum(ell(:).*t, 1);            % 1 x 2
    pos_err = XY.' - B(:);                  % 2 x 1

    % End angle
    ang_err = wrap(phi(end) - thetaB);      % 1 x 1

    ceq = [pos_err; ang_err];
    c   = [];                               % no inequalities

    if nargout>2
        % d(pos)/d(ell_j) = t_j
        dPos_dell = t.';                    % 2 x N

        % d(pos)/d(kappa_k) = sum_{i=k+1..N} ell_i * [-sin(phi_i); cos(phi_i)]
        dPos_dk = zeros(2,K);
        for k = 1:K
            idx = (k+1):N;
            if ~isempty(idx)
                dPos_dk(:,k) = [ -sum(ell(idx).*sin(phi(idx)));  sum(ell(idx).*cos(phi(idx))) ];
            end
        end

        % d(ang)/d(kappa) = 1;  d(ang)/d(ell) = 0
        dAng_dk   = ones(1,K);
        dAng_dell = zeros(1,N);

        % Build gceq as (nVars x nCeq) = (K+N) x 3  **IMPORTANT SHAPE**
        % Stack w.r.t. variables (kappa first then ell), THEN transpose.
        gceq = [dPos_dk, dPos_dell; dAng_dk, dAng_dell].';   % 3 x (K+N)  -> transpose
        gceq = gceq.';                                       % (K+N) x 3   <-- what fmincon expects
        gc   = [];
    end
end

function XY = forward_chain_varlen(A, phi, ell)
    N  = numel(phi);
    XY = zeros(N+1,2);
    XY(1,:) = A;
    for i=1:N
        XY(i+1,:) = XY(i,:) + ell(i)*[cos(phi(i)), sin(phi(i))];
    end
end

function v_tan = normal_to_tangent(v_normal, A, B)
% Convert a clamp NORMAL vector to a consistent TANGENT pointing from A to B.
    R90  = [0 -1; 1 0];                  % +90° CCW
    v1   = (R90 * v_normal(:)).';        % candidate tangent
    v2   = -v1;                          % opposite direction
    pref = (B - A);
    if dot(v1, pref) >= dot(v2, pref), v_tan = v1; else, v_tan = v2; end
    n = norm(v_tan); if n>0, v_tan = v_tan / n; end
end
