function [triangle_new, E_total] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta, t)
% Continuous deployment via 1-D θ-minimization per delta, no triangle_unit calls inside loop.
% - Left end fixed: A = B0 (from delta=0 reference).
% - vecA fixed: [cos(pi/6+beta), -sin(pi/6+beta)] (unit).
% - B(θ,δ) = G(δ) + r [cosθ, sinθ],  C(θ,δ) = G(δ) + r [cos(θ-2π/3), sin(θ-2π/3)].
% - vecB = (C - B) / l3.  L0 = l2.

%% helpers
wrap  = @(th) atan2(sin(th),cos(th));

%% Geometric derived lengths (independent of delta)
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(pi/3 - beta);
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) - cos(pi/3 - beta) .* l4;
l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 ...
    - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) );
r  = (sqrt(3)/3) * l3;   % centroid->vertex radius

% HBM constants
N = 10;  E = 4.3e11;  b = 1.0;

%% Normalize delta as row
delta = delta(:).';  nD = numel(delta);

%% Reference geometry at delta = 0 (one-time) to define indices and B0
prev_alpha_1 = pi/3;  prev_alpha_2 = 2*pi/3;
[tri0,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta(1), beta, edgeLen, l1, l4, t);

% Fixed flank nodes (rigid left side)
B0  = tri0(20,:);     % point on flank
vecB0 = tri0(20,:) - tri0(21,:); % tangent direction
%vecB0 = [cos(-pi/6+beta),sin(-pi/6+beta)];
B  = tri0(44,:);     % inner triangle vertex at delta=0 used as A = B0

% Initialize outputs
E_total   = nan(1, nD);
theta_all = nan(1, nD);

% Initial theta (given)
theta_prev = -pi + beta;
bracket_w  = pi/40;      % search range

%triangle_new = tri0; % keep the triangle_new(will update during iterations)

%% Loop over deltas (no triangle_unit calls here)
for k = 1:nD
    del  = delta(k);
    edge = edgeLen + del;

    % Centroid G(delta) (from your unit coordinates)
    G = [-sqrt(3)/6*edge, -1/2*edge];

    % 1-D minimization in theta (warm start from previous theta)
    [E_one, pack] = energy(theta_prev, N);

    % Save θ, energy
    theta_all(k) = pack.theta;
    E_total(k)   = 3 * E_one;
    % warm start update
    theta_prev = th_opt;
end

% Build optimized geometry for last delta
[~, pack, ~, ~] = energy_theta(th_opt);
[triangle_new,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta(end), beta, edgeLen, l1, l4, t);
B_new = pack.B;  C_new = pack.C;  Gcur = pack.G;
A_new = Gcur + r * [cos(th_opt - 2*pi/3), sin(th_opt - 2*pi/3)];

% Update inner triangle vertices
triangle_new(43,:) = A_new;
triangle_new(44,:) = B_new;
triangle_new(45,:) = C_new;

% Filaments (visual offsets using thickness t). Left flank nodes are from tri0 (rigid).
triangle_new(32,:) = triangle_new(20,:);
triangle_new(33,:) = triangle_new(44,:);
triangle_new(34,:) = triangle_new(33,:) + t/(norm(A_new-B_new))*(A_new-B_new);

triangle_new(37,:) = triangle_new(45,:);
triangle_new(36,:) = triangle_new(24,:);
triangle_new(38,:) = triangle_new(37,:) + t/(norm(B_new-C_new))*(B_new-C_new);

triangle_new(40,:) = triangle_new(28,:);
triangle_new(41,:) = triangle_new(43,:);
triangle_new(42,:) = triangle_new(41,:) + t/(norm(C_new-A_new))*(C_new-A_new);

%% Define function to calculate energy given th, and N
    function [E_one, pack] = energy(th, N)
        % ENERGY solve for a single ligament with clustered segmentation.
        % Inputs:
        %   th : initial guess for theta
        %   N  : number of segments with clustering near ends

        t_eff = t * sqrt(3)/2;                     % effective thickness
        EI = E * (b * t_eff^3) / 12;               % Euler-Bernoulli bending rigidity
        EA = E * (b * t_eff);                      % axial rigidity

        % Divide ligament into N segments
        a0vec = href_end_cluster(l2, N, r, 2); % end clustering
        %a0vec = (L0/N)*ones(N,1);      % evenly distrbuted
        hhinge = 0.5*(a0vec(1:end-1) + a0vec(2:end));          % K x 1

        % stiffness
        kb_vec = EI ./ hhinge;          % N-1 x 1 (bending)
        ka_vec = EA ./ a0vec;           % N x 1 (axial)

        % initial guess for optimised variables x = [phi; e; theta]
        DeltaTot = wrap(th - pi - beta);       % total rotation angle
        phi0 = (DeltaTot / N-1) * ones(N-1,1);   % distribute over free hinges
        e0   = zeros(N,1);                     % change in length of each segment
        x0   = [phi0; e0; th];

        % Define lower bounds and upper bounds of variables
        lb = [ -pi*ones(N-1,1);  a0vec.*(1 - 0.4); -pi/10];
        ub = [  pi*ones(N-1,1);  a0vec.*(1 + 0.4); pi/10];


        opts = optimoptions('fmincon', ...
            'Algorithm','interior-point', ...
            'SpecifyObjectiveGradient',true, ...
            'SpecifyConstraintGradient',true, ...
            'Display','off', ...
            'MaxIterations',200, ...
            'OptimalityTolerance',1e-10, ...
            'ConstraintTolerance',1e-10, ...
            'StepTolerance',1e-12);

        [x_opt, fval] = fmincon(@obj_fun, x0, [],[],[],[], lb, ub, @cons_fun, opts);
        % outputs
        phi   = x_opt(1:K);
        e     = x_opt(K+1:K+N);
        theta = x_opt(end);

        % geometry for pack
        R = @(a) [cos(a), -sin(a); sin(a), cos(a)];
        ex = [1;0];
        r_ = l3/sqrt(3);
        B = G + (r_ * (R(theta)            * ex)).';
        C = G + (r_ * (R(theta + 2*pi/3)   * ex)).';
        vecB = (C - B)/l3;

        E_one = fval;
        pack = struct('phi',phi,'e',e,'theta',theta, ...
            'B',B,'C',C,'G',G,'vecB',vecB, ...
            'kb',kb_vec,'ks',ks_vec,'a0vec',a0vec,'hhinge',hhinge);

        % Define objective function
        function [f, g, pack_obj] = obj_fun(x)
            % x = [phi; e; theta]
            phi = x(1:N-1);
            e   = x(N:2*N-1);
            % E = 1/2 (phi'Kb phi + e'Ks e)
            f = 0.5*( phi.'*(kb_vec.*phi) + e.'*(ka_vec.*e) );
            if nargout > 1
                g = [kb_vec.*phi; ka_vec.*e; 0];     % dE/dtheta = 0
            end
            if nargout > 2
                theta = x(end);
                B = G + r * [cos(th), sin(th)];   % r = l3/sqrt(3)
                C = G + r * [cos(th + 2*pi/3), sin(th + 2*pi/3)];
                vecB = (C - B) / l3;
                pack_obj.B = B; pack_obj.C = C; pack_obj.vecB = vecB; pack_obj.theta = theta;
            end
        end

        % Define constraints
        function [c, ceq, gc, J] = cons_fun(x)
            phi = x(1:N-1);
            e   = x(N:2*N-1);
            theta = x(end);

            % Position clousure condition
            psi = zeros(N,1);
            psi(1) = pi/6-beta;
            if N-1 > 0
                s = cumsum(phi);                  % s(i)=sum_{j=1}^i phi_j
                psi(2:end) = pi/6-beta + s(1:end);
            end
            u   = [cos(psi),  sin(psi)];          % N×2
            up  = [-sin(psi), cos(psi)];          % N×2
            res_pos = sum(((a0vec + e).*u), 1).' - ( G + (l3/sqrt(3))*[cos(theta); sin(theta)] - B0 );

            % Angle condition
            res_ang = sum(phi) - theta + pi + beta;

            ceq = [res_pos/l2; res_ang];   % 3×1
            c   = []; gc = [];

            if nargout > 3
                Jpos_phi = zeros(2, K);
                for i = 1:K
                    idx = (i+1):N;
                    Jpos_phi(:, i) = sum( (a0vec(idx)+e(idx)) .* up(idx,:), 1 ).';
                end
                uthp = [-sin(theta); cos(theta)];
                r_   = l3/sqrt(3);                % r = l3/√3
                J_theta_phi = - r_ * uthp * ones(1, K);
                Jphi = Jpos_phi + J_theta_phi;

                % d rpos / d e_m = u(ψ_m)
                Jpos_e  = u.';                    % 2×N

                % d rpos / d θ = - r u'(θ)
                Jpos_th = - r_ * uthp;            % 2×1

                % angle row
                Jang_phi = ones(1, K);
                Jang_e   = zeros(1, N);
                Jang_th  = -1;

                J = [ Jphi, Jpos_e, Jpos_th ;
                    Jang_phi, Jang_e, Jang_th ];
            end
        end

        function href = href_end_cluster(L0, N, r, p)
            % TRUE end-clustered segment lengths that sum to L0.
            % r in (0,1): clustering strength
            eps = (1 - r);
            % segment centers in (0,1)
            s = ((1:N)' - 0.5)/N;
            d = min(s, 1 - s);               % distance to nearest end
            % weight small at ends, large in middle
            w = (eps + d).^p;
            % normalize to total length L0
            href = (w / sum(w)) * L0;
        end
    end
figure('Color','w');
box on;
plot(delta, E_total);
end