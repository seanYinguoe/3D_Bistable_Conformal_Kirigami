function [triangle_new, E_total] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta, t)
% DEFORM_TRIANGLE_ISOTROPIC
% Continuous deployment using one-layer optimization per delta.
% - Left end fixed: A = B0 (from delta=0 reference)
% - B, C move on circle around centroid G(delta)
% - Optimize x = [phi; e; theta] using fmincon

%% Helper
wrap  = @(th) atan2(sin(th), cos(th));

%% Geometry constants (independent of delta)
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(pi/3 - beta);
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) - cos(pi/3 - beta) .* l4;
l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 ...
    - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) );
r_vertex  = (sqrt(3)/3) * l3;     % centroid–vertex radius (r = l3/sqrt(3))

Emod = 4.3e11;  b = 1.0;          % material constants

%% Normalize delta
delta = delta(:).';
nD = numel(delta);

%% Reference geometry (delta = 0)
prev_alpha_1 = pi/3;  prev_alpha_2 = 2*pi/3;
[tri0,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta(1), beta, edgeLen, l1, l4, t);

B_left   = tri0(20,:);                   % left anchor
alphaL = -pi/6 + beta;                % left tangent angle
E_total   = nan(1, nD);
theta_all = nan(1, nD);
theta_prev = -pi + beta;             % initial guess
r_clust = 0.85; p_clust = 2.0;        % clustering parameters

%% Loop over deltas
for k = 1:nD
    edge = edgeLen + delta(k);
    G = [-sqrt(3)/6*edge, -1/2*edge];     % centroid

    params = struct('t',t,'E',Emod,'b',b,'beta',beta,'G',G,'B_left',B_left,...
                    'r_vertex',r_vertex,'l2',l2,'l3',l3,'alphaL',alphaL);

    [E_one, pack] = energy_lig(theta_prev, 10, params); % N=10 segments
    theta_all(k) = pack.theta;
    E_total(k)   = 3 * E_one;               % 3 ligaments total
    theta_prev   = pack.theta;
    last_pack = pack;
end

%% Build final geometry
[triangle_new,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta(end), beta, edgeLen, l1, l4, t);
th_final = last_pack.theta;
G_final  = last_pack.G;

B_new = G_final + r_vertex * [cos(th_final), sin(th_final)];
C_new = G_final + r_vertex * [cos(th_final + 2*pi/3), sin(th_final + 2*pi/3)];
A_new = G_final + r_vertex * [cos(th_final - 2*pi/3), sin(th_final - 2*pi/3)];

triangle_new(43,:) = A_new;
triangle_new(44,:) = B_new;
triangle_new(45,:) = C_new;

%% Visual filaments
triangle_new(32,:) = triangle_new(20,:);
triangle_new(33,:) = triangle_new(44,:);
triangle_new(34,:) = triangle_new(33,:) + t/(norm(A_new-B_new))*(A_new-B_new);

triangle_new(37,:) = triangle_new(45,:);
triangle_new(36,:) = triangle_new(24,:);
triangle_new(38,:) = triangle_new(37,:) + t/(norm(B_new-C_new))*(B_new-C_new);

triangle_new(40,:) = triangle_new(28,:);
triangle_new(41,:) = triangle_new(43,:);
triangle_new(42,:) = triangle_new(41,:) + t/(norm(C_new-A_new))*(C_new-A_new);

%% Plot
figure('Color','w');
box on;
plot(delta, E_total, '-o');
xlabel('\delta'); ylabel('Energy');
figure()
plot_triangle(triangle_new);

% %% Define function: energy(th, N)
%     function [E_one, pack] = energy(th, N)
%         % Solve one ligament (one-layer optimization)
%         t_eff = t * sqrt(3)/2;
%         EI = Emod * (b * t_eff^3) / 12;
%         EA = Emod * (b * t_eff);
% 
%         % Segments
%         a0vec = href_end_cluster(l2, N, r_clust, p_clust);
%         K = N - 1;
%         hhinge = 0.5 * (a0vec(1:end-1) + a0vec(2:end));
%         kb_vec = EI ./ hhinge;
%         ks_vec = EA ./ a0vec;
% 
%         % Initial guess
%         DeltaTot = wrap(th + pi - beta);
%         phi0 = (DeltaTot / K) * ones(K,1);
%         e0   = zeros(N,1);
%         x0   = [phi0; e0; th];
% 
%         % Bounds
%         lb = [ -pi*ones(K,1);   -0.6*a0vec;   -pi/40 ];
%         ub = [  pi*ones(K,1);    0.6*a0vec;    +pi/40 ];
% 
%         % Objective
%         function [f, g] = obj_fun(x)
%             phi = x(1:K);
%             e   = x(K+1:K+N);
%             f = 0.5*( phi.'*(kb_vec.*phi) + e.'*(ks_vec.*e) );
%             if nargout > 1
%                 g = [kb_vec.*phi; ks_vec.*e; 0];
%             end
%         end
% 
%         % Constraints
%         function [c, ceq, gc, J] = cons_fun(x)
%             phi   = x(1:K);
%             e     = x(K+1:K+N);
%             theta = x(end);
% 
%             % Bar directions ψ_k
%             psi = zeros(N,1);
%             psi(1) = alphaL;
%             s = cumsum(phi);
%             psi(2:end) = alphaL + s(1:end);
% 
%             % Unit vectors
%             u  = [cos(psi), sin(psi)];
%             up = [-sin(psi), cos(psi)];
% 
%             % Position closure (2×1 column)
%             res_pos = sum(((a0vec + e).*u), 1).' - ( G' + r_vertex*[cos(theta); sin(theta)] - B0' );
% 
%             % Angle constraint (1×1 scalar)
%             res_ang = sum(phi) - theta - pi + beta;
% 
%             % Stack all constraints into one 3×1 column vector
%             ceq = [res_pos; res_ang];
%             c   = [];
%             gc  = [];
% 
%             % Jacobian (3×(K+N+1))
%             if nargout > 3
%                 % d rpos / d phi_i
%                 Jpos_phi = zeros(2, K);
%                 for i = 1:K
%                     idx = (i+1):N;
%                     Jpos_phi(:, i) = sum( (a0vec(idx)+e(idx)) .* up(idx,:), 1 ).' ;
%                 end
% 
%                 uthp = [-sin(theta); cos(theta)];
%                 % J_theta_phi = -r_vertex * uthp * ones(1, K);
%                 % Jphi = (Jpos_phi + J_theta_phi) / l2;
% 
%                 % d rpos / d e_m = u(ψ_m)
%                 Jpos_e  = (u.');     % (2×N)
% 
%                 % d rpos / d θ = - r u'(θ)
%                 Jpos_th = (-r_vertex * uthp);  % (2×1)
% 
%                 % Angle row derivatives
%                 Jang_phi = ones(1, K);
%                 Jang_e   = zeros(1, N);
%                 Jang_th  = -1;
% 
%                 % Combine (3×((K)+(N)+1))
%                 J = [Jpos_phi, Jpos_e, Jpos_th;
%                     Jang_phi, Jang_e, Jang_th].';
%             end
%         end
% 
%         % Optimization
%         opts = optimoptions('fmincon', ...
%             'Algorithm','interior-point', ...
%             'SpecifyObjectiveGradient',true, ...
%             'SpecifyConstraintGradient',true, ...
%             'Display','off', ...
%             'MaxIterations',500, ...
%             'OptimalityTolerance',1e-13, ...
%             'ConstraintTolerance',1e-13, ...
%             'StepTolerance',1e-12);
% 
%         [x_opt, fval] = fmincon(@obj_fun, x0, [], [], [], [], lb, ub, @cons_fun, opts);
% 
%         % Extract results
%         phi   = x_opt(1:K);
%         e     = x_opt(K+1:K+N);
%         theta = x_opt(end);
% 
%         Bp = G + r_vertex * [cos(theta), sin(theta)];
%         Cp = G + r_vertex * [cos(theta + 2*pi/3), sin(theta + 2*pi/3)];
%         vecBp = (Cp - Bp) / l3;
% 
%         E_one = fval;
%         pack = struct('phi',phi,'e',e,'theta',theta, ...
%             'B',Bp,'C',Cp,'G',G,'vecB',vecBp, ...
%             'kb',kb_vec,'ks',ks_vec,'a0vec',a0vec,'hhinge',hhinge);
%     end
% 
% %% Clustering helper
%     function href = href_end_cluster(L0, Nloc, rloc, ploc)
%         % End-clustered segment lengths that sum to L0
%         if nargin < 4, ploc = 2.0; end
%         rloc = max(1e-6, min(0.9999, rloc));
%         s = ((1:Nloc)' - 0.5)/Nloc;
%         d = min(s, 1 - s);
%         w = ((1 - rloc) + d).^ploc;
%         href = (w / sum(w)) * L0;
%     end
end
