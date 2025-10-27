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
end
