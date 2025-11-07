function [triangle_out,E_total] = deform_triangle_anisotropic(q1,q2,q3,edgeLen,l1,l4,beta,t,i_out)
% Inputs:
%   q1, q2, q3   - Unit node coordinates (1x3 vectors)
%   edgeLen      - Original triangle edge length
%   l1, l4, beta   - l1: length of flanks l4:thickness of flanks
%   beta:tilting angle
%   t            - Thickness of filaments
%   i_out        - Orientation flag (0 for upwards, 1 for downwards)
%
% Output:
%   triangle_new - Deformed triangle coordinates
%   E_total      - Deployed energy

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


%% Outer-edge lengths (current unit)
edge1 = norm(q1-q2);
edge2 = norm(q2-q3);
edge3 = norm(q1-q3);
p1 = [0,0];
p3 = [0,-edge3];
p2_y = (edge2^2 - edge1^2 - edge3^2) / (2*edge3);
p2_x = -sqrt(max(edge1^2 - p2_y^2,0));
p2 = [p2_x,p2_y];

%% Normalize delta
delta = delta(:).';
nD = numel(delta);

%% Reference geometry (delta = 0)
prev_alpha_1 = pi/3;  prev_alpha_2 = 2*pi/3;
[tri0,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, 0, beta, edgeLen, l1, l4, t);

% centroid
xG0 = -sqrt(3)/6*edgeLen;
yG0 = -1/2*edgeLen;
G = [xG0, yG0];     % centroid

E_total  = nan(1, nD);
theta_all = nan(1, nD);
theta_prev = -pi + beta;             % start with undefomred one

% first ligament
B_flank = tri0(20,:);                   % left anchor(on flank)
alphaL_B = -pi/6 + beta;                % left tangent angle

% second ligament
A_flank = tri0(28,:);
alphaL_A = -pi/6 + beta;

% third ligament
C_flank = tri0(36,:);
alphaL_C = -pi/6 + beta;

% input initial guess for optimisation
K = N-1;  % number of torsional spring
DeltaTot = wrap(theta_prev + pi - beta); % total angle difference
phi_prevB = (DeltaTot / K) * ones(K,1); % initial rotational angle of torsional springs
e_prevB   = zeros(N,1); % initial length change in ligaments
phi_prevA = (DeltaTot / K) * ones(K,1);
e_prevA   = zeros(N,1); 
phi_prevC = (DeltaTot / K) * ones(K,1);
e_prevC   = zeros(N,1); 
x_prev = [phi_prevB; e_prevB; phi_prevA; e_prevA; phi_prevC; e_prevC; theta_prev; xG0; yG0]; % initial guess

%% Loop over deltas
for k = 1:nD
    params = struct('t',t,'E',Emod,'b',b,'beta',beta,'B_flank',B_flank,...
        'A_flank',A_flank,'C_flank',C_flank,'r_vertex',r_vertex,'l2',l2,...
        'l3',l3,'alphaL_B',alphaL_B,'alphaL_A',alphaL_A,'alphaL_C',alphaL_C);

    [E_one, pack] = energy_lig(x_prev, N, params); % N=10 segments    

    % Save the reuslt
    theta_all(k) = pack.theta;
    E_total(k)   = 3 * E_one;               % 3 ligaments total
    theta_prev   = pack.theta;
    phi_prevB = pack.phiB;
    e_prevB = pack.eB;
    phi_prevA = pack.phiA;
    e_prevA = pack.eA;
    phi_prevC = pack.phiC;
    e_prevC = pack.eC;
    xG0 = pack.xG;
    yG0 = pack.yG;
    x_prev = [phi_prevB; e_prevB; phi_prevA; e_prevA; phi_prevC; e_prevC; theta_prev; xG0; yG0];
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
% Plot energy curve
figure('Color','w');
box on;
plot(delta, E_total, '-o');
xlabel('\delta'); ylabel('Energy');
figure();
plot_triangle(triangle_new);
end