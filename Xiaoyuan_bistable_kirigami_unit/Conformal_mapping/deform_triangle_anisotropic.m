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

%% Define helpers
rotation = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)]; % rotation matrix  
rotrow = @(v, ang) (rotation(ang) * v(:))';   % row -> column -> rotate -> row
nrm1   = @(v) v / norm(v);          % normalize

%% Outer-edge lengths (current unit)
edge1 = norm(q1-q2);
edge2 = norm(q2-q3);
edge3 = norm(q1-q3);

%% Geometry
Aang = pi/3 - beta;
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(Aang);      % filament length (ref)
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) - cos(pi/3 - beta) .* l4;
l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) );
R  = sqrt(3)/3 * l3;                                     % inscribed-circle radius

%% Uniformly deployed template (for topology & indices)
delta = 0;
prev_alpha_1 = pi/3; prev_alpha_2 = 2*pi/3;
[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);

%% Build reference outer triangle in a local frame (p1,p2,p3)
p1 = [0,0];
p3 = [0,-edge3];
p2_y = (edge2^2 - edge1^2 - edge3^2) / (2*edge3);
p2_x = -sqrt(max(edge1^2 - p2_y^2,0));
p2 = [p2_x,p2_y];

%% Move & rotate flanks
f_flank = [19 20 21 22; 23 24 25 26; 27 28 29 30];
flank1_t = p1 - triangle(22,:); flank1 = triangle(f_flank(1,:),:) + flank1_t;
flank2_t = p2 - triangle(26,:); flank2 = triangle(f_flank(2,:),:) + flank2_t;
flank3_t = p3 - triangle(30,:); flank3 = triangle(f_flank(3,:),:) + flank3_t;

% rotate to fit edges
vector1 = flank1(3,:) - flank1(4,:); vector2 = p2 - p1;
flank1_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));
vector1 = flank2(3,:) - flank2(4,:); vector2 = p3 - p2;
flank2_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));
vector1 = flank3(3,:) - flank3(4,:); vector2 = p1 - p3;
flank3_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));

flank1 = real((flank1-p1)*rotation(-flank1_r) + p1);
flank2 = real((flank2-p2)*rotation(-flank2_r) + p2);
flank3 = real((flank3-p3)*rotation(-flank3_r) + p3);

triangle_new = zeros(size(triangle));
triangle_new(19:30,:) = [flank1;flank2;flank3];

%% Flank "tangent" anchor pairs
d1  = flank1(2,:);  d1_ = flank1(3,:);
d2  = flank2(2,:);  d2_ = flank2(3,:);
d3  = flank3(2,:);  d3_ = flank3(3,:);

%% HBM-ANISO: define inner vertices on a circle, with independent angles
centroid0 = (p1 + p2 + p3)/3; % initial centroid

% Get the initial angle
B0 = triangle(44,:);
theta0 = atan2(B0(2)-centroid0(2), B0(1)-centroid0(1));
x0 = [centroid0(1), centroid0(2), theta0];  

% HBM
N = 10; E = 4.3e11; b = 1.0;

% Objective function: sum of energy of three ligaments
    function [E_sum, pack] = energy_x(x)
        xc = x(1); yc = x(2); th = x(3);
        centroid = [xc,yc];
        % place B on the circle around centroid
        B = centroid + R*[cos(th), sin(th)];
        % get C by rotating B around centroid by -120°
        C = centroid + (B - centroid) * rotation(-2*pi/3);
        % get A by rotating B around centroid by 120°
        A = centroid + (B - centroid) * rotation(2*pi/3);

        % tangent direnction1 = (d - d_)
        vecA1 = (d1 - d1_); vecB1 = (C  - B);   % ligament 1: d1 ↔ B
        vecA2 = (d2 - d2_); vecB2 = (A  - C);   % ligament 2: d2 ↔ C
        vecA3 = (d3 - d3_); vecB3 = (B  - A);   % ligament 3: d3 ↔ A

        % Get energy for each ligament
        [E1,XY1] = hbm_energy(l2, d1, B, vecA1, vecB1, N, E, b, t*sqrt(3)/2, 'VectorsAreNormals', false);
        [E2,XY2] = hbm_energy(l2, d2, C, vecA2, vecB2, N, E, b, t*sqrt(3)/2, 'VectorsAreNormals', false);
        [E3,XY3] = hbm_energy(l2, d3, A, vecA3, vecB3, N, E, b, t*sqrt(3)/2, 'VectorsAreNormals', false);

        E_sum = E1 + E2 + E3;

        if nargout>1
            pack.centroid = centroid;
            pack.A=A; pack.B=B; pack.C=C;
            pack.XY = {XY1,XY2,XY3};
        end
    end

% Define search ranges
dx_max = 0.08 * edgeLen;  
dy_max = 0.08 * edgeLen;
dtheta_max = pi/40;         

% Lower and upper bounds
L = [centroid0(1) - dx_max, centroid0(2) - dy_max, -dtheta_max];
U = [centroid0(1) + dx_max, centroid0(2) + dy_max,  dtheta_max];
opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','iter', ...
    'MaxIterations', 100, ...
    'OptimalityTolerance',1e-10, ...
    'StepTolerance',1e-12, ...
    'ConstraintTolerance',1e-10);
x_opt = fmincon(@(x) energy_x(x), x0, [], [], [], [], L, U, [], opts);

% ---- final geometry ----
[~, pk] = energy_x(x_opt);
A_new = pk.A; B_new = pk.B; C_new = pk.C;

%% 写回 inner triangle & filaments（与原逻辑一致）
triangle_new(43,:) = A_new;
triangle_new(44,:) = B_new;
triangle_new(45,:) = C_new;

% filaments strips
triangle_new(32,:) = triangle_new(20,:);
triangle_new(33,:) = triangle_new(44,:);
triangle_new(31,:) = triangle_new(32,:) + t/(norm(p1-p2))*(p1-p2);
triangle_new(34,:) = triangle_new(33,:) + t/(norm(A_new-B_new))*(A_new-B_new);

triangle_new(37,:) = triangle_new(45,:);
triangle_new(36,:) = triangle_new(24,:);
triangle_new(35,:) = triangle_new(36,:) + t/(norm(p2-p3))*(p2-p3);
triangle_new(38,:) = triangle_new(37,:) + t/(norm(B_new-C_new))*(B_new-C_new);

triangle_new(40,:) = triangle_new(28,:);
triangle_new(41,:) = triangle_new(43,:);
triangle_new(39,:) = triangle_new(40,:) + t/(norm(p3-p1))*(p3-p1);
triangle_new(42,:) = triangle_new(41,:) + t/(norm(C_new-A_new))*(C_new-A_new);

%% Define whether need to flip or not
if i_out == 0 % upwards
    triangle_new = [triangle_new, zeros(size(triangle_new,1),1)];
    boundary_tri = [[p1,0];[p2,0];[p3,0]];
else          % downwards
    triangle_new = [triangle_new(:,1), -triangle_new(:,2), zeros(size(triangle_new,1),1)];
    boundary_tri = [[p1,0];[p2,0];[p3,0]];
    boundary_tri = [boundary_tri(:,1), -boundary_tri(:,2), boundary_tri(:,3)];
end

%% barycentric map back to (q1,q2,q3) space
bc_out = zeros(size(triangle_new)); 
for j = 1:size(bc_out,1)
    p = triangle_new(j,:);
    bc_out(j,:) = cart2barycentric(boundary_tri,p);
end
boundary_tri_new = [q1;q2;q3];
triangle_out = zeros(size(triangle_new));
for j = 1:size(bc_out,1)
    triangle_out(j,:) = bc_out(j,:) * boundary_tri_new ;
end

%% Plot
% figure; 
% plot_triangle(triangle_new); 

end
