function [triangle_out,energy_b,energy_s] = deform_triangle(q1,q2,q3,edgeLen,l1,l4,beta,t,i_out)
% DEFORM_TRIANGLE Deforms a triangle based on input node positions and edge lengths.
%
% Inputs:
%   q1, q2, q3   - Unit node coordinates (1x3 vectors)
%   edgeLen      - Original triangle edge length
%   l1, l4, beta   - l1: length of flanks l4:thickness of flanks
%   beta:tilting angle
%   t            - Thickness of filaments
%   i_out        - Orientation flag (0 for upwards, 1 for downwards)
%
% Output:
%   triangle_out - Deformed triangle coordinates

%% Compute length and mean strain
edge1 = norm(q1-q2);
edge2 = norm(q2-q3);
edge3 = norm(q1-q3);

stretch_facs = [edge1; edge2; edge3] / edgeLen;
strain = mean(stretch_facs);

%% Calculate geometric parameters
A = pi/3-beta;
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(A); % length of filament
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) ...
     - cos(pi/3 - beta) .* l4;
l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 ...
     - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) );
%% Generate uniformaly deployed triangle
%delta = (strain-1) * edgeLen;
delta = 0;
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;

[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);
% Plot the reference(initial) triangle
% figure(1)
% colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255};
% plot_triangle(triangle_new,colour)  % Plot the results and outline triangle
% patch('Vertices', triangle_new([22,26,30],:), 'Faces', [1,2,3], ...
%     'FaceColor', 'none', 'FaceAlpha', 0.5, 'EdgeColor', 'black' ...
%     ,'LineWidth', 1.5);
% axis off

%% Define orignal boundary vertices
tri_orig = [triangle(22,:);triangle(26,:);triangle(30,:)];
tri_orig = tri_orig(:);

%% Define deformed boundary vertices
p1 = [0,0];
p3 = [0,-edge3];
p2_y = (edge2^2 - edge1^2 - edge3^2) / (2 * edge3);
p2_x = -sqrt(edge1^2 - p2_y^2);
p2 = [p2_x,p2_y];
tri_new = [p1;p2;p3];
tri_new = tri_new(:);

%% Move and Rotate Flanks to Fit Outer Triangle
triangle_new = zeros(size(triangle));

% Indecies of inner triangle
f_flank = [19 20 21 22;
    23 24 25 26;
    27 28 29 30];

% Move the flanks to fit the outer boundary(rigid conditions)
flank1_t = p1 - triangle(22,:);
flank1 = triangle(f_flank(1,:),:) + flank1_t;
flank2_t = p2 - triangle(26,:);
flank2 = triangle(f_flank(2,:),:) + flank2_t;
flank3_t = p3 - triangle(30,:);
flank3 = triangle(f_flank(3,:),:) + flank3_t;

% Rotate the flanks to fit the outer boundary
vector1 = flank1(3,:) - flank1(4,:);
vector2 = p2 - p1;
flank1_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));
vector1 = flank2(3,:) - flank2(4,:);
vector2 = p3 - p2;
flank2_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));
vector1 = flank3(3,:) - flank3(4,:);
vector2 = p1 - p3;
flank3_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));

flank1 = real((flank1-p1)*rotation(-flank1_r) + p1);
flank2 = real((flank2-p2)*rotation(-flank2_r) + p2);
flank3 = real((flank3-p3)*rotation(-flank3_r) + p3);

triangle_new(19:30,:) = [flank1;flank2;flank3];


%% Define optimise parameters
% Nodes in flanks
d1 = flank1(2,:);
d2 = flank2(2,:);
d3 = flank3(2,:);
d1_ = flank1(1,:);
d2_ = flank2(1,:);
d3_ = flank3(1,:);

% Current rotational spring connecting flanks
fk1 = [d1;d2;d3];
fk1 = fk1(:);
fk1_ = [d1_;d2_;d3_];
fk1_ = fk1_(:);

% Extract original inner triangle vertices as initial guess
A_orig = triangle(43,:);
B_orig = triangle(44,:);
C_orig = triangle(45,:);
x0 = [A_orig; B_orig; C_orig];

% The original vertices of triangle
x0 = x0(:);

% Define the optimization conditions
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'iter', ...
    'MaxIterations', 10000, ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12, ...
    'ConstraintTolerance', 1e-10);

% Run optimization
[x_opt, ~] = fmincon(@(x)objective_energy(x, fk1, fk1_, l2,t), x0, ...
    [], [], [], [], [], [], ...
    @(x)constraints(x, d1, d2, d3,d1_, d2_,d3_, l3), options);


[~,energy_b,energy_s] = objective_energy(x_opt, fk1, fk1_, l2,t);


%% Create optimised unit
% Create the new inner triangle
A_new = x_opt([1,4])';
B_new = x_opt([2,5])';
C_new = x_opt([3,6])';

triangle_new(43,:) = A_new;
triangle_new(44,:) = B_new;
triangle_new(45,:) = C_new;

% Create the new filaments
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

% Flip over the triangle, if the unit is downwards
if i_out == 0 % Upwards unit
    triangle_new = [triangle_new,zeros(size(triangle_new,1),1)];
    boundary_tri = [[p1,0];[p2,0];[p3,0]];
else  % Downwards unit
    triangle_new = [triangle_new(:,1),-triangle_new(:,2),zeros(size(triangle_new,1),1)];
    boundary_tri = [[p1,0];[p2,0];[p3,0]];
    boundary_tri = [boundary_tri(:,1),-boundary_tri(:,2),boundary_tri(:,3)];
end

%% Use barycentric coordinate system to move unit to specific grid
bc_out = zeros(size(triangle_new)); % Create barycentric coordinate
for j = 1:size(bc_out,1)
    p = triangle_new(j,:);
    bc_out(j,:) = cart2barycentric(boundary_tri,p);
end

% Covert barycentric coordinate to cartesian coordinate with new triangle
boundary_tri_new = [q1;q2;q3];
triangle_out = zeros(size(triangle_new));
for j = 1:size(bc_out,1)
    triangle_out(j,:) = bc_out(j,:) * boundary_tri_new ;
end

%% Plot the results
% Set the coulour of display
% colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255}; % The colour of void, flank, filament, Innertriangle
% % Plot the optimised triangle
% figure()
% plot_triangle(triangle_new,colour)  % Plot the results and outline triangle
% patch('Vertices', boundary_tri, 'Faces', [1,2,3], ...
%     'FaceColor', 'none', 'FaceAlpha', 0.5, 'EdgeColor', 'black' ...
%     ,'LineWidth', 1.5);
% axis off
% Plot the spatial triangle
% figure()
% plot_triangle(triangle_out,colour)  % Plot the results and outline triangle
% patch('Vertices', boundary_tri_new, 'Faces', [1,2,3], ...
%     'FaceColor', 'none', 'FaceAlpha', 0.5, 'EdgeColor', 'black' ...
%     ,'LineWidth', 1.5);
end

%% Define the constraint function
function [c, ceq] = constraints(x, d1, d2, d3,d1_, d2_, d3_, l3)
% Extract vertex coordinates
A = x([1,4])';
B = x([2,5])';
C = x([3,6])';

% Rigidity constraints (triangle edge lengths preserved)
ceq_rigidity = [
    norm(B - A) - l3;  % AB
    norm(C - B) - l3;  % BC
    norm(A - C) - l3;  % CA
    ];

% keep the same length of filaments l2(Keep the filaments elastic, the length can be changed)
% ceq_distance = [
%     norm(A - d3) - l2;  % A must be l2 from d3
%     norm(B - d1) - l2;  % B must be l2 from d1
%     norm(C - d2) - l2;  % C must be l2 from d2
%     ];

ceq = ceq_rigidity;

% No-overlapping condition
nodes1 = [d1_;d1;B];
nodes2 = [d2_;d2;C];
nodes3 = [d3_;d3;A];
res1 = ifoverlapping(nodes1); % if nonoverlapping, res<0
res2 = ifoverlapping(nodes2);
res3 = ifoverlapping(nodes3);

c = [res1; res2; res3];  % No inequality constraints
end

%% Difine the objective function to minimize the energy(stretch energy, bending energy)
function [cost,energy_b,energy_s] = objective_energy(x, fk1, fk1_,l2,t)
% Define the stretch stiffness and bend stiffness
E = 1; % large E could pollute the results
b = 1;
nu = 0.2;
I  = b*t^3/12;
A  = b*t;
G  = E/(2*(1+nu));
%k  = 5/6;                    % shear coeff for rectangle
%k_euler = E*I/(l2/2);
%K_b  = k_euler / (1 + (12*E*I)/(k*G*A*(l2/2)^2)); % Use Timoshenko hinge
K_b = 1/12*E*b*t^3/(l2/2);
K_s = E*b*t/(l2/2);

% vertices of current triangle
A_current = x([1,4])';
B_current = x([2,5])';
C_current = x([3,6])';

% vertices of current spring connecting to flanks
dc1 = fk1([1,4])';
dc2 = fk1([2,5])';
dc3 = fk1([3,6])';

dc1_ = fk1_([1,4])';
dc2_ = fk1_([2,5])';
dc3_ = fk1_([3,6])';

%% Calculate bending energy
energy_b = 0;

% rotational angle at Filament 1
vec1_new = dc1_ - dc1;
vec2_new = B_current - dc1;
vec3_new = dc1 - B_current;
vec4_new = A_current - B_current;

angle1_new = acos((vec1_new * vec2_new')/(norm(vec1_new) * norm(vec2_new)));
angle2_new = acos((vec3_new * vec4_new')/(norm(vec3_new) * norm(vec4_new)));

energy_b = energy_b + 1/2*K_b*((angle1_new-pi/3)^2 + (angle2_new-2*pi/3)^2);

% rotational angle at Filament 2
vec1_new = dc2_ - dc2;
vec2_new = C_current - dc2;
vec3_new = dc2 - C_current;
vec4_new = B_current - C_current;

angle1_new = acos((vec1_new * vec2_new')/(norm(vec1_new) * norm(vec2_new)));
angle2_new = acos((vec3_new * vec4_new')/(norm(vec3_new) * norm(vec4_new)));

energy_b = energy_b + 1/2*K_b*((angle1_new-pi/3)^2 + (angle2_new-2*pi/3)^2);

% rotational angle at Filament 3
vec1_new = dc3_ - dc3;
vec2_new = A_current - dc3;
vec3_new = dc3 - A_current;
vec4_new = C_current - A_current;

angle1_new = acos((vec1_new * vec2_new')/(norm(vec1_new) * norm(vec2_new)));
angle2_new = acos((vec3_new * vec4_new')/(norm(vec3_new) * norm(vec4_new)));

energy_b = energy_b + 1/2*K_b*((angle1_new-pi/3)^2 + (angle2_new-2*pi/3)^2);

%% Calculate stretch energy
% Calculate length change of each ligaments
energy_s = 0;
length1_new = norm(dc1 - B_current);
energy_s = energy_s + 1/2 * K_s * (length1_new - l2)^2;

length2_new = norm(dc2 - C_current);
energy_s = energy_s + 1/2 * K_s * (length2_new -l2)^2;

length3_new = norm(dc3 - A_current);
energy_s = energy_s + 1/2 * K_s * (length3_new - l2)^2;

% Calculate the total ealstic energy of stretch and bend
cost = energy_b + energy_s;

end
