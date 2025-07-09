function triangle_out = deform_triangle(q1,q2,q3,edgeLen,l1,l2,l3,t,i_out)
% DEFORM_TRIANGLE Deforms a triangle based on input node positions and edge lengths.
%
% Inputs:
%   q1, q2, q3   - Unit node coordinates (1x3 vectors)
%   edgeLen      - Original triangle edge length
%   l1, l2, l3   - l1: length of flanks l2:length of filament l3:length of inner triangle
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

%% Generate uniformaly deployed triangle
delta = (strain-1) * edgeLen;
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;

[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta,l1,l2,l3,t);
% figure(1)
% Plot the reference(initial) triangle
% colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255}; 
% plot_triangle(triangle,colour)  % Plot the results and outline triangle
% patch('Vertices', triangle([22,26,30],:), 'Faces', [1,2,3], ...
%     'FaceColor', 'none', 'FaceAlpha', 0.5, 'EdgeColor', 'black' ...
%     ,'LineWidth', 1.5);
% axis off

%tri_orig = [triangle(22,:);triangle(26,:);triangle(30,:)];
%tri_orig = tri_orig(:);

%% Compute outer triangle(boundary) vertices
p1 = [0,0];
p3 = [0,-edge3];
p2_y = (edge2^2 - edge1^2 - edge3^2) / (2 * edge3);
p2_x = -sqrt(edge1^2 - p2_y^2);
p2 = [p2_x,p2_y];
%tri_new = [p1;p2;p3];
%tri_new = tri_new(:);

%% Move and Rotate Flanks to Fit Outer Triangle
triangle_new = zeros(size(triangle));

% Indecies of inner triangle
f_flank = [19 20 21 22;
23 24 25 26;
27 28 29 30];

% Rigid conditions
% Move the flanks to fit the outer triangle(rigid conditions)
flank1_t = p1 - triangle(22,:);
flank1 = triangle(f_flank(1,:),:) + flank1_t;
flank2_t = p2 - triangle(26,:);
flank2 = triangle(f_flank(2,:),:) + flank2_t;
flank3_t = p3 - triangle(30,:);
flank3 = triangle(f_flank(3,:),:) + flank3_t;

% Calculate the rotational angle
vector1 = flank1(3,:) - flank1(4,:);
vector2 = p2 - p1;
flank1_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));
vector1 = flank2(3,:) - flank2(4,:);
vector2 = p3 - p2;
flank2_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));
vector1 = flank3(3,:) - flank3(4,:);
vector2 = p1 - p3;
flank3_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));

% Rotate the flanks to fit the outer triangle
flank1 = real((flank1-p1)*rotation(-flank1_r) + p1);
flank2 = real((flank2-p2)*rotation(-flank2_r) + p2);
flank3 = real((flank3-p3)*rotation(-flank3_r) + p3);

triangle_new(19:30,:) = [flank1;flank2;flank3];


%% Optimise inner triangle to minimize the distortion of filaments

% nodes in flank that connect to the inner triangle
d1 = flank1(2,:);
d2 = flank2(2,:);
d3 = flank3(2,:);
d1_ = flank1(1,:);
d2_ = flank2(1,:);
d3_ = flank3(1,:);
%fk1 = [d1;d2;d3];
%fk1 = fk1(:);
%fk0 = [triangle(20,:);triangle(24,:);triangle(28,:)];
%fk0 = fk0(:);

% Extract original triangle vertices
A_orig = triangle(43,:);
B_orig = triangle(44,:);
C_orig = triangle(45,:);
x0 = [A_orig; B_orig; C_orig];
x0 = x0(:);  % The original vertices of triangle

% Define the optimization conditions
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'iter', ...
    'MaxIterations', 10000, ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12, ...
    'ConstraintTolerance', 1e-10);  

% Run optimization
% l3 is the length of filaments and l4 is the length of triangle
[x_opt, ~] = fmincon(@(x)objective(x, x0), x0, ...
    [], [], [], [], [], [], ...
    @(x)constraints(x, d1, d2, d3, d1_, d2_,d3_, l3, l2), options); 
% [x_opt, ~] = fmincon(@(x)objective_energy(x, x0, fk1, fk0, tri_new, tri_orig), x0, ...
%     [], [], [], [], [], [], ...
%     @(x)constraints(x, d1, d2, d3, l3, l2), options); 


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

%% Use barycentric coordinate to change the location of triangle
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
% Plot the optimised triangle
% figure()
% plot_triangle(triangle_new,colour)  % Plot the results and outline triangle
% patch('Vertices', boundary_tri, 'Faces', [1,2,3], ...
%     'FaceColor', 'none', 'FaceAlpha', 0.5, 'EdgeColor', 'black' ...
%     ,'LineWidth', 1.5);
% axis off
% Plot the final triangle
% figure()
% plot_triangle(triangle_out,colour)  % Plot the results and outline triangle
% patch('Vertices', boundary_tri_new, 'Faces', [1,2,3], ...
%     'FaceColor', 'none', 'FaceAlpha', 0.5, 'EdgeColor', 'black' ...
%     ,'LineWidth', 1.5);
end

%% Define the constraint function
function [c, ceq] = constraints(x, d1, d2, d3,d1_, d2_,d3_, l3, l2)
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

% keep the same length of filaments l2
ceq_distance = [
    norm(A - d3) - l2;  % A must be l2 from d3
    norm(B - d1) - l2;  % B must be l2 from d1
    norm(C - d2) - l2;  % C must be l2 from d2
    ];

ceq = [ceq_rigidity; ceq_distance];
% No-overlapping condition 
nodes1 = [d1_;d1;B];
nodes2 = [d2_;d2;C];
nodes3 = [d3_;d3;A];
res1 = ifoverlapping(nodes1); % if nonoverlapping, res<0
res2 = ifoverlapping(nodes2);
res3 = ifoverlapping(nodes3);

c = [res1; res2; res3];  % No inequality constraints
end

%% Difine the objective function to minimize the energy
% function cost = objective_energy(x, x0, fk1, fk0, tri_new, tri_orig)
% % vertices of deformed triangle
% current_triangle = reshape(x, 2, 3)';
% A_current = current_triangle(1, :);
% B_current = current_triangle(2, :);
% C_current = current_triangle(3, :);
% 
% % vertices of original triangle
% original_triangle = reshape(x0, 2, 3)';
% A_orig = original_triangle(1, :);
% B_orig = original_triangle(2, :);
% C_orig = original_triangle(3, :);
% 
% % vertices of deformed flanks
% current_flanks = reshape(fk1, 2, 3)';
% dc1 = current_flanks(1, :);
% dc2 = current_flanks(2, :);
% dc3 = current_flanks(3, :);
% 
% % vertices of original flanks
% original_flanks = reshape(fk0, 2, 3)';
% do1 = original_flanks(1, :);
% do2 = original_flanks(2, :);
% do3 = original_flanks(3, :);
% 
% % Calculate rotational energy, the bending stiffness is considered
% % as 1
% energy = 0;
% current_boundary = reshape(tri_new, 2, 3)';
% p1_new = current_boundary(1,:);
% p2_new = current_boundary(2,:);
% p3_new = current_boundary(3,:);
% original_boundary = reshape(tri_orig, 2, 3)';
% p1_orig = original_boundary(1,:);
% p2_orig = original_boundary(2,:);
% p3_orig = original_boundary(3,:);
% 
% % Filament 1
% vec1_orig = p1_orig - p2_orig;
% vec2_orig = B_orig - do1;
% vec3_orig = do1 - B_orig;
% vec4_orig = A_orig - B_orig;
% 
% vec1_new = p1_new - p2_new;
% vec2_new = B_current - dc1;
% vec3_new = dc1 - B_current;
% vec4_new = A_current - B_current;
% 
% angle1_orig = acos((vec1_orig * vec2_orig')/(norm(vec1_orig) * norm(vec2_orig)));
% angle2_orig = acos((vec3_orig * vec4_orig')/(norm(vec3_orig) * norm(vec4_orig)));
% 
% angle1_new = acos((vec1_new * vec2_new')/(norm(vec1_new) * norm(vec2_new)));
% angle2_new = acos((vec3_new * vec4_new')/(norm(vec3_new) * norm(vec4_new)));
% 
% energy = energy + 1/2*1*((angle1_new-angle1_orig)^2 + (angle2_new-angle2_orig)^2);
% 
% % Filament 2
% vec1_orig = p2_orig - p3_orig;
% vec2_orig = C_orig - do2;
% vec3_orig = do2 - C_orig;
% vec4_orig = B_orig - C_orig;
% 
% vec1_new = p2_new - p3_new;
% vec2_new = C_current - dc2;
% vec3_new = dc2 - C_current;
% vec4_new = B_current - C_current;
% 
% angle1_orig = acos((vec1_orig * vec2_orig')/(norm(vec1_orig) * norm(vec2_orig)));
% angle2_orig = acos((vec3_orig * vec4_orig')/(norm(vec3_orig) * norm(vec4_orig)));
% 
% angle1_new = acos((vec1_new * vec2_new')/(norm(vec1_new) * norm(vec2_new)));
% angle2_new = acos((vec3_new * vec4_new')/(norm(vec3_new) * norm(vec4_new)));
% 
% energy = energy + 1/2*1*((angle1_new-angle1_orig)^2 + (angle2_new-angle2_orig)^2);
% 
% % Filament 3
% vec1_orig = p3_orig - p1_orig;
% vec2_orig = A_orig - do3;
% vec3_orig = do3 - A_orig;
% vec4_orig = C_orig - A_orig;
% 
% vec1_new = p3_new - p1_new;
% vec2_new = A_current - dc3;
% vec3_new = dc3 - A_current;
% vec4_new = C_current - A_current;
% 
% angle1_orig = acos((vec1_orig * vec2_orig')/(norm(vec1_orig) * norm(vec2_orig)));
% angle2_orig = acos((vec3_orig * vec4_orig')/(norm(vec3_orig) * norm(vec4_orig)));
% 
% angle1_new = acos((vec1_new * vec2_new')/(norm(vec1_new) * norm(vec2_new)));
% angle2_new = acos((vec3_new * vec4_new')/(norm(vec3_new) * norm(vec4_new)));
% 
% energy = energy + 1/2*1*((angle1_new-angle1_orig)^2 + (angle2_new-angle2_orig)^2);
% 
% cost = energy;
% end

%% Difine the objective function to minimize the difference
function cost = objective(x, x0)
cost = sum((x-x0).^2);
end


