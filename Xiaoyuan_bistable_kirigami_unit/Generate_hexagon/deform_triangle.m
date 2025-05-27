function deformed_triangle = deform_triangle(edge1,edge2,edge3,edgeLen,l3,l4)
% Input :
% edge1 : deformed edge1
% edge2 : deformed edge2
% edge3 : deformed edge3
% edgeLen  : original length
% Output:
% deformed_triangle
% The flanks, innertriangle are considered as rigid parts, and the
% filaments are considered as soft material, their length and angles can be
% changed to accomodate the geometric incompatibility due to non-uniform
% deployment

% Calculate stretch factor
stretch_facs = [edge1/edgeLen;
    edge2/edgeLen;
    edge3/edgeLen];
strain = mean(stretch_facs); % Get the mean stretch factors as expand ratio

% Get the uniformly deployed triangle
delta = (strain-1) * edgeLen;
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;

[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta,l1,l2,l3,t);
figure(1)
plot_triangle(triangle,colour)  % Plot the results and outline triangle
patch('Vertices', triangle([22,26,30],:), 'Faces', [1,2,3], ...
    'FaceColor', 'none', 'FaceAlpha', 0.5, 'EdgeColor', 'black' ...
    ,'LineWidth', 1.5);

tri_orig = [triangle(22,:);triangle(26,:);triangle(30,:)];
tri_orig = tri_orig(:);

% Apply fluctuation on each each to fit the triangle while
% Get the coordinates of triangle
p1 = [0,0];
p3 = [0,-edge3];
p2_y = (edge1^2 - edge2^2 - edge3^2) / (2 * edge3);
p2_x = -sqrt(edge2^2 - p2_y^2);
p2 = [p2_x,p2_y];
tri_new = [p1;p2;p3];
tri_new = tri_new(:);

%% Get the different parts of units
triangle_new = zeros(size(triangle));

% Indecies of inner triangle
f_void = [1  2  3  4  5  6;
 7  8  9 10 11 12;
13 14 15 16 17 18];
f_flank = [19 20 21 22;
23 24 25 26;
27 28 29 30];
f_filament = [31 32 33 34;
35 36 37 38;
39 40 41 42];
f_inner = [43 44 45];


% Move the flanks to fit the outer triangle
flank1_t = p1 - triangle(22,:);
flank1 = triangle(f_flank(1,:),:) + flank1_t;
flank2_t = p2 - triangle(26,:);
flank2 = triangle(f_flank(2,:),:) + flank2_t;
flank3_t = p3 - triangle(30,:);
flank3 = triangle(f_flank(3,:),:) + flank3_t;

% Calculate the rotational angle
vector1 = flank1(3,:) - flank1(4,:);
vector2 = p2 - p1;
flank1_r = acos(dot(vector1, vector2) / (norm(vector1) * norm(vector2)));
vector1 = flank2(3,:) - flank2(4,:);
vector2 = p3 - p2;
flank2_r = acos(dot(vector1, vector2) / (norm(vector1) * norm(vector2)));
vector1 = flank3(3,:) - flank3(4,:);
vector2 = p1 - p3;
flank3_r = acos(dot(vector1, vector2) / (norm(vector1) * norm(vector2)));

% Rotate the flanks to fit the outer triangle
flank1 = (flank1-p1)*rotation(-flank1_r) + p1;
flank2 = (flank2-p2)*rotation(-flank2_r) + p2;
flank3 = (flank3-p3)*rotation(-flank3_r) + p3;

triangle_new(19:30,:) = [flank1;flank2;flank3];

% Plot the results
figure(2)
plot_triangle(triangle_new,colour)  % Plot the results and outline triangle
% Plot boundary triangle
patch('Vertices', [p1;p2;p3], 'Faces', [1,2,3], ...
    'FaceColor', 'none', 'FaceAlpha', 0.5, 'EdgeColor', 'black' ...
    ,'LineWidth', 1.5);


%% Optimise inner triangle to minimize the distortion of filaments

% nodes in flank that connect to the inner triangle
d1 = flank1(2,:);
d2 = flank2(2,:);
d3 = flank3(2,:);
fk1 = [d1;d2;d3];
fk1 = fk1(:);
fk0 = [triangle(20,:);triangle(24,:);triangle(28,:)];
fk0 = fk0(:);

% Extract original triangle vertices
A_orig = triangle(43,:);
B_orig = triangle(44,:);
C_orig = triangle(45,:);
x0 = [A_orig; B_orig; C_orig];
x0 = x0(:);  % The original vertices of triangle

% Define the optimization problem
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'iter', ...
    'MaxIterations', 1000);

% Run optimization
% l3 is the length of filaments and l4 is the length of triangle
[x_opt, ~] = fmincon(@(x)objective(x, x0, fk1, fk0, tri_new, tri_orig), x0, ...
    [], [], [], [], [], [], ...
    @(x)constraints(x, d1, d2, d3, l3, l4), options); 

% Create the new triangle(deformed)


% Define the constraint function
    function [c, ceq] = constraints(x, d1, d2, d3, l3, l4)
        % Extract vertex coordinates
        A = x(1:2)';
        B = x(3:4)';
        C = x(5:6)';

        % Rigidity constraints (triangle edge lengths preserved)
        ceq_rigidity = [
            norm(B - A)^2 - l4^2;  % AB
            norm(C - B)^2 - l4^2;  % BC
            norm(A - C)^2 - l4^2;  % CA
            ];

        % keep the same length of filaments l3
        ceq_distance = [
            norm(A - d3) - l3;  % A must be l3 from d3
            norm(B - d1) - l3;  % B must be l3 from d1
            norm(C - d2) - l3;  % C must be l3 from d2
            ];

        ceq = [ceq_rigidity; ceq_distance];
        c = [];  % No inequality constraints
    end

% Difine the objective function to minimize the energy
    function cost = objective(x, x0, fk1, fk0, tri_new, tri_orig)
        % vertices of deformed triangle
        current_triangle = reshape(x, 3, 2);
        A_current = current_triangle(1, :);
        B_current = current_triangle(2, :);
        C_current = current_triangle(3, :);

        % vertices of original triangle
        original_triangle = reshape(x0, 3, 2);
        A_orig = original_triangle(1, :);
        B_orig = original_triangle(2, :);
        C_orig = original_triangle(3, :);

        % vertices of deformed flanks
        current_flanks = reshape(fk1, 3, 2);
        dc1 = current_flanks(1, :);
        dc2 = current_flanks(2, :);    
        dc3 = current_flanks(3, :);

        % vertices of original flanks
        original_flanks = reshape(fk0, 3, 2);
        do1 = original_flanks(1, :);
        do2 = original_flanks(2, :);    
        do3 = original_flanks(3, :);

        % Calculate rotational energy, the bending stiffness is considered
        % as 1
        energy = 0;
        current_boundary = reshape(tri_new, 3, 2);
        p1_new = current_boundary(1,:);
        p2_new = current_boundary(2,:);
        p3_new = current_boundary(3,:);
        original_boundary = reshape(tri_orig, 3, 2);
        p1_orig = original_boundary(1,:);
        p2_orig = original_boundary(2,:);
        p3_orig = original_boundary(3,:);        

        % Filament 1
        vec1_orig = p1_orig - p2_orig;
        vec2_orig = B_orig - do1;
        vec3_orig = do1 - B_orig;
        vec4_orig = A_orig - B_orig;

        vec1_new = p1_new - p2_new;
        vec2_new = B_current - dc1;
        vec3_new = dc1 - B_current;
        vec4_new = A_current - B_current;

        angle1_orig = acos((vec1_orig * vec2_orig')/(norm(vec1_orig) * norm(vec2_orig)));
        angle2_orig = acos((vec3_orig * vec4_orig')/(norm(vec3_orig) * norm(vec4_orig)));

        angle1_new = acos((vec1_new * vec2_new')/(norm(vec1_new) * norm(vec2_new)));
        angle2_new = acos((vec3_new * vec4_new')/(norm(vec3_new) * norm(vec4_new)));

        energy = energy + 1/2*1*((angle1_new-angle1_orig)^2 + (angle2_new-angle2_orig)^2); 

        % Filament 2
        vec1_orig = p2_orig - p3_orig;
        vec2_orig = C_orig - do2;
        vec3_orig = do2 - C_orig;
        vec4_orig = B_orig - C_orig;

        vec1_new = p2_new - p3_new;
        vec2_new = C_current - dc2;
        vec3_new = dc2 - C_current;
        vec4_new = B_current - C_current;

        angle1_orig = acos((vec1_orig * vec2_orig')/(norm(vec1_orig) * norm(vec2_orig)));
        angle2_orig = acos((vec3_orig * vec4_orig')/(norm(vec3_orig) * norm(vec4_orig)));

        angle1_new = acos((vec1_new * vec2_new')/(norm(vec1_new) * norm(vec2_new)));
        angle2_new = acos((vec3_new * vec4_new')/(norm(vec3_new) * norm(vec4_new)));

        energy = energy + 1/2*1*((angle1_new-angle1_orig)^2 + (angle2_new-angle2_orig)^2); 

        % Filament 3
        vec1_orig = p3_orig - p1_orig;
        vec2_orig = A_orig - do3;
        vec3_orig = do3 - A_orig;
        vec4_orig = C_orig - A_orig;

        vec1_new = p3_new - p1_new;
        vec2_new = A_current - dc3;
        vec3_new = dc3 - A_current;
        vec4_new = C_current - A_current;

        angle1_orig = acos((vec1_orig * vec2_orig')/(norm(vec1_orig) * norm(vec2_orig)));
        angle2_orig = acos((vec3_orig * vec4_orig')/(norm(vec3_orig) * norm(vec4_orig)));

        angle1_new = acos((vec1_new * vec2_new')/(norm(vec1_new) * norm(vec2_new)));
        angle2_new = acos((vec3_new * vec4_new')/(norm(vec3_new) * norm(vec4_new)));

        energy = energy + 1/2*1*((angle1_new-angle1_orig)^2 + (angle2_new-angle2_orig)^2); 

        cost = energy;
    end
end






