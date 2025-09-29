function [triangle_new,E_total] = deform_triangle_isotropic(delta,edgeLen,l1,l4,beta,t)
% DEFORM_TRIANGLE Deforms a isotropic triangle based on input node positions and edge lengths.
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
%   triangle_new - Deformed triangle coordinates

%% Define local function
rotation = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)]; % rotation matrix  

%% Calculate geometric parameters
A = pi/3-beta;
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(A); % length of filament
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) ...
     - cos(pi/3 - beta) .* l4;
l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 ...
     - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) );

%% Generate uniformaly deployed triangle as initial guess
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;
[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);
triangle_new = triangle;

%% Define optimise parameters
% Nodes in flanks
d1 = triangle(20,:);
d1_ = triangle(19,:);

% Extract original inner triangle vertices as initial guess
B_orig = triangle(44,:);
x0 = B_orig;

edge = edgeLen + delta; % deformed length of a unit
centroid = [-sqrt(3)/6*edge,-1/2*edge];

% Define the optimization conditions
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'iter', ...
    'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12, ...
    'ConstraintTolerance', 1e-10);

% Run optimization
[x_opt, ~] = fmincon(@(x)objective_energy(x, d1, d1_, l2, t, centroid, rotation), x0, ...
    [], [], [], [], [], [], ...
    @(x)constraints(x, d1, d1_, centroid,l3), options);


E_total = objective_energy(x_opt, d1, d1_, l2, t, centroid, rotation);

%% Create optimised unit
% Create the new inner triangle
B_new = x_opt;
A_new = x_opt * rotation(-2*pi/3);
C_new = x_opt * rotation(2*pi/3);

triangle_new(43,:) = A_new;
triangle_new(44,:) = B_new;
triangle_new(45,:) = C_new;

% Create the new filaments
triangle_new(32,:) = triangle_new(20,:);
triangle_new(33,:) = triangle_new(44,:);
triangle_new(34,:) = triangle_new(33,:) + t/(norm(A_new-B_new))*(A_new-B_new);

triangle_new(37,:) = triangle_new(45,:);
triangle_new(36,:) = triangle_new(24,:);
triangle_new(38,:) = triangle_new(37,:) + t/(norm(B_new-C_new))*(B_new-C_new);

triangle_new(40,:) = triangle_new(28,:);
triangle_new(41,:) = triangle_new(43,:);
triangle_new(42,:) = triangle_new(41,:) + t/(norm(C_new-A_new))*(C_new-A_new);
end

%% Define the constraint function
function [c, ceq] = constraints(x, d1, d1_,centroid,l3) 
% Extract vertex coordinates
B = x;

% Rigidity constraints (keep the triangle rigid)
ceq_rigidity = norm(B - centroid) - sqrt(3)/3*l3;

ceq = ceq_rigidity;

% No-overlapping condition
nodes1 = [d1_;d1;B];
res1 = ifoverlapping(nodes1); % if nonoverlapping, res<0
c = res1;  % No inequality constraints
end



%% Difine the objective function to minimize the energy(Hencky bar-chain model)
function E_total = objective_energy(x, dc1, dc1_, l2, t, centroid, rotation)
% Define the stretch stiffness and bend stiffness
E = 4.33e11; % large E could pollute the results
b = 1;  % the width of the sheet
N = 10; % the number of segment in HBM

% vertices of current triangle
B_current = x;
A_current = (B_current-centroid)*rotation(-2*pi/3);


vec1 = dc1_ - dc1;
vec2 = A_current - B_current;

%% Calculate bending energy
[E_total, ~, ~, ~] = hbm_energy(l2, dc1, B_current, vec1, vec2, N, E, b, t,...
    'VectorsAreNormals', true);

E_total = 3 * E_total;
end