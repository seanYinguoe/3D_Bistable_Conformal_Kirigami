function [v_deployed, c_deployed] = conformal_grid_deployment(v_out, f_out, c_out, scale_area, v_target)
% CONFORMAL_3D_DEPLOYMENT Deploys 2D grid to 3D surface with conformal mapping
% Inputs:
%   v_out: Nx3 - Original 2D vertices (Z=0)
%   f_out: Mx3 - Face connectivity
%   c_out: Mx3 - Initial face centroids (Z=0)
%   scale_area: Mx1 - Area scaling factors
%   v_target: Nx3 - Target 3D vertices
% Outputs:
%   v_deployed: Nx3 - Deployed 3D vertices
%   c_deployed: Mx3 - Preserved centroids

%% Input Validation
assert(size(v_out,2) == 3, 'v_out must be Nx3 matrix');
assert(all(abs(v_out(:,3)) < eps), 'Vertices must lie on XY plane');
assert(size(f_out,2) == 3, 'Faces must be triangles');
assert(size(c_out,2) == 3, 'Centroids must be Mx3 matrix');
assert(all(scale_area > 0), 'Scale factors must be positive');

%% Preprocessing
% Build KD-tree for target surface projection
target_tree = KDTreeSearcher(v_target(:,1:2));
num_verts = size(v_out,1);
num_faces = size(f_out,1);

% Calculate original areas and centroids
original_areas = zeros(num_faces,1);
for i = 1:num_faces
    tri = v_out(f_out(i,:), :);
    original_areas(i) = polyarea(tri(:,1), tri(:,2));
end

%% Optimization Setup
options = optimoptions('fmincon',...
    'Algorithm','interior-point',...
    'MaxIterations',500,...
    'Display','iter',...
    'UseParallel',true);

% Initialize with original positions
x0 = double(v_out(:));

%% Solve Optimization Problem
[x_sol, ~] = fmincon(@conformal_energy, x0, [], [], [], [], [], [], @deployment_constraints, options);

%% Post-processing
v_deployed = reshape(x_sol, [], 3);
c_deployed = compute_centroids(v_deployed, f_out);

%% Nested Functions
    function E = conformal_energy(x)
        v_current = reshape(x, [], 3);
        E = 0;

        % Conformal distortion energy
        for i = 1:num_faces
            tri_orig = v_out(f_out(i,:), :);
            tri_def = v_current(f_out(i,:), :);

            % Compute deformation gradient
            J = (tri_def' - tri_def(1,:)') / (tri_orig' - tri_orig(1,:)' + eps);
            [~,S,~] = svd(J);

            % Conformal distortion metric
            distortion = (S(1,1)/S(2,2) - 1)^2 + (S(2,2)/S(3,3) - 1)^2;
            E = E + distortion;
        end

        % Target surface attraction
        [~, dists] = knnsearch(target_tree, v_current(:,1:2));
        E = E + 0.1*sum(dists.^2); % Weighting factor adjustable
    end

    function [c, ceq] = deployment_constraints(x)
        v_current = reshape(x, [], 3);
        ceq = zeros(num_faces*4,1); % x,y,z centroids + area

        % Face constraints
        for i = 1:num_faces
            verts = v_current(f_out(i,:), :);

            % Centroid preservation
            ceq(3*i-2:3*i) = mean(verts) - c_out(i,:);

            % Area scaling (3D area)
            current_area = mesh_face_area(verts);
            target_area = original_areas(i) * scale_area(i);
            ceq(3*num_faces + i) = (current_area - target_area)/target_area;
        end
        c = [];
    end

    function centroids = compute_centroids(vertices, faces)
        centroids = zeros(size(faces,1),3);
        for i = 1:size(faces,1)
            tri = vertices(faces(i,:), :);
            centroids(i,:) = mean(tri);
        end
    end

    function area = mesh_face_area(verts)
        % Calculates 3D face area
        v1 = verts(2,:) - verts(1,:);
        v2 = verts(3,:) - verts(1,:);
        area = 0.5*norm(cross(v1, v2));
    end
end
