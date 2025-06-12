function angles = dihedral_angle_calculate(v, f, x)
% CALCULATE_DIHEDRAL_ANGLES Computes angles between adjacent triangular faces
% Inputs:
%   v - Mx3 matrix of vertex coordinates
%   f - Nx3 matrix of face vertex indices
%   x - Nx1 cell array of adjacent face indices
% Output:
%   angles - Px3 matrix [face1, face2, angle_rads]

% Calculate face normals
normals = zeros(size(f,1), 3);
for i = 1:size(f,1)
    verts = v(f(i,:),:);
    v1 = verts(1,:);
    v2 = verts(2,:);
    v3 = verts(3,:);

    edge1 = v2 - v1;
    edge2 = v3 - v1;
    normal = cross(edge1, edge2);

    % Handle degenerate faces
    if norm(normal) < eps
        normals(i,:) = [0 0 0];
    else
        normals(i,:) = normal / norm(normal);
    end
end

% Track every adjacent faces, and calculate dihedral angle
angles = [];
for i = 1:size(f,1)
    % Get neighbors with j > i to prevent reciprocal duplication
    neighbors = x{i}(x{i} > i);
    
    for j = neighbors'
        n1 = normals(i,:);
        n2 = normals(j,:);
        
        % Calculate angle between normals
        cos_theta = dot(n1, n2);
        cos_theta = max(min(cos_theta, 1), -1); % Clamp for precision
        theta = acos(cos_theta);
        
        % Store dihedral angle 
        angles = [angles; i j theta];
    end
end
end