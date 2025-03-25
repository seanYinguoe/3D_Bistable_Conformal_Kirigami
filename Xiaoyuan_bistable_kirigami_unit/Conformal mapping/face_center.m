function centroids = face_center(vertices, faces)
% Calculates the centroid of each face in a mesh.
%
% Inputs:
%   vertices - Nx3 matrix of vertex coordinates (each row is [x, y, z])
%   faces    - Mx3 (or MxN for N-sided polygons) matrix of face connectivity
%              (each row contains indices into the vertices matrix)
%
% Output:
%   centroids - Mx3 matrix of face centroids (each row is [cx, cy, cz])

% Initialize centroids matrix
num_faces = size(faces, 1);
centroids = zeros(num_faces, 2);

% Loop through each face
for i = 1:num_faces
    % Get the indices of the vertices for the current face
    vertex_indices = faces(i, :);
    
    % Remove any zero entries (for non-triangular faces with fewer vertices)
    vertex_indices = vertex_indices(vertex_indices > 0);
    
    % Extract the coordinates of the vertices for the current face
    face_vertices = vertices(vertex_indices, :);
    
    % Calculate the centroid of the current face
    centroids(i, :) = mean(face_vertices, 1); % Average of all vertex coordinates
end

end
