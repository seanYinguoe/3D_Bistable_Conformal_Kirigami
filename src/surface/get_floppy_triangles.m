function idE = get_floppy_triangles(f_grid)
%{
    Filter polygons that have at least one vertex that is not connected to any other triangles.
    
    f_grid - An Mx3 matrix where each row represents a triangle face, and each element in a row 
             is a vertex index.
    
    idE - An array of indices of faces that contain floppy vertices (vertices not shared by other faces).
%}

% Flatten the f_grid to get all vertex indices from the faces
f_grid_flat = f_grid(:);  % Convert the matrix into a single column vector

% Count the occurrences of each vertex in the flattened f_grid
[vertex_counts, unique_vertices] = hist(f_grid_flat, unique(f_grid_flat));

% Find vertices that appear only once (floppy vertices)
floppy_vertices = unique_vertices(vertex_counts == 1);

% Initialize the result array for faces that contain at least one floppy vertex
idE = [];

% Iterate through each face in f_grid to check if it contains any floppy vertices
for i = 1:size(f_grid, 1)
    % If any vertex in the current face is a floppy vertex, add the face index to idE
    if any(ismember(f_grid(i, :), floppy_vertices))
        idE = [idE; i];
    end
end
end
