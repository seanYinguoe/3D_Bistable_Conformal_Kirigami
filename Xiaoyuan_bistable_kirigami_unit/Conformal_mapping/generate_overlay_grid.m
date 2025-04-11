function [v_grid, f_grid, c_grid, i_grid, x_grid] = generate_overlay_grid(vt, edgeLen)
% Generates a polygon grid (triangular or square) that envelops the input mesh.
% Inputs:
%   vt   - Flattened mesh (Nx2 array of 2D points)
%   edgeLen  -  edge length of grid polygons
%
% Outputs:
%   v_grid  - Vertex coordinates of the grid
%   f_grid  - Faces (polygons) defined by vertex indices
%   c_grid  - Centroid coordinates of each polygon
%   i_grid  - Auxiliary data for polygons
%   x_grid  - Neighboring polygon IDs

factor = sqrt(3)/2;

% Compute bounding box of the input mesh-
boundv2d = [min(vt(:,1)), min(vt(:,2)); max(vt(:,1)), max(vt(:,2))];

% Dimensions of the bounding box
v2d_dim = [boundv2d(2,1) - boundv2d(1,1), boundv2d(2,2) - boundv2d(1,2)];

% Center of the bounding box
v2d_c = [v2d_dim(1)/2 + boundv2d(1,1), v2d_dim(2)/2 + boundv2d(1,2)];

% Number of columns and rows in the grid
cols = ceil(v2d_dim(1)/edgeLen+5);
rows = ceil(v2d_dim(2)/(factor*edgeLen)+5);

if rows + cols < 5
    warning('The design is very small. Check your model size and units (mm assumed).');
end

% Size of the grid and offset vector from center of bounding box to grid center
gc = [cols*edgeLen, rows*factor*edgeLen];
gd = [v2d_c(1) - gc(1)/2, v2d_c(2) - gc(2)/2];

% Generate the grid based on geometry type
[v_grid, f_grid, c_grid, i_grid, x_grid] = generate_triangular_grid(rows, cols, edgeLen, gd(1), gd(2));
end
