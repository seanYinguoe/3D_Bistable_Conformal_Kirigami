function [v_out, f_out, c_out, i_out, x_out, scale_area] = fit_grid(v_grid, f_grid, c_grid, i_grid, x_grid, v_mesh,f_mesh, scale_facs, edgeLen)
% Fits a regular triangular grid to a flattened triangular mesh and interpolates scale factors.
%
% Inputs:
%   v_grid     - Vertex coordinates of the triangular grid
%   f_grid     - Face connectivity information (triangles)
%   c_grid     - Centroid coordinates for each triangle in the grid
%   i_grid     - Orientation indicator for each triangle (up/down)
%   x_grid     - Neighboring triangles information
%   uv_mesh    - Flattened triangular mesh (UV map: obj_2D.vt, obj_2D.f.vt)
%   scale_facs - Scale factors calculated per UV triangle
%   edgeLen    - Edge length of polygons in the grid
%
% Outputs:
%   v_out      - Adjusted vertex coordinates
%   f_out      - Adjusted face connectivity
%   c_out      - Centroids of valid triangles
%   i_out      - Orientations of valid triangles
%   x_out      - Neighboring triangles relationships
%   scale_area - Interpolated scale factors for each triangle
%   uv_c       - Centroids of UV mesh faces


%% Calculate UV Mesh Face Centers
uv_c = face_center(v_mesh,f_mesh); % find the centroids of triangular mesh

%% Check regular grids validity
tri = triangulation(f_mesh, v_mesh); % Define triangulation using given faces
grid_points = c_grid(:,[1,2]);
inside = ~isnan(tri.pointLocation(grid_points)); % Check if points are inside the mesh

% Replace all the outside centroid as NaN
f_grid(~inside, :) = NaN;
c_grid(~inside, :) = NaN;
i_grid(~inside, :) = NaN;
for i = find(~inside)'
    x_grid{i} = NaN;
end

% Remove NaN
index = ~isnan(f_grid);
c_grid = c_grid(index(:,1),:);
i_grid = i_grid(index(:,1),:);
x_grid = x_grid(index(:,1),:);
f_grid = f_grid(index(:,1),:);

%% Identify "Floppy" Polygons, find the triangles that are detached from other triangles
idE = get_floppy_triangles(f_grid); 

[f_grid, idxn] = removeId(f_grid, idE); % Remove floppy polygons, replaces as NaN
c_grid(idxn,:) = NaN;
i_grid(idxn,:) = NaN;
for i = find(idxn)'
    x_grid{i} = NaN;
end

% Remove NaN and Adjust Indices
index = ~isnan(f_grid);
c_grid = c_grid(index(:,1),:);
i_grid = i_grid(index(:,1),:);
x_grid = x_grid(index(:,1),:);
f_grid = f_grid(index(:,1),:);

% Iterate over each element of x_grid (each cell represents a list of neighboring faces)
for xi = 1:numel(x_grid)
    current_neighbors = x_grid{xi};
    updated_neighbors = current_neighbors(~ismember(current_neighbors, idxn));
    x_grid{xi} = updated_neighbors;
end

%% Calculate Scale Factors Using Weighted Averaging
scale_area = zeros(size(c_grid,1),1);

for ci_idx = 1:size(c_grid,1)
    beta = 0.5;
    neighbors_idx = pointCloudClosestPoints(uv_c, c_grid(ci_idx,1:2), beta * edgeLen); 
    weights = normalizeWeights(...
    -1/(beta*edgeLen) * vecnorm(c_grid(ci_idx,1:2) - uv_c(neighbors_idx,1:2), 2, 2) + 1 ...
);
    scale_area(ci_idx) = sum(scale_facs(neighbors_idx) .* weights); % Weighted average of scale factors
end

%% Final Adjustments and Output Preparation
v_out = v_grid; % Vertex coordinates remain unchanged
f_out = f_grid; % Face connectivity remains unchanged
c_out = c_grid; % Centroids remain unchanged after filtering invalid entries
i_out = i_grid; % Orientation remains unchanged after filtering invalid entries
x_out = x_grid; % Neighboring relationships remain unchanged after filtering invalid entries

[v_out, f_out] = remove_unused_verts(v_out, f_out); 

end
