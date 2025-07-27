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


%% Calculate UV Mesh Face Centers
uv_c = face_center(v_mesh,f_mesh); % find the centroids of triangular mesh

%% Check regular grids validity(Remove the outside triangle)
tri = triangulation(f_mesh, v_mesh); % Define triangulation using given faces
grid_points = c_grid(:,[1,2]);
inside = ~isnan(tri.pointLocation(grid_points)); % Check if points are inside the mesh

% Replace all the outside centroid as NaN
f_grid(~inside, :) = NaN;
c_grid(~inside, :) = NaN;
i_grid(~inside, :) = NaN;

% Clean neighbor lists for all faces, keep the order unchanged
for xi = 1:numel(x_grid)
    if inside(xi)  % Only process valid faces
        % Remove neighbors that are invalid
        x_grid{xi} = x_grid{xi}(inside(x_grid{xi}));
        % If no neighbors left, mark as invalid
        if isempty(x_grid{xi})
            inside(xi) = true;
        end
    else
        x_grid{xi} = NaN;
    end
end


% Create index mapping
index_map = zeros(size(f_grid,1), 1);
index_map(inside) = (1:sum(inside))';  % Fixed dimension mismatch

% Update neighbor indices for x_grid
for xi = 1:numel(x_grid)
    if inside(xi)
        % Convert neighbors using index map, keep only valid ones
        x_grid{xi} = index_map(x_grid{xi}(inside(x_grid{xi})));
        x_grid{xi}(isnan(x_grid{xi})) = [];
    end
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

% Clean neighbor lists for all faces, keep the order unchanged
floppy_idx = false(size(f_grid,1),1);
floppy_idx(idE) = true;
non_floppy = ~floppy_idx;

for xi = 1:numel(x_grid)
    if non_floppy(xi)  % Only process valid faces
        % Remove neighbors that are invalid
        x_grid{xi} = x_grid{xi}(non_floppy(x_grid{xi}));
        % If no neighbors left, mark as invalid
        if isempty(x_grid{xi})
            non_floppy(xi) = true;
        end
    else
        x_grid{xi} = NaN;
    end
end

% Create index mapping
index_map = zeros(size(f_grid,1), 1);
index_map(non_floppy) = (1:sum(non_floppy))';  % Fixed dimension mismatch

% Update neighbor indices for x_grid
for xi = 1:numel(x_grid)
    if non_floppy(xi)
        % Convert neighbors using index map, keep only valid ones
        x_grid{xi} = index_map(x_grid{xi}(non_floppy(x_grid{xi})));
        x_grid{xi}(isnan(x_grid{xi})) = [];
    end
end

% Remove NaN and Adjust Indices
index = ~isnan(f_grid);
c_grid = c_grid(index(:,1),:);
i_grid = i_grid(index(:,1),:);
x_grid = x_grid(index(:,1),:);
f_grid = f_grid(index(:,1),:);

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
