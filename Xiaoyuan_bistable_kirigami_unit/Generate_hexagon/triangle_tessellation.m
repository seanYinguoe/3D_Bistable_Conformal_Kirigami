function tessellation = triangle_tessellation(num_x, num_y, beta, edgeLen, l1, l4, t)
%TRIANGLE_TESSELLATION Build a rectangular sheet using the existing
%triangular-grid and conformal-mapping tessellation pipeline.
%
%   tessellation = triangle_tessellation(num_x, num_y, beta, edgeLen, l1, l4, t)
%
% Workflow
%   1) Create a rectangular triangular grid.
%   2) Use the existing tessellated_triangle_initial function to stamp one
%      kirigami triangle unit on each grid face.
%
% Notes
%   - num_x and num_y control the rectangle size in grid cells.
%   - beta can be a scalar or one value per triangle face.

if nargin < 7
    error('triangle_tessellation requires num_x, num_y, beta, edgeLen, l1, l4, and t.');
end

% Rectangle size implied by the requested unit counts.
length_x = edgeLen * num_x;
length_y = edgeLen * num_y;

% Build a regular triangular grid over the rectangle.
% generate_triangular_grid uses rows/cols of vertices, so we add one.
rows = max(2, num_y + 1);
cols = max(2, num_x + 1);
[v_grid, f_grid, ~, i_grid, ~] = generate_triangular_grid(rows, cols, edgeLen, 0, 0);

% Center the grid so the rectangle is anchored cleanly in XY.
v_grid(:,1) = v_grid(:,1) - min(v_grid(:,1));
v_grid(:,2) = v_grid(:,2) - min(v_grid(:,2));

% Keep only triangles whose centroids lie inside the desired rectangle.
c_tri = face_center(v_grid(:,1:2), f_grid);
keep = c_tri(:,1) >= 0 & c_tri(:,1) <= length_x & ...
       c_tri(:,2) >= 0 & c_tri(:,2) <= length_y;

f_out = f_grid(keep, :);
i_out = i_grid(keep);

% Use the same parameter packing as the conformal-mapping pipeline.
params = [edgeLen; l1; l4; t];

% beta may be scalar or per-face.
if isscalar(beta)
    beta_use = repmat(beta, size(f_out,1), 1);
else
    beta_use = beta(:);
    if numel(beta_use) ~= size(f_out,1)
        error('beta must be a scalar or have one value per retained triangle face.');
    end
end

tessellation = tessellated_triangle_initial(f_out, i_out, params, v_grid, beta_use, t);
end
