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

% Build a regular triangular grid.
% generate_triangular_grid(rows, cols) uses vertex counts.
%
% The face loop runs  for col = 0:cols-2, giving (cols-1) column steps.
% Each step produces 2 triangles (one up + one down), so one row-band has
%   2*(cols-1)  triangles in x.
%
% Strategy for any num_x (even or odd):
%   cols = ceil(num_x/2) + 1
%       Even num_x → 2*(cols-1) = num_x  (exact, no trimming needed)
%       Odd  num_x → 2*(cols-1) = num_x+1 (one extra per band → trim it)
%
% In the y direction: rows = num_y + 1  →  (rows-1) = num_y row-bands.
% Total triangles = num_x * num_y.
rows = max(2, num_y + 1);
cols = max(2, ceil(num_x/2) + 1);
[v_grid, f_grid, ~, i_grid, ~] = generate_triangular_grid(rows, cols, edgeLen, 0, 0);

% Anchor the grid at the origin.
v_grid(:,1) = v_grid(:,1) - min(v_grid(:,1));
v_grid(:,2) = v_grid(:,2) - min(v_grid(:,2));

% For odd num_x the grid has num_x+1 triangles per row-band; remove the
% last face of every row-band (the rightmost triangle) to get exactly num_x.
n_bands       = rows - 1;
faces_per_band = 2 * (cols - 1);          % num_x+1 when num_x is odd
if mod(num_x, 2) ~= 0
    idx_last = (1:n_bands)' * faces_per_band;  % last face index in each band
    keep = true(size(f_grid, 1), 1);
    keep(idx_last) = false;
    f_out = f_grid(keep, :);
    i_out = i_grid(keep);
else
    f_out = f_grid;
    i_out = i_grid;
end

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
