function [v_out, f_out, c_out, i_out, x_out, UC_param, uv_scale, scale_length, errorTog] = interp_and_rescale(v_grid, f_grid, c_grid, i_grid, x_grid, scale_area, params, verbose)
% Interpolates and rescales a triangular grid based on scale factors from a flattened mesh.
% Inputs:
%   v_grid      - Vertex coordinates of the triangular grid
%   f_grid      - Face connectivity information (triangles)
%   c_grid      - Centroid coordinates for each triangle
%   i_grid      - Orientation indicator for each triangle (up/down)
%   x_grid      - Neighboring triangles information
%   scale_area  - Scale areas computed per triangle
%   params      - Struct containing parameters:
%                 params.scaleShift ('max' or other)
%                 params.scale_maximal (maximum scaling factor)
%                 params.scale_minimal (minimum scaling factor)
%                 params.DBPath (path to database file)
%                 params.edgeLen (edge length parameter)
%                 params.selectMode (selection mode for stable expansions)
%   verbose     - Boolean flag to print detailed information (default: true)

if nargin < 8
    verbose = true;
end

% Convert area to length
scale_length = sqrt(scale_area);

if verbose
    fprintf("Min factor (uv): %.4f\n", min(scale_length));
    fprintf("Max factor (uv): %.4f\n", max(scale_length));
end

% Determine scaling factor (uv_scale)
if strcmp(params.scaleShift, 'max')
    % Scale so that the largest scale matches the maximum allowed scale
    uv_scale = params.scale_maximal / max(scale_length);
else
    % Scale so that the smallest scale matches the minimum allowed scale
    uv_scale = params.scale_minimal / min(scale_length);
end

% Step 3: Apply scaling factor
scale_length = scale_length * uv_scale;

if verbose
    fprintf("Scaled by %.4f to:\n", uv_scale);
    fprintf("Scaled min factor (uv): %.4f\n", min(scale_length));
    fprintf("Scaled max factor (uv): %.4f\n", max(scale_length));
    fprintf("\n");
end

% Step 4: Check for errors in scaling limits
errorTog = false;
if max(scale_length) > params.scale_maximal * 1.05
    warning("Design exceeds maximal length scaling.");
    errorTog = true;
end
if min(scale_length) < params.scale_minimal / 1.05
    warning("Design exceeds minimal length scaling.");
    errorTog = true;
end

% Step 5: Retrieve parameters from database
UC_param = [];
data = readData(params.DBPath, params.edgeLen); % External helper function to read data
id_selected = get_most_stable_for_expansion(scale_length, data, params.selectMode); % External helper function

if ~isempty(id_selected)
    for idi = id_selected(:)'
        UC_param(end+1,:) = cell2mat(arrayfun(@(i) data{i}(idi), 1:5, 'UniformOutput', false));
    end
end

% Step 6: Prepare output variables
v_out = v_grid; % Vertex coordinates remain unchanged
f_out = f_grid; % Face connectivity remains unchanged
c_out = c_grid; % Centroids remain unchanged
i_out = i_grid; % Orientation remains unchanged
x_out = x_grid; % Neighboring triangles remain unchanged

% Step 7: Remove unused vertices from geometry
[v_out, f_out] = remove_unused_verts(v_out, f_out); % External helper function

end

