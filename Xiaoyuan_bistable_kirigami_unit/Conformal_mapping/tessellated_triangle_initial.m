function triangle_tessellation = tessellated_triangle_initial(f_out, i_out, params, v, beta, varargin)
%TESSELLATED_TRIANGLE_INITIAL Build initial (delta=0) unit on each grid cell.
%   This version does NOT run deform_triangle optimization and does NOT
%   move flanks independently. It uses triangle_unit at delta=0 and maps
%   the full unit to each target triangle by barycentric coordinates.
%
% Inputs
%   f_out  : Mx3 grid connectivity
%   i_out  : Mx1 orientation flag (0 up, 1 down)
%   params : [edgeLen; l1; l4; t_default]
%   v      : Nv x3 deployed/planar grid vertices
%   beta   : Mx1 beta per cell
%   varargin{1} (optional): scalar t or Mx1 t per cell
%
% Output
%   triangle_tessellation : cell(M,1), each cell is 45x3 unit vertices

edgeLen = params(1);
l1 = params(2);
l4 = params(3);
t = params(4);

if nargin >= 6 && ~isempty(varargin{1})
    t = varargin{1};
end

if isscalar(t)
    t = repmat(t, size(f_out,1), 1);
else
    t = t(:);
    if numel(t) ~= size(f_out,1)
        error('t must be a scalar or an Nx1 vector with one value per unit.');
    end
end

if size(v,2) == 2
    v = [v, zeros(size(v,1),1, 'like', v)];
end

delta = 0;
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;

triangle_tessellation = cell(size(f_out,1),1);
for i = 1:size(f_out,1)
    % 1) Build the canonical initial unit in local 2D coordinates.
    tri_local = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta(i), edgeLen, l1, l4, t(i));
    tri_local3 = [tri_local, zeros(size(tri_local,1),1, 'like', tri_local)];

    % 3) Local boundary used for barycentric mapping.
    boundary_local = [tri_local3(22,:); tri_local3(26,:); tri_local3(30,:)];

    % 4) Target grid triangle (same up/down ordering logic as tessellation).
    q1 = v(f_out(i,1),:);
    q2 = v(f_out(i,2),:);
    q3 = v(f_out(i,3),:);
    if i_out(i) == 0 % upwards triangle
        boundary_target = [q3; q1; q2];
    else % downwards triangle
        boundary_target = [q1; q3; q2];
    end

    % 5) Map every local unit point to target triangle.
    bc = zeros(size(tri_local3,1), 3, 'like', tri_local3);
    for j = 1:size(tri_local3,1)
        bc(j,:) = cart2barycentric(boundary_local, tri_local3(j,:));
    end

    triangle_out = bc * boundary_target;
    triangle_tessellation{i} = triangle_out;
end

end
