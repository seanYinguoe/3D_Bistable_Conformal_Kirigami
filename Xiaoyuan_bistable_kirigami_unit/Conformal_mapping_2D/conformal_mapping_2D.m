function [v_initial_mesh, v_target_mesh, f_mesh, info_map] = conformal_mapping_2D(initialShape, targetShape, opts)
%CONFORMAL_MAPPING_2D Planar boundary-controlled mapping with fixed topology.
%   [v_initial_mesh, v_target_mesh, f_mesh, info_map] = conformal_mapping_2D(initialShape, targetShape, opts)
%
% Inputs
%   initialShape : shape descriptor (string or struct)
%                  strings: 'square', 'rectangle', 'circle', 'ellipse'
%                  struct:  .boundary (Nx2), or shape params (width/height/radius/a/b)
%   targetShape  : same format as initialShape
%   opts         : options struct
%                  .n_boundary (default 200)
%                  .mesh_h     (default auto from initial size)
%                  .initial_scale (default 100)
%                  .center     (default [0 0])
%
% Outputs
%   v_initial_mesh : Nv x3 initial planar mesh vertices (z=0)
%   v_target_mesh  : Nv x3 target planar mesh vertices (z=0)
%   f_mesh         : Nf x3 triangle connectivity
%   info_map       : diagnostics

if nargin < 1 || isempty(initialShape)
    initialShape = 'square';
end
if nargin < 2 || isempty(targetShape)
    targetShape = 'circle';
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'n_boundary') || isempty(opts.n_boundary)
    opts.n_boundary = 200;
end
if ~isfield(opts, 'initial_scale') || isempty(opts.initial_scale)
    opts.initial_scale = 100;
end
if ~isfield(opts, 'center') || isempty(opts.center)
    opts.center = [0 0];
end

% 1) Construct initial and target boundaries.
bnd_initial = shape_to_boundary(initialShape, opts.n_boundary, opts.initial_scale, opts.center);
bnd_target  = shape_to_boundary(targetShape,  opts.n_boundary, opts.initial_scale, opts.center);

% 2) Build initial triangular mesh inside initial boundary.
if ~isfield(opts, 'mesh_h') || isempty(opts.mesh_h)
    bb = [min(bnd_initial,[],1); max(bnd_initial,[],1)];
    span = bb(2,:) - bb(1,:);
    opts.mesh_h = max(min(span) / 25, 1e-3);
end
[v2, f_mesh] = mesh_inside_polygon(bnd_initial, opts.mesh_h);

% 3) Extract ordered boundary vertices of initial mesh.
[bnd_vid, bnd_loop] = ordered_boundary_vertices(v2, f_mesh);

% 4) Boundary correspondence by arclength; map to target boundary.
t = normalized_arclength(bnd_loop);
bnd_target_sampled = sample_closed_polyline(bnd_target, t);

% 5) Harmonic solve for interior vertices (Tutte/Laplacian style).
v2_target = harmonic_with_fixed_boundary(v2, f_mesh, bnd_vid, bnd_target_sampled);

% 6) Pack outputs as 3D-with-zero-z for compatibility.
v_initial_mesh = [v2, zeros(size(v2,1),1)];
v_target_mesh  = [v2_target, zeros(size(v2_target,1),1)];

info_map = struct();
info_map.n_vertices = size(v2,1);
info_map.n_faces = size(f_mesh,1);
info_map.n_boundary_vertices = numel(bnd_vid);
info_map.mesh_h = opts.mesh_h;
info_map.initial_boundary = bnd_initial;
info_map.target_boundary = bnd_target;
info_map.target_boundary_sampled = bnd_target_sampled;
end

function B = shape_to_boundary(shape, n, scale0, center)
if isnumeric(shape)
    if size(shape,2) ~= 2
        error('Numeric shape must be Nx2 boundary points.');
    end
    B = shape;
    return;
end

if ischar(shape) || (isstring(shape) && isscalar(shape))
    shape = struct('type', lower(char(shape)));
elseif ~isstruct(shape)
    error('shape must be string, struct, or Nx2 numeric boundary.');
end

if isfield(shape, 'boundary') && ~isempty(shape.boundary)
    if size(shape.boundary,2) ~= 2
        error('shape.boundary must be Nx2.');
    end
    B = shape.boundary;
    return;
end

if isfield(shape, 'type')
    type = lower(char(shape.type));
else
    type = 'polygon';
end

cx = center(1);
cy = center(2);

t = linspace(0, 1, n+1).';
t(end) = [];

switch type
    case 'square'
        L = field_or(shape, 'side', scale0);
        poly = [cx-L/2, cy-L/2;
                cx+L/2, cy-L/2;
                cx+L/2, cy+L/2;
                cx-L/2, cy+L/2];
        B = sample_closed_polyline(poly, t);

    case 'rectangle'
        W = field_or(shape, 'width', scale0);
        H = field_or(shape, 'height', 0.75*scale0);
        poly = [cx-W/2, cy-H/2;
                cx+W/2, cy-H/2;
                cx+W/2, cy+H/2;
                cx-W/2, cy+H/2];
        B = sample_closed_polyline(poly, t);

    case 'circle'
        R = field_or(shape, 'radius', 0.5*scale0);
        th = 2*pi*t;
        B = [cx + R*cos(th), cy + R*sin(th)];

    case 'ellipse'
        a = field_or(shape, 'a', 0.5*scale0);
        b = field_or(shape, 'b', 0.35*scale0);
        th = 2*pi*t;
        B = [cx + a*cos(th), cy + b*sin(th)];

    case 'polygon'
        if ~isfield(shape, 'vertices') || isempty(shape.vertices)
            error('For shape.type="polygon", provide shape.vertices (Nx2).');
        end
        poly = shape.vertices;
        if size(poly,2) ~= 2
            error('shape.vertices must be Nx2.');
        end
        B = sample_closed_polyline(poly, t);

    otherwise
        error('Unsupported shape type: %s', type);
end
end

function val = field_or(s, name, default_val)
if isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default_val;
end
end

function [V, F] = mesh_inside_polygon(poly, h)
minxy = min(poly, [], 1);
maxxy = max(poly, [], 1);

xv = (minxy(1):h:maxxy(1)).';
yv = (minxy(2):h:maxxy(2)).';
[XX, YY] = meshgrid(xv, yv);
P = [XX(:), YY(:)];

inside = inpolygon(P(:,1), P(:,2), poly(:,1), poly(:,2));
P_in = P(inside,:);

P_all = [poly; P_in];
P_all = unique(round(P_all / max(h,1e-9)) * max(h,1e-9), 'rows', 'stable');

DT = delaunayTriangulation(P_all);
F_all = DT.ConnectivityList;
V_all = DT.Points;

C = (V_all(F_all(:,1),:) + V_all(F_all(:,2),:) + V_all(F_all(:,3),:)) / 3;
keep = inpolygon(C(:,1), C(:,2), poly(:,1), poly(:,2));
F = F_all(keep,:);
V = V_all;

[V, F] = remove_unused_verts_local(V, F);
end

function [V2, F2] = remove_unused_verts_local(V, F)
used = false(size(V,1),1);
used(F(:)) = true;
map = zeros(size(V,1),1);
map(used) = 1:nnz(used);
V2 = V(used,:);
F2 = map(F);
end

function [bnd_vid, bnd_xy] = ordered_boundary_vertices(V, F)
TR = triangulation(F, V);
B = freeBoundary(TR);
if isempty(B)
    error('No free boundary found in initial mesh.');
end

adj = cell(size(V,1),1);
for e = 1:size(B,1)
    a = B(e,1);
    b = B(e,2);
    adj{a}(end+1) = b; %#ok<AGROW>
    adj{b}(end+1) = a; %#ok<AGROW>
end

start = B(1,1);
bnd_vid = start;
prev = 0;
curr = start;
for k = 1:size(B,1)+5
    nbr = adj{curr};
    if isempty(nbr)
        break;
    end
    if numel(nbr) == 1
        nxt = nbr(1);
    else
        cand = nbr(nbr ~= prev);
        if isempty(cand)
            nxt = nbr(1);
        else
            nxt = cand(1);
        end
    end
    if nxt == start
        break;
    end
    bnd_vid(end+1,1) = nxt; %#ok<AGROW>
    prev = curr;
    curr = nxt;
end

if numel(bnd_vid) < 3
    error('Failed to reconstruct ordered boundary loop.');
end

bnd_xy = V(bnd_vid,:);
end

function v_target = harmonic_with_fixed_boundary(V, F, bnd_idx, bnd_target)
n = size(V,1);
is_bnd = false(n,1);
is_bnd(bnd_idx) = true;

E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
E = sort(E,2);
E = unique(E, 'rows');

ii = [E(:,1); E(:,2)];
jj = [E(:,2); E(:,1)];
A = sparse(ii, jj, 1, n, n);
L = spdiags(sum(A,2), 0, n, n) - A;

bx = zeros(n,1);
by = zeros(n,1);
bx(bnd_idx) = bnd_target(:,1);
by(bnd_idx) = bnd_target(:,2);

M = L;
I = speye(n);
M(is_bnd,:) = I(is_bnd,:);

rhsx = zeros(n,1);
rhsy = zeros(n,1);
rhsx(~is_bnd) = -L(~is_bnd,is_bnd) * bx(is_bnd);
rhsy(~is_bnd) = -L(~is_bnd,is_bnd) * by(is_bnd);
rhsx(is_bnd) = bx(is_bnd);
rhsy(is_bnd) = by(is_bnd);

x = M \ rhsx;
y = M \ rhsy;
v_target = [x, y];
end

function t = normalized_arclength(P)
if size(P,1) < 2
    t = zeros(size(P,1),1);
    return;
end
Pc = [P; P(1,:)];
d = sqrt(sum(diff(Pc,1,1).^2,2));
s = [0; cumsum(d(1:end-1))];
L = sum(d);
if L <= eps
    t = linspace(0,1,size(P,1)+1).';
    t(end) = [];
else
    t = s / L;
end
end

function Q = sample_closed_polyline(poly, t)
if size(poly,1) < 2
    Q = repmat(poly(1,:), numel(t), 1);
    return;
end
polyC = [poly; poly(1,:)];
seg = sqrt(sum(diff(polyC,1,1).^2,2));
S = [0; cumsum(seg)];
L = S(end);
if L <= eps
    Q = repmat(poly(1,:), numel(t), 1);
    return;
end

ss = mod(t,1) * L;
Q = zeros(numel(t), 2);
for i = 1:numel(t)
    s = ss(i);
    k = find(S <= s, 1, 'last');
    if k >= numel(S)
        k = numel(S)-1;
    end
    ds = seg(k);
    if ds <= eps
        a = 0;
    else
        a = (s - S(k)) / ds;
    end
    Q(i,:) = (1-a)*polyC(k,:) + a*polyC(k+1,:);
end
end
