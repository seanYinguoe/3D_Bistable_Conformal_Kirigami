function [output_path, info] = generate_svg(tessellation, varargin)
%GENERATE_SVG Export kirigami cut/engrave pattern from topology.
%   output_path = generate_svg(tessellation)
%   output_path = generate_svg(tessellation, filename)
%   output_path = generate_svg(tessellation, filename, add_engrave)
%   output_path = generate_svg(tessellation, filename, add_engrave, snap_tol)
%   output_path = generate_svg(tessellation, filename, add_engrave, snap_tol, fillet_radius)
%
% Output contains:
%   1) global outer boundary loops (from single-use cell outer edges)
%   2) void loops (direct from f_void for each cell)
%
% No point-cloud boundary()/convhull() is used.

if nargin < 1 || ~iscell(tessellation)
    error('tessellation must be a cell array.');
end

% ========================= INPUT PARSING =========================
[filename, add_engrave, snap_tol, fillet_radius] = resolve_inputs(varargin{:});
func_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(func_dir, 'output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_path = fullfile(output_dir, filename);

% Known per-cell face indexing
f_void = [1 2 3 4 5 6;
          7 8 9 10 11 12;
          13 14 15 16 17 18];

% Per-cell outer boundary vertex ids (triangle corners in unit indexing)
% (from deform_triangle mapping: 22->p1, 26->p2, 30->p3)
f_outer = [22 26; 26 30; 30 22];

% Collect all coordinates for scale
all_xy = zeros(0,2);
for i = 1:numel(tessellation)
    V = tessellation{i};
    if isempty(V), continue; end
    if size(V,2) < 2
        continue;
    end
    all_xy = [all_xy; V(:,1:2)]; %#ok<AGROW>
end
if isempty(all_xy)
    error('No valid XY data in tessellation.');
end

if isempty(snap_tol)
    span = max(max(all_xy,[],1) - min(all_xy,[],1));
    snap_tol = max(1e-8, 1e-6 * max(span, 1));
else
    snap_tol = max(1e-12, abs(snap_tol));
end

% Global snapped point table + edge counting
point_map = containers.Map('KeyType','char', 'ValueType','uint32');
sum_xy = zeros(0,2);
cnt_xy = zeros(0,1);

edge_map_outer = containers.Map('KeyType','char', 'ValueType','uint32');
edge_tbl_outer = zeros(0,3,'uint32'); % [count, i, j]
edge_map_void = containers.Map('KeyType','char', 'ValueType','uint32');
edge_tbl_void = zeros(0,3,'uint32'); % [count, i, j]

n_invalid_cells = 0;

for i = 1:numel(tessellation)
    V = tessellation{i};
    if isempty(V) || size(V,2) < 2 || size(V,1) < 30
        n_invalid_cells = n_invalid_cells + 1;
        continue;
    end
    XY = V(:,1:2);

    % 1) Void edges (topology). Shared edges are removed later.
    for r = 1:size(f_void,1)
        idx = f_void(r,:);
        P = XY(idx,:);
        P = sanitize_points(P, snap_tol);
        if size(P,1) < 2
            continue;
        end
        ids = zeros(size(P,1),1,'uint32');
        for m = 1:size(P,1)
            [ids(m), point_map, sum_xy, cnt_xy] = register_point(P(m,:), snap_tol, point_map, sum_xy, cnt_xy);
        end
        for m = 1:numel(ids)
            a = ids(m);
            b = ids(mod(m, numel(ids)) + 1);
            if a == b
                continue;
            end
            [edge_map_void, edge_tbl_void] = add_undirected_edge(a, b, edge_map_void, edge_tbl_void);
        end
    end

    % 2) Outer edges for topology counting
    for e = 1:size(f_outer,1)
        aLoc = f_outer(e,1);
        bLoc = f_outer(e,2);
        if aLoc > size(XY,1) || bLoc > size(XY,1)
            continue;
        end
        pa = XY(aLoc,:);
        pb = XY(bLoc,:);
        if ~all(isfinite(pa)) || ~all(isfinite(pb))
            continue;
        end
        [idA, point_map, sum_xy, cnt_xy] = register_point(pa, snap_tol, point_map, sum_xy, cnt_xy);
        [idB, point_map, sum_xy, cnt_xy] = register_point(pb, snap_tol, point_map, sum_xy, cnt_xy);
        if idA == idB
            continue;
        end
        [edge_map_outer, edge_tbl_outer] = add_undirected_edge(idA, idB, edge_map_outer, edge_tbl_outer);
    end
end

if isempty(sum_xy)
    error('No valid outer-boundary edges extracted.');
end
pts = sum_xy ./ max(cnt_xy, 1);

% ======================= CUT EDGE EXPORT ========================
% Keep only single-use outer edges (global boundary edges)
outer_edges = single_use_edges(edge_map_outer, edge_tbl_outer);
if isempty(outer_edges)
    error('No single-use outer edges found for global boundary.');
end
boundary_vertex_mask = false(size(pts,1),1);
boundary_vertex_mask(double(unique(outer_edges(:)))) = true;
outer_edges_d = double(outer_edges);
boundary_segA = pts(outer_edges_d(:,1), :);
boundary_segB = pts(outer_edges_d(:,2), :);

% Reconstruct ordered boundary loops
boundary_loops = reconstruct_loops(outer_edges);
boundary_loops = filter_short_loops(boundary_loops);
if isempty(boundary_loops)
    error('Failed to reconstruct global boundary loop(s).');
end

% Convert boundary loops to SVG paths
boundary_paths = cell(numel(boundary_loops),1);
for i = 1:numel(boundary_loops)
    ids = boundary_loops{i};
    P = pts(ids, :);
    boundary_paths{i} = points_to_closed_path(P);
end

% Keep only single-use void edges (shared/internal void edges removed)
void_edges = single_use_edges(edge_map_void, edge_tbl_void);
void_paths = {};
if ~isempty(void_edges)
    % Vertex-based path reconstruction from edge graph.
    [void_path_ids, void_path_closed, deg_void] = reconstruct_paths_from_edges(void_edges);
    void_paths = cell(numel(void_path_ids),1);
    for i = 1:numel(void_path_ids)
        ids = void_path_ids{i};
        P = pts(ids, :);
        if fillet_radius > 0
            % Proper corner fillets only at valid degree-2 vertices.
            void_paths{i} = points_to_path_with_vertex_fillets( ...
                P, ids, void_path_closed(i), deg_void, boundary_vertex_mask, ...
                boundary_segA, boundary_segB, snap_tol, fillet_radius);
        else
            if void_path_closed(i)
                void_paths{i} = points_to_closed_path(P);
            else
                void_paths{i} = points_to_open_path(P);
            end
        end
    end
end

% ==================== ENGRAVE EDGE DETECTION ====================
% Engrave edges are shared panel/flank boundaries across neighboring units:
% use multi-use outer edges (count >= 2), then remove any edge that is
% already exported as a cut edge.
engrave_edges = zeros(0,2,'uint32');
engrave_paths = {};
if add_engrave
    shared_outer_edges = multi_use_edges(edge_map_outer, edge_tbl_outer, 2);
    cut_keys = build_edge_key_set([outer_edges; void_edges]);
    keep = true(size(shared_outer_edges,1),1);
    for i = 1:size(shared_outer_edges,1)
        if isKey(cut_keys, edge_key(shared_outer_edges(i,1), shared_outer_edges(i,2)))
            keep(i) = false;
        end
    end
    engrave_edges = shared_outer_edges(keep,:);
    engrave_paths = edges_to_path_segments(engrave_edges, pts);
end

% SVG canvas from all data points
xmin = min(all_xy(:,1)); xmax = max(all_xy(:,1));
ymin = min(all_xy(:,2)); ymax = max(all_xy(:,2));
span2 = max([xmax-xmin, ymax-ymin, 1e-6]);
pad = 0.02 * span2;
vb = [xmin-pad, ymin-pad, (xmax-xmin)+2*pad, (ymax-ymin)+2*pad];

fid = fopen(output_path, 'w');
if fid < 0
    error('Cannot open %s for writing.', output_path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid, '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" ');
fprintf(fid, 'viewBox="%.9g %.9g %.9g %.9g">\n', vb(1), vb(2), vb(3), vb(4));
fprintf(fid, '  <title>Kirigami cut pattern</title>\n');
fprintf(fid, '  <desc>Topology-based cut paths with optional engrave paths.</desc>\n');

% ======================== CUT PATH EXPORT ========================
fprintf(fid, '  <g id="cut" fill="none" stroke="#FF0000" stroke-width="0.05" stroke-linecap="round" stroke-linejoin="round">\n');
for i = 1:numel(boundary_paths)
    fprintf(fid, '    <path d="%s" />\n', boundary_paths{i});
end
for i = 1:numel(void_paths)
    fprintf(fid, '    <path d="%s" />\n', void_paths{i});
end
fprintf(fid, '  </g>\n');

if add_engrave
    % ====================== ENGRAVE EDGE EXPORT =====================
    fprintf(fid, '  <g id="engrave" fill="none" stroke="#0000FF" stroke-width="0.05" stroke-linecap="round" stroke-linejoin="round">\n');
    for i = 1:numel(engrave_paths)
        fprintf(fid, '    <path d="%s" />\n', engrave_paths{i});
    end
    fprintf(fid, '  </g>\n');
end

fprintf(fid, '</svg>\n');

info = struct();
info.output_path = output_path;
info.n_boundary_loops = numel(boundary_paths);
info.n_void_paths = numel(void_paths);
info.n_engrave_paths = numel(engrave_paths);
info.n_outer_edges = size(outer_edges,1);
info.n_engrave_edges = size(engrave_edges,1);
info.n_invalid_cells = n_invalid_cells;
info.snap_tol = snap_tol;
info.add_engrave = add_engrave;
info.fillet_radius = fillet_radius;
info.viewBox = vb;
end

function [id, point_map, sum_xy, cnt_xy] = register_point(p, tol, point_map, sum_xy, cnt_xy)
key = snap_key(p, tol);
if isKey(point_map, key)
    id = point_map(key);
else
    id = uint32(size(sum_xy,1) + 1);
    point_map(key) = id;
    sum_xy(end+1,:) = [0, 0];
    cnt_xy(end+1,1) = 0;
end
sum_xy(id,:) = sum_xy(id,:) + p;
cnt_xy(id,1) = cnt_xy(id,1) + 1;
end

function [edge_map, edge_tbl] = add_undirected_edge(a, b, edge_map, edge_tbl)
if a < b
    i = a; j = b;
else
    i = b; j = a;
end
key = edge_key(i, j);
if isKey(edge_map, key)
    idx = edge_map(key);
    edge_tbl(idx,1) = edge_tbl(idx,1) + 1;
else
    idx = uint32(size(edge_tbl,1) + 1);
    edge_map(key) = idx;
    edge_tbl(idx,:) = [1, i, j];
end
end

function E = single_use_edges(edge_map, edge_tbl)
keys = edge_map.keys;
E = zeros(0,2,'uint32');
for i = 1:numel(keys)
    idx = edge_map(keys{i});
    if edge_tbl(idx,1) == 1
        E(end+1,:) = edge_tbl(idx,2:3); %#ok<AGROW>
    end
end
end

function E = multi_use_edges(edge_map, edge_tbl, min_count)
if nargin < 3 || isempty(min_count)
    min_count = 2;
end
keys = edge_map.keys;
E = zeros(0,2,'uint32');
for i = 1:numel(keys)
    idx = edge_map(keys{i});
    if edge_tbl(idx,1) >= min_count
        E(end+1,:) = edge_tbl(idx,2:3); %#ok<AGROW>
    end
end
end

function loops = reconstruct_loops(E)
loops = {};
if isempty(E), return; end

E = double(E);
nE = size(E,1);
maxV = max(E(:));
adj = cell(maxV,1);
for e = 1:nE
    a = E(e,1); b = E(e,2);
    adj{a}(end+1) = e;
    adj{b}(end+1) = e;
end

used = false(nE,1);
for e0 = 1:nE
    if used(e0), continue; end
    used(e0) = true;
    a = E(e0,1); b = E(e0,2);
    loop = [a; b];
    prev = a;
    curr = b;

    while true
        cand = adj{curr};
        nextEdge = 0;
        nextV = 0;
        for k = 1:numel(cand)
            ee = cand(k);
            if used(ee), continue; end
            u = E(ee,1); v = E(ee,2);
            nv = u + v - curr;
            if nv == prev
                continue;
            end
            nextEdge = ee;
            nextV = nv;
            break;
        end
        if nextEdge == 0
            break;
        end
        used(nextEdge) = true;
        loop(end+1,1) = nextV; %#ok<AGROW>
        prev = curr;
        curr = nextV;
        if curr == loop(1)
            break;
        end
    end

    loop = clean_loop(loop);
    if numel(loop) >= 3
        loops{end+1} = loop; %#ok<AGROW>
    end
end
end

function loops = filter_short_loops(loops)
if isempty(loops), return; end
keep = true(numel(loops),1);
for i = 1:numel(loops)
    if numel(loops{i}) < 3
        keep(i) = false;
    end
end
loops = loops(keep);
end

function P = sanitize_points(P, tol)
if isempty(P) || size(P,2) ~= 2
    P = [];
    return;
end
P = P(all(isfinite(P),2), :);
if isempty(P), return; end

if tol > 0
    P = round(P / tol) * tol;
end

if size(P,1) >= 2
    d = hypot(diff(P(:,1)), diff(P(:,2)));
    thr = max(tol, 1e-12);
    keep = [true; d > thr];
    P = P(keep,:);
end
end

function loop = clean_loop(loop)
if isempty(loop), return; end
keep = [true; diff(loop) ~= 0];
loop = loop(keep);
if numel(loop) >= 2 && loop(1) == loop(end)
    loop(end) = [];
end
end

function d = points_to_closed_path(P)
d = sprintf('M %.9g %.9g', P(1,1), P(1,2));
for i = 2:size(P,1)
    d = sprintf('%s L %.9g %.9g', d, P(i,1), P(i,2)); %#ok<AGROW>
end
d = sprintf('%s Z', d);
end

function d = points_to_open_path(P)
d = sprintf('M %.9g %.9g', P(1,1), P(1,2));
for i = 2:size(P,1)
    d = sprintf('%s L %.9g %.9g', d, P(i,1), P(i,2)); %#ok<AGROW>
end
end

function d = points_to_path_with_vertex_fillets(P, ids, is_closed, deg_map, boundary_vertex_mask, boundary_segA, boundary_segB, snap_tol, r_in)
% Replace each valid corner by exactly one tangent circular fillet arc.
n = size(P,1);
if n < 2 || ~isfinite(r_in) || r_in <= 0
    if is_closed
        d = points_to_closed_path(P);
    else
        d = points_to_open_path(P);
    end
    return;
end

hasF = false(n,1);
Tin  = nan(n,2);
Tout = nan(n,2);
rEff = nan(n,1);
sweep = zeros(n,1);

for i = 1:n
    if ~is_closed && (i == 1 || i == n)
        continue; % never fillet open endpoints
    end
    vid = ids(i);
    if vid < 1 || vid > numel(deg_map) || deg_map(vid) ~= 2
        continue; % only corners shared by exactly two edges
    end
    if vid <= numel(boundary_vertex_mask) && boundary_vertex_mask(vid)
        continue; % do not fillet any corner touching global outer boundary
    end
    % Also block fillets for points lying on boundary segments
    if is_point_on_any_segment(P(i,:), boundary_segA, boundary_segB, max(snap_tol, 1e-8))
        continue;
    end

    im1 = i - 1; ip1 = i + 1;
    if is_closed
        if im1 < 1, im1 = n; end
        if ip1 > n, ip1 = 1; end
    end
    if im1 < 1 || ip1 > n
        continue;
    end

    A = P(im1,:); B = P(i,:); C = P(ip1,:);
    v1 = A - B; v2 = C - B;
    L1 = norm(v1); L2 = norm(v2);
    if L1 < 1e-12 || L2 < 1e-12
        continue;
    end
    u1 = v1 / L1; u2 = v2 / L2;

    cang = dot(u1, u2);
    cang = min(max(cang, -1), 1);
    phi = acos(cang);
    if phi < 1e-5 || abs(pi - phi) < 1e-5
        continue;
    end

    dtrim = r_in / tan(phi/2);
    dtrim = min([dtrim, 0.45*L1, 0.45*L2]);
    if ~isfinite(dtrim) || dtrim <= 1e-12
        continue;
    end
    re = dtrim * tan(phi/2);

    T1 = B + u1 * dtrim;
    T2 = B + u2 * dtrim;
    bvec = u1 + u2;
    nb = norm(bvec);
    if nb < 1e-12
        continue;
    end
    bvec = bvec / nb;
    O = B + bvec * (re / sin(phi/2));

    w1 = T1 - O;
    w2 = T2 - O;
    zc = w1(1)*w2(2) - w1(2)*w2(1);

    hasF(i) = true;
    Tin(i,:) = T1;
    Tout(i,:) = T2;
    rEff(i) = re;
    sweep(i) = zc > 0;
end

if is_closed
    startPt = P(1,:);
    if hasF(1), startPt = Tout(1,:); end
    d = sprintf('M %.9g %.9g', startPt(1), startPt(2));
    curPt = startPt;
    for i = 1:n
        j = i + 1;
        if j > n, j = 1; end
        endPt = P(j,:);
        if hasF(j), endPt = Tin(j,:); end
        if norm(endPt - curPt) > 1e-12
            d = sprintf('%s L %.9g %.9g', d, endPt(1), endPt(2)); %#ok<AGROW>
        end
        if hasF(j)
            d = sprintf('%s A %.9g %.9g 0 0 %d %.9g %.9g', ...
                d, rEff(j), rEff(j), sweep(j), Tout(j,1), Tout(j,2)); %#ok<AGROW>
            curPt = Tout(j,:);
        else
            curPt = endPt;
        end
    end
    d = sprintf('%s Z', d);
else
    d = sprintf('M %.9g %.9g', P(1,1), P(1,2));
    for i = 2:n-1
        if hasF(i)
            d = sprintf('%s L %.9g %.9g', d, Tin(i,1), Tin(i,2)); %#ok<AGROW>
            d = sprintf('%s A %.9g %.9g 0 0 %d %.9g %.9g', ...
                d, rEff(i), rEff(i), sweep(i), Tout(i,1), Tout(i,2)); %#ok<AGROW>
        else
            d = sprintf('%s L %.9g %.9g', d, P(i,1), P(i,2)); %#ok<AGROW>
        end
    end
    d = sprintf('%s L %.9g %.9g', d, P(end,1), P(end,2));
end
end

function [paths, is_closed, deg] = reconstruct_paths_from_edges(E)
% Reconstruct open/closed connected paths from undirected edges.
paths = {};
is_closed = false(0,1);
if isempty(E)
    deg = zeros(0,1);
    return;
end
E = double(E);
nE = size(E,1);
maxV = max(E(:));
adj = cell(maxV,1);
deg = zeros(maxV,1);
for e = 1:nE
    a = E(e,1); b = E(e,2);
    adj{a}(end+1) = e;
    adj{b}(end+1) = e;
    deg(a) = deg(a) + 1;
    deg(b) = deg(b) + 1;
end
used = false(nE,1);

starts = find(deg > 0 & deg ~= 2);
for s = starts(:).'
    while true
        e0 = first_unused_incident(s, adj, used);
        if e0 == 0, break; end
        [ids, used] = trace_path_from_start(s, e0, E, adj, used);
        if numel(ids) >= 2
            paths{end+1,1} = uint32(ids); %#ok<AGROW>
            is_closed(end+1,1) = false; %#ok<AGROW>
        end
    end
end

for e0 = 1:nE
    if used(e0), continue; end
    a = E(e0,1);
    [ids, used] = trace_path_from_start(a, e0, E, adj, used);
    if numel(ids) >= 3
        if ids(end) == ids(1)
            ids(end) = [];
        end
        paths{end+1,1} = uint32(ids); %#ok<AGROW>
        is_closed(end+1,1) = true; %#ok<AGROW>
    end
end
end

function tf = is_point_on_any_segment(P, A, B, tol)
tf = false;
if isempty(A) || isempty(B)
    return;
end
for i = 1:size(A,1)
    if point_on_segment(P, A(i,:), B(i,:), tol)
        tf = true;
        return;
    end
end
end

function tf = point_on_segment(P, A, B, tol)
AB = B - A;
AP = P - A;
LAB2 = dot(AB, AB);
if LAB2 < 1e-16
    tf = norm(P - A) <= tol;
    return;
end
t = dot(AP, AB) / LAB2;
if t < -1e-9 || t > 1+1e-9
    tf = false;
    return;
end
t = min(max(t, 0), 1);
Q = A + t*AB;
tf = norm(P - Q) <= tol;
end

function e0 = first_unused_incident(v, adj, used)
e0 = 0;
if v < 1 || v > numel(adj), return; end
cand = adj{v};
for k = 1:numel(cand)
    if ~used(cand(k))
        e0 = cand(k);
        return;
    end
end
end

function [ids, used] = trace_path_from_start(v_start, e_start, E, adj, used)
ids = v_start;
used(e_start) = true;
a = E(e_start,1); b = E(e_start,2);
if a == v_start
    curr = b;
else
    curr = a;
end
prev = v_start;
ids(end+1,1) = curr; %#ok<AGROW>
while true
    cand = adj{curr};
    nextEdge = 0; nextV = 0;
    for k = 1:numel(cand)
        ee = cand(k);
        if used(ee), continue; end
        u = E(ee,1); v = E(ee,2);
        nv = u + v - curr;
        if nv == prev
            continue;
        end
        nextEdge = ee;
        nextV = nv;
        break;
    end
    if nextEdge == 0
        break;
    end
    used(nextEdge) = true;
    prev = curr;
    curr = nextV;
    ids(end+1,1) = curr; %#ok<AGROW>
    if curr == v_start
        break;
    end
end
end

function paths = edges_to_path_segments(E, pts)
paths = cell(size(E,1),1);
for i = 1:size(E,1)
    a = double(E(i,1));
    b = double(E(i,2));
    pa = pts(a,:);
    pb = pts(b,:);
    paths{i} = sprintf('M %.9g %.9g L %.9g %.9g', pa(1), pa(2), pb(1), pb(2));
end
end

function key = snap_key(p, tol)
sx = round(p(1) / tol);
sy = round(p(2) / tol);
key = sprintf('%d_%d', sx, sy);
end

function key = edge_key(i, j)
key = sprintf('%u_%u', i, j);
end

function M = build_edge_key_set(E)
M = containers.Map('KeyType','char', 'ValueType','logical');
for i = 1:size(E,1)
    a = E(i,1);
    b = E(i,2);
    if a < b
        key = edge_key(a, b);
    else
        key = edge_key(b, a);
    end
    M(key) = true;
end
end

function [filename, add_engrave, snap_tol, fillet_radius] = resolve_inputs(varargin)
filename = 'tessellation_cut_pattern.svg';
add_engrave = false;
snap_tol = [];
fillet_radius = 0;

if nargin >= 1 && ~isempty(varargin{1})
    arg = varargin{1};
    if isstring(arg) && isscalar(arg)
        filename = char(arg);
    elseif ischar(arg)
        filename = arg;
    end
end

if nargin >= 2 && ~isempty(varargin{2})
    arg = varargin{2};
    if islogical(arg) && isscalar(arg)
        add_engrave = arg;
    elseif isnumeric(arg) && isscalar(arg)
        % Backward-compatible fallback: treat 3rd argument as snap_tol.
        snap_tol = double(arg);
    end
end

if nargin >= 3 && ~isempty(varargin{3})
    arg = varargin{3};
    if isnumeric(arg) && isscalar(arg) && isfinite(arg)
        snap_tol = double(arg);
    end
end
if nargin >= 4 && ~isempty(varargin{4})
    arg = varargin{4};
    if isnumeric(arg) && isscalar(arg) && isfinite(arg) && arg >= 0
        fillet_radius = double(arg);
    end
end
[~, name, ext] = fileparts(filename);
if isempty(name)
    filename = 'tessellation_cut_pattern.svg';
elseif isempty(ext)
    filename = [filename '.svg'];
end
end
