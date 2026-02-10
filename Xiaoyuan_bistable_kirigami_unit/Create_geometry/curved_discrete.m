function curved_discrete(path, filename)
% PLOT_CURVED_DISCRETE  Plot a discrete curved OBJ mesh in the same style as curved_surface
%
% Usage:
%   plot_curved_discrete('/path/to/folder/', 'curved_discrete.obj')
%
% Requirements:
%   - You already have readObj(path,filename) on your MATLAB path.
%
% Style:
%   - Lavender translucent surface
%   - Black mesh edges (decimated)
%   - Thick boundary outline
%   - Orthographic view + soft lighting

close all; clc;

%% ===================== 1) Read OBJ =====================
path = '/Users/sean/Desktop/Project 2/3D_Bistable_Conformal_Kirigami/Xiaoyuan_bistable_kirigami_unit/Create_geometry/';
filename = 'curved_discrete.obj';  % Specify the filename for the OBJ file;
obj = readObj(path, filename);

%% ===================== 2) Extract vertices + faces (robust) =====================
% Prefer obj.v for 3D geometry
if isfield(obj,'v') && ~isempty(obj.v)
    V = double(obj.v);
elseif isfield(obj,'vt') && ~isempty(obj.vt)
    % fallback (some OBJs store geometry in vt)
    V = double(obj.vt);
    if size(V,2)==2, V(:,3)=0; end
else
    error('readObj output has no vertices: expected obj.v or obj.vt');
end

% Faces: prefer obj.f.v
if ~isfield(obj,'f') || isempty(obj.f)
    error('readObj output has no face struct obj.f');
end

if isfield(obj.f,'v') && ~isempty(obj.f.v)
    F = double(obj.f.v);
elseif isfield(obj.f,'vt') && ~isempty(obj.f.vt)
    F = double(obj.f.vt);
else
    error('readObj output has no faces: expected obj.f.v or obj.f.vt');
end

% Keep triangles only (fan triangulate if needed)
if size(F,2) > 3
    F = fan_triangulate(F);
end

% Remove invalid faces
F = F(all(isfinite(F),2),:);
F = F(all(F>=1,2),:);
F = F(all(F<=size(V,1),2),:);

% If vertices have >3 cols, keep XYZ
if size(V,2) > 3
    V = V(:,1:3);
end

%% ===================== 3) Plot surface =====================
fig = figure('Color','w','Position',[200 200 620 460]);
ax  = axes(fig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');

lav = [0.82 0.83 0.92];

hSurf = trisurf(F, V(:,1), V(:,2), V(:,3), ...
    'EdgeColor','none', ...
    'FaceColor', lav, ...
    'FaceAlpha', 0.70);

%% ===================== 4) Draw mesh edges (grid) =====================
% Build unique edges
E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
E = sort(E,2);
E = unique(E,'rows');

% Decimate edges so figure stays clean
targetEdges = 2200;  % tune: 1500–4000 depending on mesh density
step = max(1, round(size(E,1)/targetEdges));
Eplot = E(1:step:end,:);

gridCol = [0 0 0 0.85];
gridLW  = 0.8;

% Plot edges as line segments
% (vectorized plotting would be faster, but this is robust & clear)
for i = 1:size(Eplot,1)
    i1 = Eplot(i,1); i2 = Eplot(i,2);
    plot3(ax, [V(i1,1) V(i2,1)], ...
              [V(i1,2) V(i2,2)], ...
              [V(i1,3) V(i2,3)], ...
          'Color', gridCol, 'LineWidth', gridLW);
end

%% ===================== 5) Boundary outline (thicker) =====================
TR = triangulation(F, V(:,1:3));
B  = freeBoundary(TR);   % boundary edges (unordered)

% Trace boundary loops (handles multiple components)
loops = boundary_loops_from_edges(B);

outlineLW = 2.0;
for k = 1:numel(loops)
    idx = loops{k};
    plot3(ax, V(idx,1), V(idx,2), V(idx,3), 'k-', 'LineWidth', outlineLW);
end

%% ===================== 6) Camera + lighting (paper style) =====================
% Similar to your curved_surface view
view(ax, 315, 18);
camproj(ax,'orthographic');

material(ax,'dull');
lighting(ax,'gouraud');

camlight(ax,'headlight');
L = light(ax,'Style','infinite');
L.Position = [-1 1 1];

hSurf.AmbientStrength  = 0.55;
hSurf.DiffuseStrength  = 0.60;
hSurf.SpecularStrength = 0.06;
hSurf.SpecularExponent = 25;

% Tight bounds
pad = 0.02;
xlim(ax, padlim(V(:,1), pad));
ylim(ax, padlim(V(:,2), pad));
zlim(ax, padlim(V(:,3), pad));

end

%% ===================== helpers =====================

function Ft = fan_triangulate(F)
% Convert polygon faces to triangles via fan (row-wise).
Ft = [];
for i = 1:size(F,1)
    row = F(i,:);
    row = row(isfinite(row) & row>0);
    if numel(row) < 3, continue; end
    if numel(row) == 3
        Ft = [Ft; row]; %#ok<AGROW>
    else
        v1 = row(1);
        for k = 2:(numel(row)-1)
            Ft = [Ft; v1 row(k) row(k+1)]; %#ok<AGROW>
        end
    end
end
end

function loops = boundary_loops_from_edges(B)
% B: nb x 2 boundary edges (unordered). Returns cell array of vertex loops.

adj = containers.Map('KeyType','int32','ValueType','any');
for i = 1:size(B,1)
    a = int32(B(i,1)); b = int32(B(i,2));
    if ~isKey(adj,a), adj(a) = int32([]); end
    if ~isKey(adj,b), adj(b) = int32([]); end
    adj(a) = [adj(a) b];
    adj(b) = [adj(b) a];
end

unused = true(size(B,1),1);
loops = {};

while any(unused)
    eidx = find(unused,1);
    a0 = B(eidx,1); b0 = B(eidx,2);
    unused(eidx) = false;

    loop = [a0 b0];
    prev = a0; curr = b0;

    while true
        neigh = double(adj(int32(curr)));
        if isempty(neigh), break; end

        if numel(neigh) == 1
            nxt = neigh(1);
        else
            if neigh(1) ~= prev
                nxt = neigh(1);
            else
                nxt = neigh(2);
            end
        end

        loop(end+1) = nxt; %#ok<AGROW>
        if nxt == loop(1), break; end

        % mark edge used
        mask = unused & ((B(:,1)==curr & B(:,2)==nxt) | (B(:,1)==nxt & B(:,2)==curr));
        if any(mask)
            unused(find(mask,1)) = false;
        end

        prev = curr; curr = nxt;

        if numel(loop) > 200000
            error('Boundary tracing got stuck (unexpected boundary topology).');
        end
    end

    loops{end+1} = loop; %#ok<AGROW>
end
end

function lim = padlim(x, padFrac)
xmin = min(x); xmax = max(x);
d = xmax - xmin;
if d < eps, d = 1; end
lim = [xmin - padFrac*d, xmax + padFrac*d];
end