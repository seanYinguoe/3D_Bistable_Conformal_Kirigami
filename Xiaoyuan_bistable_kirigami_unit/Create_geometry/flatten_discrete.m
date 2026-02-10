function flatten_discrete(path, filename)
% FLATTEN_DISCRETE  Plot a flattened OBJ mesh (2D) using vt + f.vt
% Style consistent with curved_surface:
%   - Lavender translucent fill
%   - Black mesh edges (decimated)
%   - Thick boundary outline
%   - Orthographic view

close all; clc;

%% ===================== 1) Read OBJ =====================
% Use input arguments if provided, otherwise fall back to your defaults
if nargin < 1 || isempty(path)
    path = '/Users/sean/Desktop/Project 2/3D_Bistable_Conformal_Kirigami/Xiaoyuan_bistable_kirigami_unit/Create_geometry/';
end
if nargin < 2 || isempty(filename)
    filename = 'flatten_discrete.obj';
end

obj = readObj(path, filename);

%% ===================== 2) Extract vt vertices + vt faces =====================
if ~isfield(obj,'vt') || isempty(obj.vt)
    error('Expected flattened vertices in obj.vt, but obj.vt is missing/empty.');
end
V = double(obj.vt);

if ~isfield(obj,'f') || isempty(obj.f) || ~isfield(obj.f,'vt') || isempty(obj.f.vt)
    error('Expected faces indexing vt in obj.f.vt, but obj.f.vt is missing/empty.');
end
F = double(obj.f.vt);

% Triangulate if needed
if size(F,2) > 3
    F = fan_triangulate(F);
end

% V is Nx2 -> make Nx3 for trisurf/plot3 convenience
if size(V,2) == 2
    V = [V, zeros(size(V,1),1)];
else
    V(:,3) = 0; % force planar display
end

% Clean faces
F = F(all(isfinite(F),2),:);
F = F(all(F>=1,2),:);
F = F(all(F<=size(V,1),2),:);

%% ===================== 3) Plot surface =====================
fig = figure('Color','w','Position',[200 200 620 460]);
ax  = axes(fig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');

lav = [0.82 0.83 0.92];

hSurf = trisurf(F, V(:,1), V(:,2), V(:,3), ...
    'EdgeColor','none', ...
    'FaceColor', lav, ...
    'FaceAlpha', 0.70);

%% ===================== 4) Draw mesh edges (grid) =====================
E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
E = sort(E,2);
E = unique(E,'rows');

targetEdges = 2200;               % tune as needed
step = max(1, round(size(E,1)/targetEdges));
Eplot = E(1:step:end,:);

gridCol = [0 0 0 0.85];
gridLW  = 0.8;

for i = 1:size(Eplot,1)
    i1 = Eplot(i,1); i2 = Eplot(i,2);
    plot(ax, [V(i1,1) V(i2,1)], [V(i1,2) V(i2,2)], ...
        'Color', gridCol, 'LineWidth', gridLW);
end

%% ===================== 5) Boundary outline (thicker) =====================
TR = triangulation(F, V(:,1:2));
B  = freeBoundary(TR);   % boundary edges (unordered)
loops = boundary_loops_from_edges(B);

outlineLW = 2.0;
for k = 1:numel(loops)
    idx = loops{k};
    plot(ax, V(idx,1), V(idx,2), 'k-', 'LineWidth', outlineLW);
end

%% ===================== 6) Camera / projection =====================
view(ax, 2);
camproj(ax,'orthographic');

% Lighting is optional on a flat surface; keep minimal to match style
material(ax,'dull');
lighting(ax,'gouraud');
camlight(ax,'headlight');

hSurf.AmbientStrength  = 0.55;
hSurf.DiffuseStrength  = 0.60;
hSurf.SpecularStrength = 0.06;
hSurf.SpecularExponent = 25;

pad = 0.02;
xlim(ax, padlim(V(:,1), pad));
ylim(ax, padlim(V(:,2), pad));

end

%% ===================== helpers =====================
function Ft = fan_triangulate(F)
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