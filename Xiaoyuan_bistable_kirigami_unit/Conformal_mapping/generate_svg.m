function [output_path, info] = generate_svg(tessellation, varargin)
%GENERATE_SVG Export laser-cut SVG from compact-state tessellation.
%   output_path = generate_svg(tessellation)
%   output_path = generate_svg(tessellation, filename)
%   [output_path, info] = generate_svg(...)
%
% Path order in SVG:
%   1) outer boundary cut(s)
%   2) internal void cut(s)

if nargin < 1 || ~iscell(tessellation)
    error('tessellation must be a cell array of unit vertex sets.');
end

filename = resolve_svg_filename(varargin{:});
func_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(func_dir, 'output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_path = fullfile(output_dir, filename);

% Face connectivity consistent with plot_triangle.m
f_void = [1 2 3 4 5 6; 7 8 9 10 11 12; 13 14 15 16 17 18];
f_flank = [19 20 21 22; 23 24 25 26; 27 28 29 30];
f_filament = [31 32 33 34; 35 36 37 38; 39 40 41 42];
f_inner = [43 44 45];

solid_all = polyshape();
void_all = polyshape();
has_solid = false;
has_void = false;
all_xy = zeros(0,2);

for k = 1:numel(tessellation)
    V = tessellation{k};
    if isempty(V)
        continue;
    end
    if size(V,2) < 2 || size(V,1) < 45
        error('tessellation{%d} must be at least 45x2.', k);
    end
    XY = V(:,1:2);
    all_xy = [all_xy; XY]; %#ok<AGROW>

    p_solid = polyshape();
    for j = 1:size(f_flank,1)
        p_solid = union(p_solid, polyshape(XY(f_flank(j,:),1), XY(f_flank(j,:),2), ...
            'Simplify', true, 'KeepCollinearPoints', true));
    end
    for j = 1:size(f_filament,1)
        p_solid = union(p_solid, polyshape(XY(f_filament(j,:),1), XY(f_filament(j,:),2), ...
            'Simplify', true, 'KeepCollinearPoints', true));
    end
    p_solid = union(p_solid, polyshape(XY(f_inner,1), XY(f_inner,2), ...
        'Simplify', true, 'KeepCollinearPoints', true));

    p_void = polyshape();
    for j = 1:size(f_void,1)
        p_void = union(p_void, polyshape(XY(f_void(j,:),1), XY(f_void(j,:),2), ...
            'Simplify', true, 'KeepCollinearPoints', true));
    end

    if ~has_solid
        solid_all = p_solid;
        has_solid = true;
    else
        solid_all = union(solid_all, p_solid);
    end
    if ~has_void
        void_all = p_void;
        has_void = true;
    else
        void_all = union(void_all, p_void);
    end
end

if isempty(all_xy)
    error('tessellation is empty.');
end

% Boundary cut(s): only outer contours of the solid sheet
outer_sheet = rmholes(solid_all);

% Void cut(s): internal void polygons clipped to sheet
void_cut = intersect(void_all, outer_sheet);
if (~has_void) || isempty(void_cut) || area(void_cut) <= 0
    % fallback: extract holes directly from solid geometry
    [xh, yh] = holes(solid_all);
    void_path_data = loops_to_paths(xh, yh);
else
    [xv, yv] = boundary(void_cut);
    void_path_data = loops_to_paths(xv, yv);
end

[xb, yb] = boundary(outer_sheet);
boundary_path_data = loops_to_paths(xb, yb);

if isempty(boundary_path_data)
    error('No boundary extracted from tessellation.');
end

% SVG canvas
xmin = min(all_xy(:,1)); xmax = max(all_xy(:,1));
ymin = min(all_xy(:,2)); ymax = max(all_xy(:,2));
span = max([xmax-xmin, ymax-ymin, 1e-6]);
pad = 0.02 * span;
vb_x = xmin - pad;
vb_y = ymin - pad;
vb_w = (xmax - xmin) + 2*pad;
vb_h = (ymax - ymin) + 2*pad;

fid = fopen(output_path, 'w');
if fid < 0
    error('Cannot open file for writing: %s', output_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid, '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" ');
fprintf(fid, 'viewBox="%.9g %.9g %.9g %.9g">\n', vb_x, vb_y, vb_w, vb_h);
fprintf(fid, '  <title>Kirigami cut pattern</title>\n');
fprintf(fid, '  <desc>Path order: boundary first, then void cuts.</desc>\n');

% Cut style: red stroke, no fill (common laser workflow)
fprintf(fid, '  <g id="boundary_cut" fill="none" stroke="#ff0000" stroke-width="0.05">\n');
for i = 1:numel(boundary_path_data)
    fprintf(fid, '    <path d="%s" />\n', boundary_path_data{i});
end
fprintf(fid, '  </g>\n');

fprintf(fid, '  <g id="void_cut" fill="none" stroke="#ff0000" stroke-width="0.05">\n');
for i = 1:numel(void_path_data)
    fprintf(fid, '    <path d="%s" />\n', void_path_data{i});
end
fprintf(fid, '  </g>\n');
fprintf(fid, '</svg>\n');

info = struct();
info.n_boundary_paths = numel(boundary_path_data);
info.n_void_paths = numel(void_path_data);
info.output_path = output_path;
info.viewBox = [vb_x, vb_y, vb_w, vb_h];
end

function paths = loops_to_paths(x, y)
paths = {};
if isempty(x) || isempty(y)
    return;
end
nanBreak = isnan(x) | isnan(y);
idx = [0; find(nanBreak); numel(x)+1];
for k = 1:numel(idx)-1
    s = idx(k)+1;
    e = idx(k+1)-1;
    if e - s + 1 < 3
        continue;
    end
    xx = x(s:e);
    yy = y(s:e);
    if hypot(xx(end)-xx(1), yy(end)-yy(1)) < 1e-12
        xx(end) = [];
        yy(end) = [];
    end
    if numel(xx) < 3
        continue;
    end
    d = sprintf('M %.9g %.9g', xx(1), yy(1));
    for i = 2:numel(xx)
        d = sprintf('%s L %.9g %.9g', d, xx(i), yy(i)); %#ok<AGROW>
    end
    d = sprintf('%s Z', d);
    paths{end+1} = d; %#ok<AGROW>
end
end

function filename = resolve_svg_filename(varargin)
filename = 'tessellation_pattern.svg';
for k = numel(varargin):-1:1
    arg = varargin{k};
    if isstring(arg) && isscalar(arg)
        filename = char(arg);
        break;
    end
    if ischar(arg)
        filename = arg;
        break;
    end
end
[~, name, ext] = fileparts(filename);
if isempty(name)
    filename = 'tessellation_pattern.svg';
elseif isempty(ext)
    filename = [filename '.svg'];
end
end
