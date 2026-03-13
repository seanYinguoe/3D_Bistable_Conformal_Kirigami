function output_path = generate_gif(tessellation, tessellation_target, varargin)
%GENERATE_GIF Create and save deployment GIF next to this function.
%   output_path = generate_gif(tessellation, tessellation_target)
%   output_path = generate_gif(tessellation, tessellation_target, filename)
%   output_path = generate_gif(vt_mesh, v_mesh, filename)
%   output_path = generate_gif(vt_mesh, v_mesh, f_mesh, filename)
%
% The GIF is saved to an "output" folder in the same directory as this
% function. Last char/string varargin is treated as filename.

filename = resolve_filename(varargin{:});
[faces_mesh, is_mesh_mode] = resolve_faces_and_mode(tessellation, tessellation_target, varargin{:});

func_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(func_dir, 'output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_path = fullfile(output_dir, filename);

numFrames = 50;

% Keep figure size appropriate to the current screen.
screen_size = get(groot, 'ScreenSize'); % [left bottom width height]
fig_w = min(max(round(0.72 * screen_size(3)), 900), 1400);
fig_h = min(max(round(0.72 * screen_size(4)), 650), 1000);
fig_x = max(1, round((screen_size(3) - fig_w) / 2));
fig_y = max(1, round((screen_size(4) - fig_h) / 2));
fig = figure('Name', 'Deployment Control', ...
             'Units', 'pixels', ...
             'Position', [fig_x, fig_y, fig_w, fig_h], ...
             'Resize', 'off');

% Fixed axis bounds to avoid frame-to-frame rescaling/zooming.
if is_mesh_mode
    allV = [ensure_xyz(tessellation); ensure_xyz(tessellation_target)];
else
    allV = collect_tessellation_vertices(tessellation, tessellation_target);
end
vMin = min(allV, [], 1);
vMax = max(allV, [], 1);
span = max(vMax - vMin, 1e-9);
pad = 0.08 * max(span);
xLimFix = [vMin(1)-pad, vMax(1)+pad];
yLimFix = [vMin(2)-pad, vMax(2)+pad];
zLimFix = [min(0, vMin(3)-pad), vMax(3)+pad];

for k = 0:numFrames
    alpha = k / numFrames;

    clf(fig);
    set(fig, 'Color', 'w');
    hold on;

    if is_mesh_mode
        v_deploy = interpolate_mesh_vertices(tessellation, tessellation_target, alpha);
        patch('Vertices', v_deploy, 'Faces', faces_mesh, ...
              'FaceColor', 'none', 'EdgeColor', 'black', 'FaceAlpha', 0.6);
    else
        colour = {'white', [0.9216 0.8863 0.4235], [0.7059 0.9608 0.4118], [0.9216 0.8863 0.4235]};
        for m = 1:numel(tessellation)
            x_deploy = (1-alpha) * tessellation{m}(:,1) + alpha * tessellation_target{m}(:,1);
            y_deploy = (1-alpha) * tessellation{m}(:,2) + alpha * tessellation_target{m}(:,2);
            z_deploy = alpha * tessellation_target{m}(:,3);
            tessellation_deploy = [x_deploy, y_deploy, z_deploy];
            plot_triangle(tessellation_deploy, colour);
        end
    end

    grid off;
    axis equal;
    xlim(xLimFix);
    ylim(yLimFix);
    zlim(zLimFix);
    axis vis3d;
    axis off;
    set(gca, 'Color', 'w');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['Deployment Progress: ' num2str(alpha*100, '%.1f') '%']);
    view([45, 45]);
    %view([45, 25]);
    hold off;

    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    if k == 0
        imwrite(imind, cm, output_path, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, output_path, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

close(fig);
end

function filename = resolve_filename(varargin)
filename = 'deployment.gif';

if nargin == 0
    return;
end

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
    filename = 'deployment.gif';
elseif isempty(ext)
    filename = [filename '.gif'];
end
end

function [faces_mesh, is_mesh_mode] = resolve_faces_and_mode(tessellation, tessellation_target, varargin)
faces_mesh = [];
is_mesh_mode = isnumeric(tessellation) && isnumeric(tessellation_target);

if ~is_mesh_mode
    return;
end

if size(tessellation,2) == 2
    tessellation = [tessellation, zeros(size(tessellation,1),1)];
end
if size(tessellation_target,2) == 2
    tessellation_target = [tessellation_target, zeros(size(tessellation_target,1),1)];
end
if size(tessellation,2) ~= 3 || size(tessellation_target,2) ~= 3
    error('Mesh mode requires Nx3 vertices for start and target.');
end
if size(tessellation,1) ~= size(tessellation_target,1)
    error('Mesh mode requires same number of vertices in start/target.');
end

for k = 1:numel(varargin)
    arg = varargin{k};
    if isnumeric(arg) && size(arg,2) == 3 && size(arg,1) >= 1 && all(isfinite(arg(:)))
        faces_mesh = round(arg);
        break;
    end
end

if isempty(faces_mesh)
    faces_mesh = delaunay(tessellation(:,1), tessellation(:,2));
end
end

function v_deploy = interpolate_mesh_vertices(v0, v1, alpha)
v0 = ensure_xyz(v0);
v1 = ensure_xyz(v1);
x_deploy = (1-alpha) * v0(:,1) + alpha * v1(:,1);
y_deploy = (1-alpha) * v0(:,2) + alpha * v1(:,2);
z_deploy = alpha * v1(:,3);
v_deploy = [x_deploy, y_deploy, z_deploy];
end

function V = ensure_xyz(V)
if size(V,2) == 2
    V = [V, zeros(size(V,1),1)];
end
end

function allV = collect_tessellation_vertices(t0, t1)
allV = [];
for i = 1:numel(t0)
    a = ensure_xyz(t0{i});
    b = ensure_xyz(t1{i});
    allV = [allV; a; b]; %#ok<AGROW>
end
end
