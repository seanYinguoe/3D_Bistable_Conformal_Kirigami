function output_path = generate_gif(tessellation, tessellation_target, varargin)
%GENERATE_GIF Create and save deployment GIF next to this function.
%   output_path = generate_gif(tessellation, tessellation_target)
%   output_path = generate_gif(tessellation, tessellation_target, filename)
%
% The GIF is saved to an "output" folder in the same directory as this
% function. If additional legacy inputs are passed, the last char/string
% input is treated as the filename.

filename = resolve_filename(varargin{:});
func_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(func_dir, 'output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_path = fullfile(output_dir, filename);

% Number of frames for the GIF
numFrames = 50;

% Create a figure for plotting
fig = figure('Name', 'Deployment Control', 'Position', [100, 100, 800, 600]);

for k = 0:numFrames
    alpha = k / numFrames;
    
    % Compute interpolated deployment shape
    tessellation_deploy = cell(size(tessellation));
    for m = 1:size(tessellation,1)
        x_deploy = (1-alpha) * tessellation{m}(:,1) + alpha * tessellation_target{m}(:,1);
        y_deploy = (1-alpha) * tessellation{m}(:,2) + alpha * tessellation_target{m}(:,2);
        z_deploy = alpha * tessellation_target{m}(:,3);
        tessellation_deploy{m} = [x_deploy, y_deploy, z_deploy];
    end

    % Clear previous plot and create new one
    clf(fig);
    set(fig, 'Color', 'w')
    %colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255};
    colour = {'white', [0.9216    0.8863    0.4235
], [ 0.7059    0.9608    0.4118], [0.9216    0.8863    0.4235]};
    hold on;
    for m = 1:size(tessellation,1)
        plot_triangle(tessellation_deploy{m},colour);
    end
    grid off;
    axis equal;
    axis off;
    set(gca, 'Color', 'w');      % Set axes background to white
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['Deployment Progress: ' num2str(alpha*100, '%.1f') '%']);
    view([45, 45]);
    hold off;

    % Capture the plot as an image
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Write to the GIF file
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
