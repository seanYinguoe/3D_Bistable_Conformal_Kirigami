% Deploy a flattened surface onto the deployed surface (conformal mapping)
% using interpolation and analysis.
% Input: obj_2D(flattened_surface and deployed_surface)
function [flattened_surface,deployed_surface] = mesh_deployment(obj_2D)
% Input flattened and deployed surface
flattened_surface = obj_2D.vt;  % Flattened 2D vertex positions
flattened_surface = [flattened_surface,zeros(size(flattened_surface,1),1)];
deployed_surface = obj_2D.v;    % Deployed 3D vertex positions

% Automatic rigid transform: remove global translation/rotation and align to XY plane
[flattened_surface, deployed_surface] = model_transform(flattened_surface, deployed_surface);
faces = obj_2D.f.v;             % Face connectivity

% Precompute fixed view bounds so window/axes do not jump during updates
allV = [flattened_surface; deployed_surface];
vMin = min(allV, [], 1);
vMax = max(allV, [], 1);
span = max(vMax - vMin, 1e-9);
pad = 0.08 * max(span); % isotropic padding for stable framing
xLimFix = [vMin(1)-pad, vMax(1)+pad];
yLimFix = [vMin(2)-pad, vMax(2)+pad];
zLimFix = [min(0, vMin(3)-pad), vMax(3)+pad];

% Create the figure and UI components (fixed, readable size)
fig = figure('Name', 'Deployment Control', 'Color', 'w', ...
             'Units', 'pixels', 'Position', [120 80 980 720]);
slider = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [190 18 600 22], ...
                   'Min', 0, 'Max', 1, 'Value', 0, ...
                   'Callback', @(src,event) updatePlot(src, flattened_surface, deployed_surface, faces, xLimFix, yLimFix, zLimFix));
updatePlot(slider, flattened_surface, deployed_surface, faces, xLimFix, yLimFix, zLimFix);


% Function to update the plot based on slider value
function updatePlot(src, flattened_surface, deployed_surface, faces, xLimFix, yLimFix, zLimFix)

alpha = get(src, 'Value');

% Compute interpolated deployment shape
X_deploy = (1-alpha) * flattened_surface(:,1) + alpha * deployed_surface(:,1);
Y_deploy = (1-alpha) * flattened_surface(:,2) + alpha * deployed_surface(:,2);
Z_deploy = alpha * deployed_surface(:,3);  % Gradually lift into 3D

V_deploy = [X_deploy, Y_deploy, Z_deploy];

% Clear previous plot and create new one
cla();
patch('Vertices', V_deploy, 'Faces', faces, ...
    'FaceColor', 'none', 'EdgeColor', 'black', 'FaceAlpha', 0.6);

grid off;
axis equal;
xlim(xLimFix);
ylim(yLimFix);
zlim(zLimFix);
axis vis3d;
axis off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title(['Deployment Progress: ' num2str(alpha*100, '%.1f') '%']);
view([45, 45]);
end
end
