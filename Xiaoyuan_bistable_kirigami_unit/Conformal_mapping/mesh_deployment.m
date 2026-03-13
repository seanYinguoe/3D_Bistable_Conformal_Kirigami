% Deploy a flattened surface onto the deployed surface (conformal mapping)
% using interpolation and analysis.
% Input: obj_2D(flattened_surface and deployed_surface)
function [flattened_surface,deployed_surface] = mesh_deployment(obj_2D)
% Input flattened and deployed surface
flattened_surface = obj_2D.vt;  % Flattened 2D vertex positions
flattened_surface = [flattened_surface,zeros(size(flattened_surface,1),1)];
deployed_surface = obj_2D.v;    % Deployed 3D vertex positions

% Keep flattened platform fixed, and rigidly align deployed mesh to it.
% This avoids model-dependent global rotation drift during interpolation.
[R, t] = kabsch_rigid(deployed_surface, flattened_surface);
deployed_surface = deployed_surface * R + t;

% Keep deployed surface above the flattened plane for clearer deployment.
if mean(deployed_surface(:,3)) < 0
    deployed_surface(:,3) = -deployed_surface(:,3);
end

% Recenter platform: move flattened centroid to (0,0) in XY,
% and apply the same XY shift to deployed surface.
c_flat_xy = mean(flattened_surface(:,1:2), 1);
flattened_surface(:,1:2) = flattened_surface(:,1:2) - c_flat_xy;
deployed_surface(:,1:2) = deployed_surface(:,1:2) - c_flat_xy;

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

function [R, t] = kabsch_rigid(A, B)
% Solve rigid map A*R + t ~= B (no scaling, no reflection).
if size(A,2) ~= 3 || size(B,2) ~= 3 || size(A,1) ~= size(B,1)
    error('kabsch_rigid expects A,B as Nx3 with matching N.');
end

cA = mean(A, 1);
cB = mean(B, 1);
A0 = A - cA;
B0 = B - cB;

H = A0' * B0;
[U, ~, V] = svd(H);
R = U * V';
if det(R) < 0
    U(:,end) = -U(:,end);
    R = U * V';
end
t = cB - cA * R;
end
