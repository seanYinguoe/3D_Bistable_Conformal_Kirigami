function [v_target, v_out] = grid_deployment(v_target, v_out, f_out)
%GRID_DEPLOYMENT Plot deployment between flattened and deployed grids.
%   [v_target, v_out] = grid_deployment(v_target, v_out, f_out)

if size(v_out,2) == 2
    v_out = [v_out, zeros(size(v_out,1),1)];
elseif size(v_out,2) ~= 3
    error('v_out must be Nx2 or Nx3.');
end
if size(v_target,2) ~= 3
    error('v_target must be Nx3.');
end
if size(f_out,2) ~= 3
    error('f_out must be Mx3.');
end

allV = [v_out; v_target];
vMin = min(allV, [], 1);
vMax = max(allV, [], 1);
span = max(vMax - vMin, 1e-9);
pad = 0.08 * max(span);
xLimFix = [vMin(1)-pad, vMax(1)+pad];
yLimFix = [vMin(2)-pad, vMax(2)+pad];
zLimFix = [min(0, vMin(3)-pad), vMax(3)+pad];

fig = figure('Name', 'Deployment Control', 'Color', 'w', ...
             'Units', 'pixels', 'Position', [120 80 980 720]);
slider = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [190 18 600 22], ...
                   'Min', 0, 'Max', 1, 'Value', 0, ...
                   'Callback', @(src,event) updatePlot(src, v_out, v_target, f_out, xLimFix, yLimFix, zLimFix));
updatePlot(slider, v_out, v_target, f_out, xLimFix, yLimFix, zLimFix);

end

function updatePlot(src, v_out, v_target, faces, xLimFix, yLimFix, zLimFix)
alpha = get(src, 'Value');
x_deploy = (1-alpha) * v_out(:,1) + alpha * v_target(:,1);
y_deploy = (1-alpha) * v_out(:,2) + alpha * v_target(:,2);
z_deploy = alpha * v_target(:,3);
v_deploy = [x_deploy, y_deploy, z_deploy];

cla();
patch('Vertices', v_deploy, 'Faces', faces, ...
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
