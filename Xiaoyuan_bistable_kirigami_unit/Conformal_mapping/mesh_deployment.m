% Deploy a flattened surface onto the deployed surface (conformal mapping)
% using interpolation and analysis.
% Input: obj_2D(flattened_surface and deployed_surface)
function mesh_deployment(obj_2D,modelname)
% Input flattened and deployed surface
flattened_surface = obj_2D.vt;  % Flattened 2D vertex positions
flattened_surface = [flattened_surface,zeros(size(flattened_surface,1),1)];
deployed_surface = obj_2D.v;    % Deployed 3D vertex positions

% Transform 2D and 3D configutation to let them on the same platform
[flattened_surface, deployed_surface] = model_rotate(modelname, flattened_surface, deployed_surface);
%deployed_surface = [deployed_surface(:,1), -deployed_surface(:,3),deployed_surface(:,2)]; % rotate the deployed surface 90 around x axis
T = mean(flattened_surface) - mean(deployed_surface); 
flattened_surface = flattened_surface - [T(:,[1,2]),0]; % Move to match the coorespoind node
faces = obj_2D.f.v;             % Face connectivity

% Create the figure and UI components
fig = figure('Name', 'Deployment Control', 'Position', [100 100 800 600]);
slider = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [150 20 500 20], ...
                   'Min', 0, 'Max', 1, 'Value', 0, ...
                   'Callback', @(src,event) updatePlot(src, flattened_surface, deployed_surface,faces));
updatePlot(slider, flattened_surface, deployed_surface,faces);


% Function to update the plot based on slider value
function updatePlot(src, flattened_surface, deployed_surface,faces)

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
% xlim([-70,70]);
% ylim([-70,70]);
% zlim([0,60]);
axis off
xlabel('X');
ylabel('Y');
zlabel('Z');
title(['Deployment Progress: ' num2str(alpha*100, '%.1f') '%']);
view([45, 45]);
end
end