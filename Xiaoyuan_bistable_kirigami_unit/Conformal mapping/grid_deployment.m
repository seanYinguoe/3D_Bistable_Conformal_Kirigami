function grid_deployment(obj_2D,v_out,f_out)
% Deploy a flattened grids onto the deployed surface
% using interpolation and analysis.
% Input: flattened_surface (vertices, faces)
%        deployed_surface (vertices, faces)

%% Input flattened mesh surface and deployed mesh surface
% Input flattened and deployed surface
flattened_surface = obj_2D.vt;  % Flattened 2D vertex positions
flattened_surface = [flattened_surface,zeros(size(flattened_surface,1),1)];
deployed_surface = obj_2D.v;    % Deployed 3D vertex positions

% Move and rotate the 3D object to match the boundary
deployed_surface = [deployed_surface(:,1), -deployed_surface(:,3),deployed_surface(:,2)]; % rotate the deployed surface 90 around x axis
T = mean(flattened_surface) - mean(deployed_surface); 
flattened_surface = flattened_surface - [T(:,[1,2]),0]; % Move to match the coorespoind node

%% Build the cooresponding between grids point and mesh points
% Find the closest mesh point to each grid point
v_out = v_out - [T(:,[1,2]),0]; % Move the grids to match the centroid of flattened surface
dis = pdist2(v_out, flattened_surface, 'euclidean');
[~, closest_indices] = min(dis, [], 2);

%% Deploy the flattened grids surface on to deployed target surface
% Create the figure and UI components
fig = figure('Name', 'Deployment Control', 'Position', [100 100 800 600]);
slider = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [150 20 500 20], ...
                   'Min', 0, 'Max', 1, 'Value', 0, ...
                   'Callback', @(src,event) updatePlot(src, flattened_surface, deployed_surface,v_out,f_out));

% Function to update the plot based on slider value
function updatePlot(src, flattened_surface, deployed_surface,v_out, faces)

alpha = get(src, 'Value');

% Compute interpolated deployment shape
x_deploy = (1-alpha) * flattened_surface(:,1) + alpha * deployed_surface(:,1);
y_deploy = (1-alpha) * flattened_surface(:,2) + alpha * deployed_surface(:,2);
z_deploy = alpha * deployed_surface(:,3);  % lift into 3D

v_deploy = [x_deploy, y_deploy, z_deploy];

% Define the grid points to interpolate 
delta_v = v_deploy - flattened_surface;

% Interpolate each displacement component using griddata instead of scatteredInterpolant
delta_v_out_x = griddata(v_deploy(:,1), v_deploy(:,2), delta_v(:,1), v_out(:,1), v_out(:,2), 'nearest');
delta_v_out_y = griddata(v_deploy(:,1), v_deploy(:,2), delta_v(:,2), v_out(:,1), v_out(:,2), 'nearest');
delta_v_out_z = griddata(v_deploy(:,1), v_deploy(:,2), delta_v(:,3), v_out(:,1), v_out(:,2), 'nearest');

% Compute the final deployed grid points in 3D
deployed_v_out = v_out + [delta_v_out_x, delta_v_out_y, delta_v_out_z];

% Clear previous plot and create new one
cla();
patch('Vertices', deployed_v_out, 'Faces', faces, ...
    'FaceColor', 'none', 'EdgeColor', 'black', 'FaceAlpha', 0.6);

grid off;
axis equal;
xlim([-90,90]);
ylim([-90,90]);
zlim([0,60]);
xlabel('X');
ylabel('Y');
zlabel('Z');
title(['Deployment Progress: ' num2str(alpha*100, '%.1f') '%']);
view([45, 45]);
end
end