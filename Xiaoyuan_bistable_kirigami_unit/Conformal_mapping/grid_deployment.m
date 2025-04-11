function [T,v_target] = grid_deployment(obj_2D,c_mesh,v_out,f_out)
% Deploy a flattened grids onto the deployed surface
% using barycentric interpolation.

%% Input flattened mesh surface and deployed mesh surface
% Input flattened and deployed surface
flattened_surface = obj_2D.vt;  % Flattened 2D vertex positions
flattened_surface = [flattened_surface,zeros(size(flattened_surface,1),1)];
deployed_surface = obj_2D.v;    % Deployed 3D vertex positions
c_mesh = [c_mesh zeros(size(c_mesh,1),1)];
f_mesh = obj_2D.f.v;             % Face connectivity

% Move and rotate the 3D object to match the boundary
deployed_surface = [deployed_surface(:,1), -deployed_surface(:,3),deployed_surface(:,2)]; % rotate the deployed surface 90 around x axis
T = mean(flattened_surface) - mean(deployed_surface); 
flattened_surface = flattened_surface - [T(:,[1,2]),0]; % Move to match the coorespoind node
c_mesh = c_mesh - [T(:,[1,2]),0]; % Move to match the coorespoind node
v_out = v_out - [T(:,[1,2]),0]; % Move the grids to match the centroid of flattened surface

%% Build the cooresponding between grids point and mesh points for each grid triangle
cp = cell(size(f_out,1),1);
for i = 1:size(f_out,1)
    closest_indices = [];
    for j = 1:size(f_out,2)
        index = f_out(i,j);
        dis = pdist2(v_out(index,:), c_mesh, 'euclidean');
        [~, nearest] = min(dis, [], 2);
        closest_indices = [closest_indices,nearest];
    end
    cp{i} = closest_indices; % cell n represent, nth grid face(three grid points) and their nearest mesh centroid.
end

%% Convert coordinates
% Convert 2D grid surface from cart cooridinate system to barycentric coordinate system based on mesh surface
bc_out = cell(size(f_out,1)*size(f_out,2),1);
k = 1;
for i = 1:size(f_out,1)
    for j = 1:size(f_out,2)
        index = cp{i}(j);
        tri = flattened_surface(f_mesh(index,:),:);
        p = v_out(f_out(i,j),:);
        bc_out{k} = cart2barycentric(tri,p);
        k = k+1;
    end
end

% find the cooresponding 3D cartesian coordinate
cart_mesh = cell(size(f_out,1)*size(f_out,2),1);
k = 1;
for i = 1:size(cp,1)
    for j = 1:size(f_mesh,2)
        index = cp{i}(j);
        cart_mesh{k} = deployed_surface(f_mesh(index,:),:);
        k = k+1;
    end
end

% Convert 2D grid surface to 3D grid surface based on barycentric
% coordinate and cooresponding 3D cartesian coordiante
cart_out = []; % cartesian coordinates of deployed grid surface
for i = 1:size(bc_out,1)
    cart_grid = 0;
    for j = 1:size(f_out,2)
        cart_grid = cart_grid + cart_mesh{i}(j,:) * bc_out{i}(j);
    end
    cart_out(i,:) = cart_grid;
end

% reorder v_target to v_out order
v_target = zeros(size(v_out));
for i = 1:size(f_out,1)
    for j = 1:size(f_out,2)
        v_target(f_out(i,j),:) = cart_out(3*(i-1)+j,:);
    end
end




%% Deploy the flattened grids surface on to deployed target surface
% Create the figure and UI components
fig = figure('Name', 'Deployment Control', 'Position', [100 100 800 600]);
slider = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [150 20 500 20], ...
                   'Min', 0, 'Max', 1, 'Value', 0, ...
                   'Callback', @(src,event) updatePlot(src, v_out, v_target, f_out));

% Function to update the plot based on slider value
function updatePlot(src, v_out, v_target, faces)

alpha = get(src, 'Value');

% Compute interpolated deployment shape
x_deploy = (1-alpha) * v_out(:,1) + alpha * v_target(:,1);
y_deploy = (1-alpha) * v_out(:,2) + alpha * v_target(:,2);
z_deploy = alpha * v_target(:,3); 

v_deploy = [x_deploy, y_deploy, z_deploy];


% Clear previous plot and create new one
cla();
patch('Vertices', v_deploy, 'Faces', faces, ...
    'FaceColor', 'none', 'EdgeColor', 'black', 'FaceAlpha', 0.6);

grid off;
axis equal;
xlim([-90,90]);
ylim([-90,90]);
zlim([-10,60]);
xlabel('X');
ylabel('Y');
zlabel('Z');
title(['Deployment Progress: ' num2str(alpha*100, '%.1f') '%']);
view([45, 45]);
end
end