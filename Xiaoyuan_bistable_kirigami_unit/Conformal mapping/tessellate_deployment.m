function tessellation_target = tessellate_deployment(tessellation,v_out,v_target,f_out)
% Deploy a flattened tessellation onto the target tessellation
% using barycentric interpolation.
%% expanse tessellation matrix to Nx3 from Nx2
for i = 1:size(tessellation,1)
    tessellation{i} = [tessellation{i} zeros(size(tessellation{i},1),1)];
end

%% Convert coordinate system from global cartesian cooridinate to local barycentric coordinate
bc_out = cell(size(tessellation));
for i = 1:size(tessellation,1)
    for j = 1:size(tessellation{i},1)
        tri = v_out(f_out(i,:),:);
        p = tessellation{i}(j,:);
        bc_out{i}(j,:) = cart2barycentric(tri,p);
        k = k+1;
    end
end

% Convert 2D tessellation surface to 3D tessellation surface based on barycentric
% coordinate and cooresponding 3D grid cartesian coordiante
tessellation_target = cell(size(tessellation)); % cartesian coordinates of deployed grid surface
for i = 1:size(tessellation,1)
    for j = 1:size(f_out,2)
        cart_grid = bc_out{i} * v_target(f_out(i,:),:);
    end
    tessellation_target{i} = cart_grid;
end


%% Deploy the flattened grids surface on to deployed target surface
% Create the figure and UI components
fig = figure('Name', 'Deployment Control', 'Position', [100 100 800 600]);
slider = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [150 20 500 20], ...
                   'Min', 0, 'Max', 1, 'Value', 0, ...
                   'Callback', @(src,event) updatePlot(src, tessellation, tessellation_target, f_out));

% Function to update the plot based on slider value
function updatePlot(src, v_out, v_target, faces)

alpha = get(src, 'Value');

% Compute interpolated deployment shape
v_deploy = cell(size(tessllation));
for m = 1:size(tessellation,1)
    x_deploy = (1-alpha) * v_out{m}(:,1) + alpha * v_target{m}(:,1);
    y_deploy = (1-alpha) * v_out{m}(:,2) + alpha * v_target{m}(:,2);
    z_deploy = alpha * v_target{m}(:,3);

    v_deploy{m} = [x_deploy, y_deploy, z_deploy];
end


% Clear previous plot and create new one
cla();
patch('Vertices', v_deploy, 'Faces', faces, ...
    'FaceColor', 'none', 'EdgeColor', 'black', 'FaceAlpha', 0.6);
%plot_triangle(tessellation{i},colour);
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