function tessellation_target = tessellation_deployment(tessellation,v_out,v_target,f_out,T,stretch_facs,params)
% Deploy a flattened tessellation onto the target tessellation
% using barycentric interpolation.
%% Modify 2D closed tessellation configuration
for i = 1:size(tessellation,1)
    tessellation{i} = [tessellation{i} zeros(size(tessellation{i},1),1)]; % expanse tessellation matrix to Nx3 from Nx2
    tessellation{i} = tessellation{i} - [T(:,[1,2]),0];% Move it to the central point   
end

v_out = v_out - [T(:,[1,2]),0]; % Move the grids to match the centroid of flattened surface 

%% Get 2D clsosed tessellation configuration based on stretch_facs and rescale to the same edgeLen
% Get the parameters
edgeLen = params(1);
l1 = params(2);
l2 = params(3);
l3 = params(4);
t = params(5);
delta = (stretch_facs - 1) * edgeLen; % should be within the delta range
triangle = triangle_unit(pi/3, 2*pi/3, delta,l1,l2,l3,t);
triangle = triangle / stretch_facs; % rescale triangle to the same length
plot_triangle(triangle)


%% Get 3D open tessellation configuration
% Convert coordinate system from global cartesian cooridinate to local barycentric coordinate
bc_out = cell(size(tessellation));
for i = 1:size(tessellation,1)
    for j = 1:size(tessellation{i},1)
        tri = v_out(f_out(i,:),:);
        p = tessellation{i}(j,:);
        bc_out{i}(j,:) = cart2barycentric(tri,p);
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
                   'Callback', @(src,event) updatePlot(src, tessellation, tessellation_target));

% Function to update the plot based on slider value
function updatePlot(src, tessellation, tessellation_target)

alpha = get(src, 'Value');

% Compute interpolated deployment shape
tessellation_deploy = cell(size(tessellation));
for m = 1:size(tessellation,1)
    x_deploy = (1-alpha) * tessellation{m}(:,1) + alpha * tessellation_target{m}(:,1);
    y_deploy = (1-alpha) * tessellation{m}(:,2) + alpha * tessellation_target{m}(:,2);
    z_deploy = alpha * tessellation_target{m}(:,3);
    tessellation_deploy{m} = [x_deploy, y_deploy, z_deploy];
end


% Clear previous plot and create new one
cla();
colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255}; % The colour of void, flank, filament, Innertriangle
for m = 1:size(tessellation,1)
    plot_triangle(tessellation_deploy{m},colour);
end
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