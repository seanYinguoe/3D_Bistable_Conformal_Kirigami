function tessellation_target = tessellation_deployment(tessellation,tessellation_close,v_out,v_target,f_out,T,def_facs)
% Deploy a flattened tessellation onto the target tessellation
% using barycentric interpolation.
%% Modify 2D closed tessellation configuration
for i = 1:size(tessellation,1)
    tessellation{i} = [tessellation{i} zeros(size(tessellation{i},1),1)]; % expanse tessellation matrix to Nx3 from Nx2
    tessellation{i} = tessellation{i} - [T(:,[1,2]),0];% Move it to the central point   
end

v_out = v_out - [T(:,[1,2]),0]; % Move the grids to match the centroid of flattened surface 

%% Get 2D clsosed tessellation configuration based on stretch_facs and rescale to the same edgeLen
for i = 1:size(tessellation_close,1)
    tessellation_close{i} = [tessellation_close{i} zeros(size(tessellation_close{i},1),1)]; % expanse tessellation matrix to Nx3 from Nx2
    tessellation_close{i} = tessellation_close{i} - [T(:,[1,2]),0];% Move it to the central point   
end

%% Get 3D open tessellation configuration
% Convert coordinate system from global cartesian cooridinate to local barycentric coordinate
bc_out = cell(size(tessellation_close));
for i = 1:size(tessellation_close,1)
    for j = 1:size(tessellation_close{i},1)
        tri = v_out(f_out(i,:),:);
        p = tessellation_close{i}(j,:);
        bc_out{i}(j,:) = cart2barycentric(tri,p);
        bc_out{i}(j,:) = [bc_out{i}(j,1)/def_facs(1),...
            bc_out{i}(j,2)/def_facs(2),...
            bc_out{i}(j,3)/def_facs(3)];% Use deforme factor to modify the carycentric coordinate
        bc_out{i}(j,:) = bc_out{i}(j,:)/sum(bc_out{i}(j,:));% enforce sum of bc_out{1} equal to 1
    end
end

% Convert 2D tessellation surface to 3D tessellation surface based on barycentric
% coordinate and cooresponding 3D grid cartesian coordiante
tessellation_target = cell(size(tessellation_close)); % cartesian coordinates of deployed grid surface
for i = 1:size(tessellation_close,1)
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