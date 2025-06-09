function tessellation_target = tessellation_deployment(tessellation,tessellation_target,T,modelname)
% Deploy a flattened tessellation onto the target tessellation
% using barycentric interpolation.
%% Move 2D closed tessellation configuration(should I rotate it as well?)
for i = 1:size(tessellation,1)
    [tessellation{i},~] = model_rotate(modelname, tessellation{i}, ones(1,3));
    tessellation{i} = tessellation{i} - [T(:,[1,2]),0];% Move it to the central point   
end

%% Deploy the flattened grids surface on to deployed target surface
% Create the figure and UI components
fig = figure('Name', 'Deployment Control', 'Position', [100 100 800 600]);
slider = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Position', [150 20 500 20], ...
    'Min', 0, 'Max', 1, 'Value', 0, ...
    'Callback', @(src,event) updatePlot(src, tessellation, tessellation_target));

updatePlot(slider, tessellation, tessellation_target);

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
% xlim([-90,90]);
% ylim([-90,90]);
% zlim([-10,60]);
axis off
xlabel('X');
ylabel('Y');
zlabel('Z');
title(['Deployment Progress: ' num2str(alpha*100, '%.1f') '%']);
view([45, 45]);
end
end