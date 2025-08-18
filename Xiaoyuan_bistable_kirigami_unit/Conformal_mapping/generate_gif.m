function generate_gif(tessellation, tessellation_target, T, modelname, filename)
% Number of frames for the GIF
numFrames = 50;

% Rotate the tessellation
for i = 1:size(tessellation,1)
    [tessellation{i},~] = model_rotate(modelname, tessellation{i}, ones(1,3));
    tessellation{i} = tessellation{i} - [T(:,[1,2]),0];% Move it to the central point   
end

% Create a figure for plotting
fig = figure('Name', 'Deployment Control', 'Position', [100, 100, 800, 600]);

for k = 0:numFrames
    alpha = k / numFrames;
    
    % Compute interpolated deployment shape
    tessellation_deploy = cell(size(tessellation));
    for m = 1:size(tessellation,1)
        x_deploy = (1-alpha) * tessellation{m}(:,1) + alpha * tessellation_target{m}(:,1);
        y_deploy = (1-alpha) * tessellation{m}(:,2) + alpha * tessellation_target{m}(:,2);
        z_deploy = alpha * tessellation_target{m}(:,3);
        tessellation_deploy{m} = [x_deploy, y_deploy, z_deploy];
    end

    % Clear previous plot and create new one
    clf(fig);
    set(fig, 'Color', 'w')
    colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255};
    hold on;
    for m = 1:size(tessellation,1)
        plot_triangle(tessellation_deploy{m},colour);
    end
    grid off;
    axis equal;
    axis off;
    set(gca, 'Color', 'w');      % Set axes background to white
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['Deployment Progress: ' num2str(alpha*100, '%.1f') '%']);
    view([45, 45]);
    hold off;

    % Capture the plot as an image
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Write to the GIF file
    if k == 0
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

close(fig);
end