function create_gif(tessellation, tessellation_target)
% Define the filename for the GIF
filename = 'deployment_progress.gif';

% Number of frames for the GIF
numFrames = 50;

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
    colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255};
    hold on;
    for m = 1:size(tessellation,1)
        plot_triangle(tessellation_deploy{m},colour);
    end
    grid off;
    axis equal;
    axis off;
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
