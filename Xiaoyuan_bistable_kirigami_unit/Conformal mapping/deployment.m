% Deploy a flattened surface onto the deployed surface (conformal mapping)
% using interpolation and analysis.
% Input: flattened_surface (vertices, faces)
%        deployed_surface (vertices, faces)

% Example Input (Replace with your actual data)
flattened_surface = obj_2D.vt;  % Flattened 2D vertex positions
deployed_surface = obj_3D.v;    % Deployed 3D vertex positions
faces = obj_2D.f.v;             % Face connectivity


% Move the deployed surface and flattned surface onto the same platform


num = 20; % Number of steps in deployment
alpha = linspace(0,1,num); % Deployment parameter

% Step 4: Visualization Loop for Deployment Animation
figure;
for k = 1:num
    % Compute interpolated deployment shape
    X_deploy = (1-alpha(k)) * flattened_surface(:,1) + alpha(k) * deployed_surface(:,1);
    Y_deploy = (1-alpha(k)) * flattened_surface(:,2) + alpha(k) * deployed_surface(:,2);
    Z_deploy = alpha(k) * deployed_surface(:,3);  % Gradually lift into 3D

    V_deploy = [X_deploy, Y_deploy, Z_deploy];

    % Plotting
    patch('Vertices', V_deploy, 'Faces', faces, ...
          'FaceColor', 'cyan', 'EdgeColor', 'blue', 'FaceAlpha', 0.6); % Transparent surface
    
    hold on;
    plot3(X_deploy, Y_deploy, Z_deploy, 'k.'); % Mesh points
    hold off;
    
    axis equal; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['Deployment Step: ' num2str(k) ' / ' num2str(num)]);
    
    % Rotate view dynamically
    view(mod(k * 5, 360), 30);  
    
    colormap jet; shading interp;
    drawnow;
end

disp('Deployment animation complete.');
