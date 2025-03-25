% We do the parameters sweep on l1,l2,l3,t to plot how the geometric
% conditions affect energy gap

disp_min = 0;
disp_max = 0.5;

% We set the length of triangular as 1, here is the constraint condition: 
% l4 = l1 - 2*l2 - l3;
% L_ = 2*l4+l2+l1+t = 1;

num_l2 = 20; % Number of steps for l2
num_t = 20;  % Number of steps for t

l1 = 1;
l3 = 0.45;

maximum_energy = NaN(num_l2,num_t);
bistablity = NaN(num_l2,num_t);
for i = 1:num_l2+1     % The length of ligament
    for j = 1:num_t+1      % The thickness of ligament
        l2 = 0.15*1 + (0.25-0.15)/num_l2*(i-1);
        t  = 0.05*1 + (0.15-0.05)/num_t*(j-1);
        [maximum_energy(i,j), bistablity(i,j)] = energy_displacement(disp_min,disp_max,l1,l2,l3,t);
    end
end

% Create meshgrid for l1 and l2 values
l2_values = linspace(0.15, 0.25, num_l2+1);
t_values = linspace(0.05, 0.15, num_t+1);
[X, Y] = meshgrid(l2_values, t_values);

% Define finer meshgrid for interpolation
l2_fine = linspace(0.15, 0.25, 50);
t_fine = linspace(0.05, 0.15, 50);
[X_fine, Y_fine] = meshgrid(l2_fine, t_fine);

% Interpolate bistability values
bistability_fine = interp2(X, Y, bistablity, X_fine, Y_fine, 'cubic');

% Plot the surface
figure;
contourf(X_fine, Y_fine, bistability_fine, 20, 'LineColor', 'none'); % 20 contour levels
colorbar; % Add color bar to visualize magnitude of stability clearly
ylim([0.10,0.15])
xlabel('length of ligament','FontSize', 14);
ylabel('thickness of ligament','FontSize', 14);
zlabel('Stability','FontSize', 14);
set(gca, 'FontSize', 14)
grid off
colorbar; % Add color bar to visualize magnitude of energy
    