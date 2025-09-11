% Plot the angle_displacement curve in deployment process
%function angle_displacement(disp_min,disp_max)
disp_min = 0;
disp_max = 0.55*edgeLen;
% Define the parameters
alpha_1 = [];
alpha_2 = [];
delta = [];

% Define the initial alpha values for optimisation
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;

% Get the alpha_1 and alpha_2 from the displacement range
for i = 1:500
    pre_delta = disp_min + (disp_max - disp_min)/500*(i-1);
    [~,alpha_1_optimal,alpha_2_optimal] = triangle_unit(prev_alpha_1, prev_alpha_2, pre_delta, beta, edgeLen,l1,l4,t);
    delta(i) = pre_delta;
    alpha_1(i) = alpha_1_optimal;
    alpha_2(i) = alpha_2_optimal;
end

% Covert radius to degree
alpha_1 = alpha_1 * 180/pi;
alpha_2 = alpha_2 * 180/pi;

% Plot the angle-displacement curves
figure;
plot(delta, alpha_1, 'b', 'LineWidth', 1.5); % Alpha_1 vs Delta (blue curve)
hold on;
plot(delta, alpha_2, 'r', 'LineWidth', 1.5); % Alpha_2 vs Delta (red curve)
hold off;

% Add labels, title, and legend
xlabel('Displacement (\delta)', 'FontSize', 14);
ylabel('Angle (\alpha)', 'FontSize', 14);
title('Angle-Displacement Curve', 'FontSize', 16);
legend({'\alpha_1 vs \delta', '\alpha_2 vs \delta'}, 'Location', 'southeast', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 14)

%end
    