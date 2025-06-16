%% Generate the deploying configuration and elastic energy
%tessellation_deploy = cell(size(tessellation));
load case1.mat;
n = 100;
E_rotation = cell(n+1,1);
E_bending  = cell(n+1,1);
tessellation_deploy = cell(n+1,1);
i = 0;
for alpha = 0:1/n:1
    i = i+1;
    for m = 1:size(tessellation,1)
        x_deploy = (1-alpha) * tessellation{m}(:,1) + alpha * tessellation_target{m}(:,1);
        y_deploy = (1-alpha) * tessellation{m}(:,2) + alpha * tessellation_target{m}(:,2);
        z_deploy = alpha * tessellation_target{m}(:,3);
        tessellation_deploy{i}{m,1} = [x_deploy, y_deploy, z_deploy];
        v_deploy = (1-alpha) * v_out + alpha * v_target;
    end
    [E_rotation{i},E_bending{i}] = energy_calculate(tessellation, tessellation_deploy{i}, v_deploy, f_out, x_out);
end

%% Calculate the rotation and bending energy during deployment
E_rt = zeros(n+1,1);
E_bt = zeros(n+1,1);
E_total = zeros(n+1,1);
for j = 1:n+1
    E_rt(j) = sum(E_rotation{j}(:));
    E_bt(j) = sum(E_bending{j}(:));
    E_total(j) = E_rt(j) + E_bt(j);
end

%% Plot the total energy during deployment curve
alpha_values = linspace(0, 1, n+1);

figure('Name','Energy Curve','Position',[200 200 900 600]);
hold on
plot(alpha_values, E_total, 'b-', 'LineWidth', 2, 'DisplayName', 'Total Energy');
plot(alpha_values, E_rt, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Rotation Energy');
plot(alpha_values, E_bt, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Bending Energy');
xlabel('Deployment Parameter \alpha');
ylabel('Energy');
title('Energy During Deployment');
legend('Location', 'northwest');
grid off
set(gca, 'FontSize', 16);
hold off