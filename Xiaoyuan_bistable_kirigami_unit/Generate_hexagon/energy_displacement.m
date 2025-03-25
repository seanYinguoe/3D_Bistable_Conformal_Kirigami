% Plot the energy_displacement curve in deployment process, regarding
% filaments, flanks and inner triangle as rigid body, and the connecting
% points between filaments and flanks as rotational springs. Using Eular
% beam theory to calculate the energy
function [maximum_energy, bistability] = energy_displacement(disp_min,disp_max,l1,l2,l3,t)
% load simulation.mat

iter_max = 2500;
% Define the parameters
alpha_1 = size(iter_max);
alpha_2 = size(iter_max);
delta = size(iter_max);
energy = size(iter_max);

% Define the initial alpha values for optimisation
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;

% Get the alpha_1 and alpha_2 from the displacement range
for i = 1:iter_max
    pre_delta = disp_min + (disp_max - disp_min)/2000*(i-1);
    [~,alpha_1_optimal,alpha_2_optimal] = triangle_unit(prev_alpha_1, prev_alpha_2, pre_delta,l1,l2,l3,t);
    delta(i) = pre_delta;
    alpha_1(i) = alpha_1_optimal;
    alpha_2(i) = alpha_2_optimal;
    energy(i) = 3*(1/2*(alpha_1(i)-prev_alpha_1).^2 + 1/2*(alpha_2(i)-prev_alpha_2).^2);
end

% Get the maximum energy
maximum_energy = max(energy);

% Numerical differentiation to find critical points
dE_dDelta = diff(energy) ./ diff(delta); % First derivative of energy w.r.t. delta

bistable_energy = maximum_energy;
bistability = NaN;
for i = 1:length(dE_dDelta)-1
    if dE_dDelta(i) < 0 && dE_dDelta(i+1) > 0
        bistable_energy = energy(i); % Store bistable energy value
        bistability = (maximum_energy-bistable_energy)/maximum_energy;
        break; 
    end
end

% Plot the energy-displacement curves
% figure;
% hold on;
% plot(delta, energy,'black', 'LineWidth', 1.5); % Energy vs Delta
% %plot(simulation(:,1), simulation(:,2)/3e6,'r', 'LineWidth', 1.5); % Energy vs Delta
% 
% %Add labels, title, and legend
% xlabel('Displacement (\delta)', 'FontSize', 18);
% ylabel('E/K', 'FontSize', 18);
% title('Energy-Displacement Curve', 'FontSize', 18);
% legend({'Theory'}, 'Location', 'southeast', 'FontSize', 18);
% grid on;
% set(gca, 'FontSize', 18)

end
    