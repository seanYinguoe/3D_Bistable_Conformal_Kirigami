%% Read the data from the txt
fname = 'pi:40G0.txt';  % <- change to your actual file
raw = fileread(fname);

% capture lines that contain exactly two numbers (delta, energy)
expr = '(?:^|\r?\n)\s*([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)\s+([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)\s*(?=\r?\n|$)';
tokens = regexp(raw, expr, 'tokens');

if isempty(tokens)
    error('No numeric (delta, energy) lines found. Check file formatting.');
end

% Convert tokens to numeric arrays
nums = cellfun(@(t)[str2double(t{1}), str2double(t{2})], tokens, 'UniformOutput', false);
data = vertcat(nums{:});
delta_fem = data(:,1);             % displacement (m)
U_fem     = data(:,2);             % total stored energy (J)



%% Get the energy-displacement curve based on simulation results
E = 4.33e11;
eps_fem = delta_fem ./ edgeLen;

%% Get the energy-displacement curve based on rigid body model
% disp_min = 0;
% disp_max = edgeLen*0.57;
% beta = pi/40;
% [delta1, energy1] = energy_displacement(disp_min,disp_max,l1,l4,t,beta,edgeLen);
% A = pi/3-beta;
% l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(A); % length of filament(effective length)
% K_r = E*(t^3/12)/(l2);
% eps_r = delta1 / edgeLen;
% U_r = K_r * energy1';

%% Get the energy-displacement curve based on hybrid spring model
%[eps_hybrid, U_hybrid] = hybrid_energy_displacement(l1,l4,t,beta,edgeLen);
[eps_hbm, U_hbm] = hbm_energy_displacement(edgeLen*0.57,l1,l4,t,beta,edgeLen);

%% Plot: strain vs dimensionless energy
figure('Color','w'); 
hold on; box on;

% rotational spring model
% plot(eps_r, U_r, 'r-', 'LineWidth', 2, ...
%      'DisplayName','Rotational spring model');

% hybrid spring model
plot(eps_hbm, U_hbm, 'g-', 'LineWidth', 2, ...
     'DisplayName','Hybrid spring model');

% FEM
plot(eps_fem, U_fem, 'b-', 'LineWidth', 2, ...
     'DisplayName','FEM simulation');
xlabel('\epsilon (strain)', 'Interpreter','tex', ...
       'FontSize',16);
ylabel('U', 'Interpreter','tex', ...
       'FontSize',16);
title('Strain vs Energy', 'Interpreter','none', ...
       'FontSize',18);

set(gca, 'FontSize',16); 
legend('Location','best');  
grid on;