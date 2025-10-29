edgeLen = 15;
l1 = edgeLen * 0.85;
l4 = edgeLen * 0.05;
t  = edgeLen * 0.015;
beta = pi/40;
%% Run convergence study on N
% N = 5:1:20;
% E_cluster = zeros(size(N));
% for i = 1:size(N,2)
%     Ni = N(i);
%     [~, E_cluster(i)] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta, t, Ni); % segments cluster
% end
% % Plot figure
% figure; hold on; box on; grid off;
% plot(N, E_cluster, 'v', 'MarkerSize', 8, ...
%     'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', ...
%     'LineWidth', 1.2);
% xlabel('Number of Segments N','FontSize', 16);
% ylabel('Total Energy','FontSize', 16);
% title('Convergence Study on N');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

%% Compare semi-analytical result with FEA
% Read FEA result from txt delta 0:0.1:0.55*edgeLen
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
E = 4.33e11;
eps_fem = delta_fem ./ edgeLen;

% Run semi-analyical result
delta = 0:0.1:0.57*edgeLen;
N = 10;
[~, U_hbm] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta, t, N);
eps_hbm = delta./edgeLen;

% Plot reslt
figure('Color','w'); 
hold on; box on;

% semi-analytical
plot(eps_hbm, U_hbm, 'g-', 'LineWidth', 2, ...
     'DisplayName','HBM model');

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

%% sensetivity study on tilting angle beta, and thickness