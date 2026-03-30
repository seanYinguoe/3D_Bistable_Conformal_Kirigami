% %% Run convergence study on N
% edgeLen = 15;
% E=4.3e11; b=1.0; t=0.015*edgeLen*(sqrt(3)/2);I = b*t^3/12;EI = E*I;
% N = 4:2:100;
% E_cluster = zeros(size(N));
% E_even = zeros(size(N));
% L0 = 1;
% A0=[0,0]; B0=[1,0];
% A1=[0,0]; B1=[1,0];
% vecA0 = [0,1];
% vecB0 = [0,1];
% vecA1=[cos(pi/10),sin(pi/10)]; vecB1=[cos(-pi/10),sin(-pi/10)];      % clamped end directions (normals)
% for i = 1:size(N,2)
%     Ni = N(i);
%     [E_cluster(i), ~, ~] = hbm_energy(L0, A1,B1, vecA1,vecB1, Ni, E,b, t, 'EndCluster',true);
%     [E_even(i), ~, ~] = hbm_energy(L0, A1,B1, vecA1,vecB1, Ni, E,b, t, 'EndCluster',false);
% end
% % for i = 1:size(N,2)
% %     Ni = N(i);
% %     [~, E_cluster(i)] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta, t, Ni); % segments cluster
% % end
% % Plot figure 
% figure; hold on; box on; grid on;
% plot(N, E_cluster, '-o', 'Color',[0 0.45 0.74], ...
%     'MarkerFaceColor',[0 0.45 0.74], 'LineWidth',1.6, ...
%     'DisplayName','End-clustered');
% 
% % Even
% plot(N, E_even, '-s', 'Color',[0.85 0.33 0.10], ...
%     'MarkerFaceColor','none', 'LineWidth',1.6, ...
%     'DisplayName','Even-divided');
% 
% xlabel('Number of segments  N','FontSize',20);
% ylabel('Total energy  E_{total}','FontSize',20);
% title('Convergence study on N','FontSize',20);
% legend('Location','northeast','Box','off', 'Fontsize',20);
% set(gca,'FontSize',20, ...
%     'FontName','Times New Roman',...
%     'LineWidth',1.0,'TickDir','out');
% axis square;  

%% Compare semi-analytical result with FEA
% Read FEA result from txt delta 0:0.1:0.55*edgeLen
edgeLen = 15;
l1 = edgeLen * 0.85;
l4 = edgeLen * 0.05;
t  = edgeLen * 0.015;
beta = pi/45;
fname = 'pi40_015.txt';  % <- change to your actual file
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
delta = 0:0.05:0.55*edgeLen;
N = 10;
[~, U_hbm] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta, t, N);
eps_hbm = delta./edgeLen;
U_hbm_max = max(U_hbm);
U_hbm_norm = U_hbm ./ max(U_hbm_max, eps);

% Plot the configuration(beta)

% Plot reslt
figure('Color','w'); 
hold on; box on;

% semi-analytical
plot(eps_hbm, U_hbm_norm, '-','Color',[0 0.45 0.74], 'LineWidth', 3, ...
     'DisplayName','HBM');

% % FEM
% plot(eps_fem, U_fem, '-', 'Color',[0.85 0.33 0.10],'LineWidth', 1, ...
%      'DisplayName','FEM');
xlabel('Strain', 'Interpreter','tex', ...
       'FontSize',24);
ylabel('U/U_{max}', 'Interpreter','tex', ...
       'FontSize',24);
set(gca, 'FontName','Times New Roman','FontSize',30,'LineWidth',2); 
%legend('Location','northwest','Box','off', 'Fontsize',18);  
grid off;
axis square;


% %% sensetivity study on tilting angle beta and l2 or (l4)
% % Scaning the titltign angle to get energy barrer and bistable delta
% %delta = 0: 0.1: l3*(1-sin(beta))*1.12;
% %[~, E_total] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta, t, 8);
% beta = 0:pi/300:pi/15;
% %l4 = edgeLen * (0.05:0.005:0.1);
% energy_barrier = zeros(size(beta));
% d_bist = zeros(size(beta));
% for i = 1:size(beta,2)
%     delta = 0: 0.05 : 0.61*edgeLen*(1-sin(beta(i)));
%     [~, U_hbm] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta(i), t, 8);
%     [energy_barrier(i), d_bist(i)] = find_bistable(U_hbm, delta);
% end
% 
% % Plot results
% figure()
% box on; grid on;
% plot(beta, energy_barrier, 'r-', 'LineWidth', 2);
% xlabel('Cut inclination \beta', 'FontSize',16);
% ylabel('Energy barrier', 'FontSize',16);
% set(gca,'FontName','Times New Roman','FontSize',18,'LineWidth',1.0,'TickDir','out');
% axis square;
% 
% 
% figure()
% box on; grid on;
% plot(beta, d_bist, 'r-', 'LineWidth', 2);
% xlabel('Cut inclination \beta', 'FontSize',16);
% ylabel('bistable strain', 'FontSize',16);
% set(gca,'FontName','Times New Roman','FontSize',18,'LineWidth',1.0,'TickDir','out');
% axis square;
% 
% 
% function [energy_barrier, d_bist] = find_bistable(U, delta)
% dU = diff(U);
% idx_max = find(dU < 0, 1, 'first');
% if isempty(idx_max)
%     energy_barrier = NaN;
%     d_bist = NaN;
%     fprintf('This unit is not bistable');
% else
%     U_max = U(idx_max);
%     idx_bist = find(dU(idx_max+1:end)>0,1,'first');
%     U_bist = U(idx_bist+idx_max);
%     d_bist = delta(idx_bist+idx_max);
%     energy_barrier = (U_max - U_bist)/U_max;
% end
% end


%% Plot anisotropic and isotropic deployment
% isotropic
% plot(alpha, E_total_iso, '-','Color',[0.85 0.33 0.10], 'LineWidth', 1, ...
%      'DisplayName','Isotropic');
% 
% % anisotropic
% plot(alpha, E_total_an, '-', 'Color',[0 0.45 0.74],'LineWidth', 1, ...
%      'DisplayName','Anisotropic');
% xlabel('Deployment', 'Interpreter','tex', ...
%        'FontSize',28);
% ylabel('Strain Energy(N/mm^2)', 'Interpreter','tex', ...
%        'FontSize',28);
% set(gca, 'FontName','Times New Roman','FontSize',28); 
% legend('Location','northwest','Box','off', 'Fontsize',28);  
% grid off;
% axis square;
