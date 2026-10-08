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
% % Plot figure using accuracy metric relative to FEM reference:
% %   accuracy(%) = |E - E_ref| / E * 100
% % E_ref can be scalar (single FEM reference) or same-size vector as N.
% E_ref = 2.021674068294217e+09;
% 
% acc_cluster = abs(E_cluster - E_ref) ./ max(abs(E_cluster), eps) * 100;
% acc_even = abs(E_even - E_ref) ./ max(abs(E_even), eps) * 100;
% 
% figure('Color','w');
% hold on; box on;
% semilogy(N, acc_cluster, '-o', 'Color',[0 0.45 0.74], ...
%     'MarkerFaceColor',[0 0.45 0.74], 'LineWidth',2, ...
%     'DisplayName','End-clustered');
% 
% semilogy(N, acc_even, '-s', 'Color',[0.85 0.33 0.10], ...
%     'MarkerFaceColor','none', 'LineWidth',2, ...
%     'DisplayName','Even-divided');
% xticks(0:20:100);          % or 4:1:100 for every integer
% xlim([min(N), max(N)]);
% xlabel('Number of segments  N', 'Interpreter','tex', 'FontSize',24);
% ylabel('Error  (%)', 'Interpreter','tex', 'FontSize',24);
% set(gca, 'FontName','Times New Roman', 'FontSize',30, 'LineWidth',2, 'YScale','log');
% legend('Location','northeast','Box','off', 'FontSize',18);
% grid off;
% axis square;
% 
% %% Compare semi-analytical result with FEA
% % Read FEA result from txt delta 0:0.1:0.55*edgeLen
% edgeLen = 15;
% l1 = edgeLen * 0.85;
% l4 = edgeLen * 0.05;
% t  = edgeLen * 0.015;
% beta = pi/40;
% fname = 'pi40_015.txt';  % <- change to your actual file
% raw = fileread(fname);
% % capture lines that contain exactly two numbers (delta, energy)
% expr = '(?:^|\r?\n)\s*([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)\s+([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)\s*(?=\r?\n|$)';
% tokens = regexp(raw, expr, 'tokens');
% if isempty(tokens)
%     error('No numeric (delta, energy) lines found. Check file formatting.');
% end
% % Convert tokens to numeric arrays
% nums = cellfun(@(t)[str2double(t{1}), str2double(t{2})], tokens, 'UniformOutput', false);
% data = vertcat(nums{:});
% delta_fem = data(:,1);             % displacement (m)
% U_fem     = data(:,2);             % total stored energy (J)
% E = 4.33e11;
% eps_fem = delta_fem ./ edgeLen;
% U_fem_max = max(U_fem);
% U_fem_norm = U_fem ./U_fem_max;
% 
% % Run semi-analyical result
% delta = 0:0.05:0.57*edgeLen;
% N = 10;
% [~, U_hbm] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta, t, N);
% eps_hbm = delta./edgeLen;
% U_hbm_max = max(U_hbm);
% U_hbm_norm = U_hbm ./ max(U_hbm_max, eps);
% 
% % Plot the configuration(beta)
% 
% % Plot reslt
% figure('Color','w'); 
% hold on; box on;
% 
% % semi-analytical
% plot(eps_hbm, U_hbm_norm, '-','Color',[0 0.45 0.74], 'LineWidth', 2, ...
%      'DisplayName','HBM');
% 
% % % FEM
% plot(eps_fem, U_fem_norm, '-', 'Color',[0.85 0.33 0.10],'LineWidth', 2, ...
%      'DisplayName','FEM');
% xlabel('Strain', 'Interpreter','tex', ...
%        'FontSize',24);
% ylabel('E/E_{max}', 'Interpreter','tex', ...
%        'FontSize',24);
% set(gca, 'FontName','Times New Roman','FontSize',30,'LineWidth',2); 
% legend('Location','northwest','Box','off', 'Fontsize',18);  
% grid off;
% axis square;
% 
% ylim([0, 1.05]);          % slightly larger than 1
% 
% 
% % %% sensetivity study on tilting angle beta and l2 or (l4)
% % % Scaning the titltign angle to get energy barrer and bistable delta
% % %delta = 0: 0.1: l3*(1-sin(beta))*1.12;
% % %[~, E_total] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta, t, 8);
% % beta = 0:pi/300:pi/15;
% % %l4 = edgeLen * (0.05:0.005:0.1);
% % energy_barrier = zeros(size(beta));
% % d_bist = zeros(size(beta));
% % for i = 1:size(beta,2)
% %     delta = 0: 0.05 : 0.61*edgeLen*(1-sin(beta(i)));
% %     [~, U_hbm] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta(i), t, 8);
% %     [energy_barrier(i), d_bist(i)] = find_bistable(U_hbm, delta);
% % end
% % 
% % % Plot results
% % figure()
% % box on; grid on;
% % plot(beta, energy_barrier, 'r-', 'LineWidth', 2);
% % xlabel('Cut inclination \beta', 'FontSize',16);
% % ylabel('Energy barrier', 'FontSize',16);
% % set(gca,'FontName','Times New Roman','FontSize',18,'LineWidth',1.0,'TickDir','out');
% % axis square;
% % 
% % 
% % figure()
% % box on; grid on;
% % plot(beta, d_bist, 'r-', 'LineWidth', 2);
% % xlabel('Cut inclination \beta', 'FontSize',16);
% % ylabel('bistable strain', 'FontSize',16);
% % set(gca,'FontName','Times New Roman','FontSize',18,'LineWidth',1.0,'TickDir','out');
% % axis square;
% % 
% % 
% % function [energy_barrier, d_bist] = find_bistable(U, delta)
% % dU = diff(U);
% % idx_max = find(dU < 0, 1, 'first');
% % if isempty(idx_max)
% %     energy_barrier = NaN;
% %     d_bist = NaN;
% %     fprintf('This unit is not bistable');
% % else
% %     U_max = U(idx_max);
% %     idx_bist = find(dU(idx_max+1:end)>0,1,'first');
% %     U_bist = U(idx_bist+idx_max);
% %     d_bist = delta(idx_bist+idx_max);
% %     energy_barrier = (U_max - U_bist)/U_max;
% % end
% % end
% 
% 
% %% Plot anisotropic and isotropic deployment
% % isotropic
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


%% Plot energy profile with different beta on same plot
edgeLen = 15;
l1      = edgeLen * 0.85;
l4      = edgeLen * 0.05;
t       = edgeLen * 0.015;
N       = 12;
delta   = 0:0.05:0.65*edgeLen;
eps_hbm = delta ./ edgeLen;

% Beta values from 0 to pi/20 (6 evenly spaced values)
n_beta   = 6;
beta_vec = linspace(0, pi/20, n_beta);

% Colour map: blue (beta=0) to red (beta=pi/20)
cmap = [linspace(0,0.85,n_beta)', linspace(0.45,0.10,n_beta)', linspace(0.74,0.10,n_beta)'];

results = struct();

for k = 1:n_beta
    [~, U_k] = deform_triangle_isotropic(delta, edgeLen, l1, l4, beta_vec(k), t, N);
    [eb_k, db_k] = find_bistable_iso(U_k, delta);

    results(k).beta           = beta_vec(k);
    results(k).delta          = delta;
    results(k).eps            = eps_hbm;
    results(k).U              = U_k;
    results(k).delta_bist     = db_k;
    results(k).eps_bist       = db_k / edgeLen;
    results(k).energy_barrier = eb_k;
    results(k).bistable       = ~isnan(db_k);

    fprintf('beta = %.4f  ->  eps_bist = %.4f  energy_barrier = %.4f\n', ...
        beta_vec(k), db_k/edgeLen, eb_k);
end

% Plot each curve up to 1.05 * bistable strain (or full range if not bistable)
figure('Color','w');
hold on; box on;

for k = 1:n_beta
    if results(k).bistable
        eps_plot_max = 1.05 * results(k).eps_bist;
    else
        eps_plot_max = max(eps_hbm);
    end
    mask = results(k).eps <= eps_plot_max;
    plot(results(k).eps(mask), results(k).U(mask), '-', ...
         'Color',     cmap(k,:), ...
         'LineWidth', 2, ...
         'DisplayName', sprintf('\\beta = %.4g', beta_vec(k)));
end

xlabel('Strain  \delta / L',          'Interpreter','tex', 'FontSize',24);
ylabel('Strain Energy (N{\cdot}mm)',   'Interpreter','tex', 'FontSize',24);
set(gca, 'FontName','Times New Roman', 'FontSize',24, 'LineWidth',2);
legend('Location','northwest', 'Box','off', 'FontSize',16);
grid off;
axis square;

save('beta_sweep_results.mat', 'results', 'beta_vec', 'edgeLen', 'l1', 'l4', 't', 'N');
fprintf('Saved beta sweep data to beta_sweep_results.mat\n');

function [energy_barrier, d_bist] = find_bistable_iso(U, delta)
    dU = diff(U);
    idx_max = find(dU < 0, 1, 'first');
    if isempty(idx_max)
        energy_barrier = NaN;  d_bist = NaN;  return;
    end
    idx_min_rel = find(dU(idx_max+1:end) > 0, 1, 'first');
    if isempty(idx_min_rel)
        energy_barrier = NaN;  d_bist = NaN;  return;
    end
    idx_min        = idx_min_rel + idx_max;
    energy_barrier = (U(idx_max) - U(idx_min)) / max(abs(U(idx_max)), eps);
    d_bist         = delta(idx_min);
end
