%% Sensitivity Study: Influence of t and l2 on Bistability
% Analyze how ligament width t and length l2 influence bistability
% We fix r4 = 0.05, r1 = 0.85
clc; clear;

edgeLen = 15;  % Base triangle edge length

% Sweep ranges for beta(beta is the tilting angle)
%beta_vec = linspace(0, pi/15, 200);         % r2 = l2 / edgeLen

beta_vec = pi/40;

results = struct();
counter = 1;

% Loop over design parameters
for i = 1:length(beta_vec)
    % Convert to physical lengths
    l1 = 0.85 * edgeLen;
    l4 = 0.05 * edgeLen;
    t = 0.02 * edgeLen;
    beta = beta_vec(i);
    % Calculate te length of inner triangle
    A = pi/3-beta;
    l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(A); % length of filament
    l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) ...
        - cos(pi/3 - beta) .* l4;
    l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
    l3 = l6 - l5 - 1.5*l2 ...
        - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) );

    % Target final shape (boundary nodes)
    scale = (1 + l3/edgeLen) * 1.1;
    q1 = [0, 0, 0] * scale;
    q3 = [0, -edgeLen, 0] * scale;
    q2 = [-sqrt(3)/2 * edgeLen, -0.5 * edgeLen, 0] * scale;

    % Run bistability analysis
    try
        [strain_bist, bistability] = bistability_analysis(l1, l4, beta, t, edgeLen, q1, q2, q3);
    catch
        strain_bist = NaN;
        bistability = NaN;
    end

    % Store result
    results(counter).beta = beta;
    results(counter).l1 = l1;
    results(counter).l4 = l4;
    results(counter).strain_bist = strain_bist;
    results(counter).bistability = bistability;
    results(counter).is_bistable = ~isnan(strain_bist);
    counter = counter + 1;
end

%% Save results to T
T = struct2table(results);
writetable(T, 'bistability_sweep.csv'); 


%% Visualize Feasible Design Space
% ==== Extract directly from T ====
beta       = T.beta;         % tilt angle
strain_bist= T.strain_bist;  % bistable strain
bistability= T.bistability;  % energy barrier

% ==== Clean and sort ====
mask = ~isnan(beta) & ~isnan(strain_bist) & ~isnan(bistability);
beta        = beta(mask);
strain_bist = strain_bist(mask);
bistability = bistability(mask);

[beta, idx] = sort(beta);
strain_bist = strain_bist(idx);
bistability = bistability(idx);

% ==== Plot 1: bistable strain vs beta ====
figure; hold on; grid on;
plot(beta, strain_bist, '-', 'LineWidth', 2);
xlabel('\beta (rad)');
ylabel('Bistable strain \epsilon^*');
title('Bistable strain vs \beta');

% ==== Plot 2: bistability (energy barrier) vs beta ====
figure; hold on; grid on;
plot(beta, bistability, '-', 'LineWidth', 2);
xlabel('\beta (rad)');
ylabel('Energy barrier');
title('Bistability vs \beta');
