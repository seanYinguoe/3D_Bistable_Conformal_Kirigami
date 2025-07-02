%% Sensitivity Study with Normalized Ratios
edgeLen = 15;  % Total edge length

% Define grid for independent ratios
r2_vec = linspace(0, 0.3, 40);        % l2/edgeLen (0 to 0.3)
r3_vec = linspace(0.2, 0.8, 40);      % l3/edgeLen (0.2 to 0.8)

% Initialize results storage
results = struct();
counter = 1;

%% Parameter Sweep
for i = 1:length(r2_vec)
    for j = 1:length(r3_vec)
        % Calculate normalized ratios
        r2 = r2_vec(i)
        r3 = r3_vec(j)
        r4 = 1 - r3 - 3*r2;  % Derived from constraint

        % Check feasibility
        if r4 <= 0 || r2 <= 0 || r3 <= 0
            continue;  % Skip infeasible designs
        end

        % Convert to physical lengths
        l2 = r2 * edgeLen; 
        l3 = r3 * edgeLen;
        l4 = r4 * edgeLen;

        % Set displacement range (proportional to system size)
        delta_range = linspace(0, 1.5*l3, 2000);

        % Run bistability analysis
        [disp_bist, bistability] = bistability_analysis(l2, l3, delta_range);

        % Store results
        results(counter).r2 = r2;
        results(counter).r3 = r3;
        results(counter).r4 = r4;
        results(counter).l2 = l2;
        results(counter).l3 = l3;
        results(counter).l4 = l4;
        results(counter).disp_bist = disp_bist;
        results(counter).bistability = bistability;
        results(counter).is_bistable = ~isempty(disp_bist);

        counter = counter + 1;
    end
end

%% Visualization: Feasible Design Space(bistable, monostable)
figure;
% Plot feasible region
r2_grid = linspace(0, 0.3, 100);
r3_grid = linspace(0.2, 0.8, 100);
[R2, R3] = meshgrid(r2_grid, r3_grid);
R4 = 1 - R3 - 3*R2;
feasible = R4 > 0 & R2 > 0 & R3 > 0;

contourf(R2, R3, double(feasible), [0.5, 0.5]);
colormap([1 0.8 0.8; 0.8 1 0.8]);
hold on;

% Highlight bistable designs
bistable_idx = [results.is_bistable];
scatter([results(bistable_idx).r2], [results(bistable_idx).r3], 50, ...
    [results(bistable_idx).bistability], 'filled');
colorbar;
xlabel('r_2 = l_2 / edgeLen');
ylabel('r_3 = l_3 / edgeLen');
title('Bistability in Normalized Design Space');
grid on;

%% Visualization: Bistability Response
figure;
% Triangular contour plot
tri = delaunay([results.r2], [results.r3]);
trisurf(tri, [results.r2], [results.r3], [results.bistability], ...
    'EdgeColor', 'none', 'FaceColor', 'interp');
xlabel('r_2');
ylabel('r_3');
zlabel('Bistability (Energy Barrier)');
title('Bistability Sensitivity');
colorbar;
view(-30, 30);
grid on;

%% Statistical Analysis
if any([results.is_bistable])
    bistable_results = results([results.is_bistable]);

    fprintf('Sensitivity Analysis Results (edgeLen = %.2f):\n', edgeLen);
    fprintf('================================================\n');
    fprintf('Bistable designs: %d/%d\n', sum([results.is_bistable]), numel(results));
    fprintf('Average bistability: %.4f ± %.4f\n', ...
        mean([bistable_results.bistability]), std([bistable_results.bistability]));

    % Correlation analysis
    [corr_r2, pval_r2] = corr([bistable_results.r2]', [bistable_results.bistability]');
    [corr_r3, pval_r3] = corr([bistable_results.r3]', [bistable_results.bistability]');

    fprintf('\nCorrelation with r_2: %.4f (p=%.4f)\n', corr_r2, pval_r2);
    fprintf('Correlation with r_3: %.4f (p=%.4f)\n', corr_r3, pval_r3);
end
