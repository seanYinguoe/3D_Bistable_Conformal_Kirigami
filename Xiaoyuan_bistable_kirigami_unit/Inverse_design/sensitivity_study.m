%% Sensitivity Study: Influence of t and l2 on Bistability
% Analyze how ligament width t and length l2 influence bistability
% We fix r4 = 0.05, r1 = 0.8
edgeLen = 15;  % Base triangle edge length

% Sweep ranges for l2 (via r2) and filament thickness t (via rt)
r2_vec = linspace(0.05, 0.3, 5);         % r2 = l2 / edgeLen
rt_vec = linspace(0.01, 0.04, 5);        % rt = t / edgeLen

% Initialize result storage
results = struct();
counter = 1;

% Loop over design parameters
for i = 1:length(r2_vec)
    for j = 1:length(rt_vec)
        r2 = r2_vec(i);
        rt = rt_vec(j);

        % Compute r3 from geometric constraint: r3 + 2*r2 = 0.85
        r3 = (0.85 - 2*r2);

        % Compute remaining parameters
        r4 = (1 - r3 - 3*r2)/3;  % from constraint
        r1 = r3 + 2*r2 + r4;

        % Skip infeasible configurations
        if r3 <= 0 || r4 <= 0 || r1 <= 0
            continue;
        end

        % Convert to physical lengths
        l1 = r1 * edgeLen;
        l2 = r2 * edgeLen;
        l3 = r3 * edgeLen;
        t = rt * edgeLen;

        % Target final shape (boundary nodes)
        scale = (1 + r3) * 1.05;
        q1 = [0, 0, 0] * scale;
        q3 = [0, -edgeLen, 0] * scale;
        q2 = [-sqrt(3)/2 * edgeLen, -0.5 * edgeLen, 0] * scale;

        % Run bistability analysis
        try
            [strain_bist, bistability] = bistability_analysis(l1, l2, l3, t, edgeLen, q1, q2, q3);
        catch
            strain_bist = NaN;
            bistability = NaN;
        end

        % Store result
        results(counter).r1 = r1;
        results(counter).r2 = r2;
        results(counter).r3 = r3;
        results(counter).l1 = l1;
        results(counter).l2 = l2;
        results(counter).l3 = l3;
        results(counter).t = t;
        results(counter).strain_bist = strain_bist;
        results(counter).bistability = bistability;
        results(counter).is_bistable = ~isnan(strain_bist);

        counter = counter + 1;
    end
end

% %% Visualize Feasible Design Space
figure;
feas_idx = [results.is_bistable];
scatter([results(feas_idx).r2], [results(feas_idx).t], 50, [results(feas_idx).bistability], 'filled');
colorbar;
xlabel('r_2 = l_2 / edgeLen');
ylabel('t = t / edgeLen');
title('Bistability Map (Normalized Geometry)');
grid on;

% %% Triangular Surface for Energy Barrier
% figure;
% tri = delaunay([results.r2], [results.r3]);
% trisurf(tri, [results.r2], [results.r3], [results.bistability], 'EdgeColor', 'none');
% xlabel('r_2'); ylabel('r_3'); zlabel('Energy Barrier');
% title('Energy Barrier Surface');
% colorbar;
% view(-30, 30);
% grid on;
% 
% %% Print Summary
% fprintf('Total tested configurations: %d\n', numel(results));
% fprintf('Bistable configurations: %d\n', sum([results.is_bistable]));
% 
% %% Save results to table
% T = struct2table(results);
% writetable(T, 'bistability_sweep.mat');
