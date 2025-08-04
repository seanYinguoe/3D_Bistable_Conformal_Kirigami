%% Sensitivity Study with Normalized Ratios
% Analyze how l2 and l3 influence bistability using deform_triangle model

edgeLen = 15;  % Base triangle edge length

% Normalized parameter ranges(we fix r1 and to keep the geometry consistency)
r2_vec = linspace(0.05, 0.15, 10);     % l2/edgeLen 0.1
r3_vec = linspace(0.4, 0.7, 30);      % l3/edgeLen 0.55

% Filament width
t = 0.2;

% Initialize result storage
results = struct();
counter = 1;

% Loop over design parameters
for i = 1:length(r2_vec)
    for j = 1:length(r3_vec)
        r2 = r2_vec(i);
        r3 = r3_vec(j);
        r4 = (1 - r3 - 3*r2)/3;  % From geometric constraint
        r1 = r3 + 2*r2 + r4;
        if r4 <= 0
            continue;  % skip infeasible geometry
        end

        % Convert to physical lengths
        l1 = r1 * edgeLen;
        l2 = r2 * edgeLen;
        l3 = r3 * edgeLen;

        % Target final shape (boundary nodes)
        q1 = [0, 0, 0] * (1+r3) * 1.02;
        q3 = [0, -edgeLen, 0] * (1+r3) * 1.02;
        q2 = [-sqrt(3)/2 * edgeLen, -0.5 * edgeLen, 0] * (1+r3) * 1.02;

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
        results(counter).strain_bist = strain_bist;
        results(counter).bistability = bistability;
        results(counter).is_bistable = ~isnan(strain_bist);

        counter = counter + 1;
    end
end

%% Visualize Feasible Design Space
figure;
feas_idx = [results.is_bistable];
scatter([results(feas_idx).r2], [results(feas_idx).r3], 50, [results(feas_idx).bistability], 'filled');
colorbar;
xlabel('r_2 = l_2 / edgeLen');
ylabel('r_3 = l_3 / edgeLen');
title('Bistability Map (Normalized Geometry)');
grid on;

%% Triangular Surface for Energy Barrier
figure;
tri = delaunay([results.r2], [results.r3]);
trisurf(tri, [results.r2], [results.r3], [results.bistability], 'EdgeColor', 'none');
xlabel('r_2'); ylabel('r_3'); zlabel('Energy Barrier');
title('Energy Barrier Surface');
colorbar;
view(-30, 30);
grid on;

%% Print Summary
fprintf('Total tested configurations: %d\n', numel(results));
fprintf('Bistable configurations: %d\n', sum([results.is_bistable]));

%% Save results to table
T = struct2table(results);
writetable(T, 'bistability_sweep.csv');
