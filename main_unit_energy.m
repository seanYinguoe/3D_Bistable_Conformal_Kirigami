function result = main_unit_energy(cfg)
%MAIN_UNIT_ENERGY Calculate unit deployment energy and detect bistability.
%   result = main_unit_energy() uses config/unit_config.m.
%   result = main_unit_energy(cfg) uses edited geometry and edge stretches.
%   Curves stop when the source detector reaches a classification or the
%   target configuration. MAT/CSV files retain unnormalised model energies.

root = setup_project;
if nargin < 1
    cfg = unit_config;
end
[q1, q2, q3, valid, reason] = scale_facs_to_q(cfg.edge_stretches, cfg.edge_length);
assert(valid, 'kirigami:UnitTarget', 'Invalid target triangle: %s', reason);
response = cell(numel(cfg.beta), 1);
strain = nan(size(cfg.beta));
eta = strain;
for k = 1:numel(cfg.beta)
    [strain(k), eta(k), response{k}] = bistability_analysis( ...
        cfg.flank_length, cfg.flank_width, cfg.beta(k), cfg.ligament_width, ...
        cfg.edge_length, q1, q2, q3, false);
    fprintf('beta %.4g: %s; detected strain %.4g, eta %.4g\n', ...
        cfg.beta(k), response{k}.status, strain(k), eta(k));
end

out = tempname(fullfile(root, 'results'));
mkdir(out);
result = struct('config', cfg, 'strain_bist', strain, 'eta', eta, ...
    'response', {response}, 'output_dir', out, 'matlab_version', version);
save(fullfile(out, 'unit_energy.mat'), 'result');
for k = 1:numel(response)
    curve = response{k};
    data = table(curve.alpha_history, curve.strain_history, curve.E_history, ...
        curve.exitflags, curve.max_constraints, 'VariableNames', ...
        {'deployment', 'strain', 'energy', 'exitflag', 'constraint_residual'});
    writetable(data, fullfile(out, sprintf('beta_%02d.csv', k)));
end
if cfg.make_plot
    fig = figure('Color', 'w', 'Position', [100 100 750 430]);
    hold on;
    for k = 1:numel(response)
        curve = response{k};
        plot(curve.strain_history, curve.E_history / max(max(curve.E_history), eps), ...
            'LineWidth', 1.6, 'DisplayName', sprintf('beta = %.2f deg (%s)', ...
            cfg.beta(k) * 180/pi, strrep(curve.status, '_', ' ')));
    end
    xlabel('Reference-edge strain');
    ylabel('Energy / maximum of each curve');
    legend('Location', 'best');
    box off;
    exportgraphics(fig, fullfile(out, 'unit_energy.png'), 'Resolution', 160);
end
fprintf('Saved energy curves and solver diagnostics to %s\n', out);
end
