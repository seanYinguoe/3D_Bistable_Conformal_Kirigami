function result = main_unit_library(cfg)
%MAIN_UNIT_LIBRARY Sample, clean and refine a numerical unit-response library.
%   result = main_unit_library() runs the dense sweep in library_config.
%   Use small angle/beta arrays while developing a study; see docs/workflow.md.
%   Completed beta batches and raw detector diagnostics are saved separately
%   from interpolated/refined values. The supplied library is never replaced.

root = setup_project;
if nargin < 1
    cfg = library_config;
end
out = tempname(fullfile(root, 'results'));
mkdir(out);
save(fullfile(out, 'settings.mat'), 'cfg');
workers = 0;
if cfg.use_parallel
    assert(license('test', 'Distrib_Computing_Toolbox'), ...
        'Parallel Computing Toolbox is required.');
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('threads');
    end
    workers = pool.NumWorkers;
end

%% Sample ordered target-triangle angles
[a2_grid, a3_grid] = meshgrid(cfg.alpha2, cfg.alpha3);
a1_grid = pi - a2_grid - a3_grid;
tol = 1e-8;
mask = (a1_grid > -tol) & (a1_grid >= a2_grid - tol) & (a2_grid >= a3_grid - tol);
a1 = a1_grid(mask);
a2 = a2_grid(mask);
a3 = a3_grid(mask);
n_config = numel(a1);
n_beta = numel(cfg.beta);
assert(n_config > 0 && n_beta > 0, 'No valid angle/beta samples were selected.');

% Sine-law edge ratios, with the shortest edge scaled by cfg.scale.
edge = cfg.edge_length;
L1 = sin(a1) ./ sin(a3) * cfg.scale * edge;
L2 = sin(a2) ./ sin(a3) * cfg.scale * edge;
L3 = ones(n_config, 1) * cfg.scale * edge;
q3_y = (L1.^2 - L2.^2 - L3.^2) ./ (2*L3);
inside = L2.^2 - q3_y.^2;
valid = inside > 0;
q3_x = nan(n_config, 1);
q3_x(valid) = -sqrt(inside(valid));
l1 = edge * cfg.flank_length_ratio;
l4 = edge * cfg.flank_width_ratio;
t = edge * cfg.ligament_width_ratio;

%% Solve and save each beta batch
batch_dir = fullfile(out, 'beta_batches_refine');
mkdir(batch_dir);
eps_all = nan(n_beta*n_config, 1);
eta_all = nan(n_beta*n_config, 1);
for ib = 1:n_beta
    beta = cfg.beta(ib);
    eps_beta = nan(n_config, 1);
    eta_beta = nan(n_config, 1);
    detector_status = repmat("invalid_geometry", n_config, 1);
    max_constraint = nan(n_config, 1);
    parfor (i = 1:n_config, workers)
        if ~valid(i)
            continue;
        end
        q1 = [0 0];
        q2 = [0 -L3(i)];
        q3 = [q3_x(i) q3_y(i)];
        [eps_bist, eta_val, diagnostics] = bistability_analysis( ...
            l1, l4, beta, t, edge, q1, q2, q3, false);
        eps_beta(i) = eps_bist;
        eta_beta(i) = eta_val;
        detector_status(i) = string(diagnostics.status);
        max_constraint(i) = max(diagnostics.max_constraints);
    end
    save(fullfile(batch_dir, sprintf('beta_%03d.mat', ib)), ...
        'ib', 'beta', 'eps_beta', 'eta_beta', 'detector_status', 'max_constraint');
    rows = (ib-1)*n_config + (1:n_config);
    eps_all(rows) = eps_beta;
    eta_all(rows) = eta_beta;
    fprintf('Saved beta %d/%d batch.\n', ib, n_beta);
end

anisotropy_study_refine = table(repmat(a1, n_beta, 1), repmat(a2, n_beta, 1), ...
    repmat(a3, n_beta, 1), eps_all, eta_all, kron(cfg.beta(:), ones(n_config, 1)), ...
    'VariableNames', {'a1', 'a2', 'a3', 'eps_bist', 'eta_val', 'beta'});
save(fullfile(out, 'anisotropy_study_refine.mat'), 'anisotropy_study_refine');

%% Preserve raw data separately from cleaned and interpolated estimates
if cfg.postprocess
    [anisotropy_clean, reportClean] = clean_isolated_outliers(anisotropy_study_refine);
    [anisotropy_filled_refine, reportFill] = clean_and_fill_anisotropy(anisotropy_clean);
    finite = isfinite(anisotropy_filled_refine.eta_val) & isfinite(anisotropy_filled_refine.eps_bist);
    anisotropy_filter = anisotropy_filled_refine(finite, :);
    save(fullfile(out, 'anisotropy_filter.mat'), 'anisotropy_filter');
    save(fullfile(out, 'anisotropy_filled_refine.mat'), ...
        'anisotropy_filled_refine', 'anisotropy_clean', 'reportClean', 'reportFill');
else
    anisotropy_filter = anisotropy_study_refine;
end

%% Refine the detected bistable state toward the target endpoint
opts_refine = struct('edgeLen', edge, 'l1', l1, 'l4', l4, 't', t, ...
    'etaThreshold', 0.03, 'alphaTarget', 0.99, 'maxIter', 6, ...
    'epsFloor', 1e-4, 'epsCeil', 2.0, 'alphaTol', 1e-3);
if cfg.refine
    anisotropy_filter_refine = refine_eps_bist_to_endpoint(anisotropy_filter, opts_refine);
else
    anisotropy_filter_refine = anisotropy_filter;
end
save(fullfile(out, 'anisotropy_filter_refine.mat'), 'anisotropy_filter_refine', 'opts_refine');
result = struct('table', anisotropy_filter_refine, 'config', cfg, ...
    'output_dir', out, 'matlab_version', version);
save(fullfile(out, 'library_result.mat'), 'result');
fprintf('Saved unit-library results to %s\n', out);
end
