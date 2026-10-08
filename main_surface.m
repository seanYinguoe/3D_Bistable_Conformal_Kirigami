function result = main_surface(cfg)
%MAIN_SURFACE Map a target surface to a flat kirigami cut pattern.
%   result = main_surface() uses config/surface_config.m.
%   result = main_surface(cfg) uses an edited configuration struct.
%   Each run saves its settings, geometry, assignment diagnostics and SVG
%   in a new results/ folder. Inspect diagnostics before fabrication.

root = setup_project;
if nargin < 1
    cfg = surface_config;
end
assert(exist('fmincon', 'file') == 2 && exist('fsolve', 'file') == 2, ...
    'Optimization Toolbox is required.');
validateattributes(cfg.edge_length, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});

%% Read the supplied curved mesh and its planar UV map
[input_dir, name, ext] = fileparts(cfg.mesh);
obj = readObj(input_dir, [name ext]);
obj.vt = vertice_sort(obj.vt, obj.f.v, obj.f.vt);
faces = obj.f.v;
uv = obj.vt;
flat_area = triangle_area_2D(uv, faces);
curved_area = triangle_area_3D(obj.v, faces);
assert(all(flat_area > 0) && all(curved_area > 0), 'Input mesh has degenerate faces.');
area_scale = sqrt(curved_area ./ flat_area);

%% Fit a regular grid, then map and reparameterise it on the target
[v, f, centres, orientation, neighbours] = generate_overlay_grid(uv, cfg.edge_length);
[v, f, centres, orientation, neighbours, grid_area_scale] = fit_grid( ...
    v, f, centres, orientation, neighbours, uv, faces, area_scale, cfg.edge_length);
assert(~isempty(f), 'No grid cells fit the mesh; reduce edge_length.');
[flat_mesh, curved_mesh] = mesh_deployment(obj);
fprintf('Reparameterising %d grid cells...\n', size(f, 1));
[flat, target, max_angle, initial_target, reparam] = reparameterization( ...
    v, f, flat_mesh, curved_mesh, faces);

%% Match directional stretch demands to the supplied unit library
saved = load(cfg.library, 'anisotropy_filter_refine');
library = saved.anisotropy_filter_refine;
library = library(library.eta_val > cfg.library_eta_min & isfinite(library.eps_bist), :);
assert(~isempty(library), 'The selected library has no admissible rows.');
stretch = calculate_scale_facs(flat, target, f);
rescale = (1 + max(library.eps_bist)) / max(min(stretch, [], 2));
[target, stretch, rescale_info] = rescale_target_edges(flat, f, obj, target, rescale, struct());
[level, ~] = scale_facs_to_angles(stretch);
mesh = triangulation(f, flat(:, 1:2));
rim_nodes = unique(freeBoundary(mesh));
is_rim = any(ismember(f, rim_nodes), 2);
assignment_opts = cfg.assignment;
assignment_opts.rimMask = is_rim;
[beta, eta, assignment] = assign_opt_beta(level, library, assignment_opts);

%% Build compact and deployed unit geometries
edge_length = cfg.edge_length;
thickness = assign_thickness(eta, edge_length * cfg.ligament_width_ratios(1), ...
    edge_length * cfg.ligament_width_ratios(2), cfg.thickness_threshold);
params = [edge_length; edge_length * cfg.flank_length_ratio; ...
    edge_length * cfg.flank_width_ratio; edge_length * cfg.ligament_width_ratios(1)];
compact = tessellated_triangle_initial(f, orientation, params, flat, beta, thickness);
deployed = {};
geometry_info = [];
if cfg.build_deployed
    [deployed, geometry_info] = tessellated_triangle(f, orientation, params, target, beta, thickness);
end

%% Save a self-contained run and explicit assignment-quality information
out = tempname(fullfile(root, 'results'));
mkdir(out);
result = struct('config', cfg, 'output_dir', out, 'matlab_version', version, ...
    'faces', f, 'flat_vertices', flat, 'target_vertices', target, ...
    'initial_target', initial_target, 'orientation', orientation, ...
    'centres', centres, 'neighbours', {neighbours}, 'area_scale', grid_area_scale, ...
    'edge_stretch', stretch, 'angle_level', level, 'max_angle', max_angle, ...
    'reparameterization', reparam, 'rescaling', rescale_info, 'params', params, ...
    'beta', beta, 'eta_prediction', eta, 'thickness', thickness, 'rim', is_rim, ...
    'assignment', assignment, 'compact', {compact}, 'deployed', {deployed}, ...
    'geometry_diagnostics', geometry_info);
% The inherited "eps_matched" flag means a nearest match was selected.
% It does not by itself mean the mismatch is within the requested tolerance.
result.assignment_within_tolerance = ~is_rim & assignment.flag == "eps_matched" & ...
    assignment.beta_error <= assignment_opts.strainTol;
unit_table = table((1:size(f, 1))', is_rim, beta, thickness, eta, assignment.flag, ...
    assignment.beta_error, result.assignment_within_tolerance, ...
    'VariableNames', {'unit', 'rim', 'beta', 'ligament_width', 'eta_predicted', ...
    'assignment_flag', 'strain_mismatch', 'within_strain_tolerance'});
writetable(unit_table, fullfile(out, 'unit_assignments.csv'));
if cfg.export_svg
    [result.svg_path, result.svg_info] = generate_svg( ...
        compact, 'cut_pattern.svg', false, [], cfg.fillet_radius, out);
end
save(fullfile(out, 'surface_result.mat'), 'result');
if cfg.evaluate_energy
    [result.energy_total, result.energy_units, result.energy_info] = calculate_global_energy( ...
        params, beta, thickness, stretch, false, cfg.energy);
    save(fullfile(out, 'surface_result.mat'), 'result');
end
if cfg.make_plot
    plot_surface_result(result);
end
fprintf('%d cells; %d rim; %d/%d interior assignments within strain tolerance.\n', ...
    size(f, 1), nnz(is_rim), nnz(result.assignment_within_tolerance), nnz(~is_rim));
fprintf('Saved surface results to %s\n', out);
if any(~result.assignment_within_tolerance & ~is_rim)
    warning('kirigami:AssignmentMismatch', ...
        'Some interior units need review; inspect unit_assignments.csv.');
end
end
