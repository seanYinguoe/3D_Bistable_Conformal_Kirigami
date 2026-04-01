%MAIN_2D Planar shape-morphing pipeline (2D) using kirigami inverse design.
% Run this script directly (not a function).
%
% This follows the same downstream workflow as main_3D, but replaces the
% 3D/BFF stage with a purely planar mapping stage:
%   conformal_mapping_2D(initialShape, targetShape, opts)
%
% Example (default square -> circle):
%   run('main_2D.m')
%
% Output variables are left in workspace; summary is in `out`.

clc;
clear;

% ----------------------- Options -----------------------
opts = struct();
opts.initialShape = 'square';
opts.targetShape = 'circle';
opts.edgeLen = 8;
opts.l1_ratio = 0.85;
opts.l4_ratio = 0.05;
opts.t_min_ratio = 0.035;
opts.t_max_ratio = 0.040;
opts.plot_flag = true;
opts.export_svg = false;
opts.svg_name = 'shape_morph_2d_pattern';
opts.compute_selected_energy = false;
opts.selected_unit_index = 1;
opts.compute_global_energy = false;

% Options for conformal_mapping_2D
opts.conformal = struct();
opts.conformal.n_boundary = 200;
opts.conformal.initial_scale = 100;
opts.conformal.center = [0 0];
% opts.conformal.mesh_h = 3.0; % optional fixed mesh size

% ---------------- 1) Load admissible anisotropy range -----------------
data_file = fullfile(repo_root, 'Inverse_design', 'anisotropy_filter_refine.mat');
if ~isfile(data_file)
    error('Missing admissible range file: %s', data_file);
end
S = load(data_file);
anisotropy_filter_refine = S.anisotropy_filter_refine;
mask = (anisotropy_filter_refine.eta_val > 0.10) & isfinite(anisotropy_filter_refine.eps_bist);
anisotropy_filter = anisotropy_filter_refine(mask,:);

% ---------------- 2) 2D conformal-like mapping stage ------------------
[v_initial_mesh, v_target_mesh, f_mesh, info_map] = conformal_mapping_2D( ...
    opts.initialShape, opts.targetShape, opts.conformal);

v_mesh = v_initial_mesh(:,1:2);

% ---------------- 3) Overlay and fit regular triangular grid ----------
[v_grid, f_grid, c_grid, i_grid, x_grid] = generate_overlay_grid(v_mesh, opts.edgeLen);

% fit_grid expects mesh-wise scale_facs values. For pure 2D domain fitting,
% use ones so geometry filtering is driven only by point-in-domain tests.
scale_facs_dummy = ones(size(f_mesh,1), 1);
[v_out, f_out, ~, i_out, ~, scale_area] = fit_grid( ...
    v_grid, f_grid, c_grid, i_grid, x_grid, v_mesh, f_mesh, scale_facs_dummy, opts.edgeLen);

% ---------------- 4) Interpolate target grid (purely planar) ----------
v_target_xy = map_points_barycentric_2D(v_out(:,1:2), v_initial_mesh(:,1:2), f_mesh, v_target_mesh(:,1:2));
v_initial = [v_out(:,1:2), zeros(size(v_out,1),1)];
v_target  = [v_target_xy, zeros(size(v_target_xy,1),1)];

% ---------------- 5) Edge-wise scale factors --------------------------
scale_facs = calculate_scale_facs(v_initial, v_target, f_out);
min_scale_factor = min(scale_facs(:), [], 'omitnan');
max_scale_factor = max(scale_facs(:), [], 'omitnan');

disp("Scale area min: " + num2str(min(scale_area, [], 'omitnan')) + ...
     ", max: " + num2str(max(scale_area, [], 'omitnan')));
disp("Min scale factor: " + num2str(min_scale_factor) + ...
     ", Max scale factor: " + num2str(max_scale_factor));

% ---------------- 6) Convert to anisotropy descriptors ----------------
[anisotropy_level, ~] = scale_facs_to_angles(scale_facs);

% ---------------- 7) Assign beta and thickness ------------------------
[opt_beta, bistability, assign_info] = assign_opt_beta(anisotropy_level, anisotropy_filter);

l1 = opts.edgeLen * opts.l1_ratio;
l4 = opts.edgeLen * opts.l4_ratio;
t_min = opts.edgeLen * opts.t_min_ratio;
t_max = opts.edgeLen * opts.t_max_ratio;
opt_t = assign_thickness(bistability, t_min, t_max);

params = [opts.edgeLen; l1; l4; t_min];

% ---------------- 8) Build tessellation in both states ----------------
tessellation = tessellated_triangle_initial(f_out, i_out, params, v_initial, opt_beta, opt_t);
tessellation_target = tessellated_triangle(f_out, i_out, params, v_target, opt_beta, opt_t);

% ---------------- 9) Visualize (similar style to main_3D) -------------
if opts.plot_flag
    figure('Name', 'Initial / Target Mesh (2D)', 'Color', 'w');
    subplot(1,2,1);
    patch('Vertices', v_initial_mesh(:,1:2), 'Faces', f_mesh, ...
        'FaceColor', 'none', 'EdgeColor', [0.2 0.2 0.2]);
    axis equal; axis off; title('Initial mesh');

    subplot(1,2,2);
    patch('Vertices', v_target_mesh(:,1:2), 'Faces', f_mesh, ...
        'FaceColor', 'none', 'EdgeColor', [0.2 0.2 0.2]);
    axis equal; axis off; title('Target mesh');

    figure('Name', 'Overlaid Grid on Initial Domain', 'Color', 'w');
    hold on;
    patch('Vertices', v_initial_mesh(:,1:2), 'Faces', f_mesh, ...
        'FaceColor', [0.95 0.95 0.98], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.5);
    plotgrid(f_out, v_out);
    axis equal; axis off; hold off;

    grid_deployment(v_target, v_initial, f_out);
    plot_edge_stretch(scale_facs, v_initial, f_out);

    colour = {'white', [0.9216 0.8863 0.4235], [0.7059 0.9608 0.4118], [0.9216 0.8863 0.4235]};
    figure('Name', 'Initial Tessellation', 'Color', 'w');
    hold on;
    for i = 1:size(f_out,1)
        plot_triangle(tessellation{i}, colour);
    end
    hold off; axis equal; axis off;

    tessellation_deployment(tessellation, tessellation_target);
end

% ---------------- 10) Optional energy evaluations ----------------------
selected_energy = struct();
if opts.compute_selected_energy
    i_sel = min(max(1, round(opts.selected_unit_index)), size(f_out,1));
    nD = 150;
    Nseg = 8;
    [q1, q2, q3, is_valid, reason] = scale_facs_to_q(scale_facs(i_sel,:), opts.edgeLen);
    selected_energy.index = i_sel;
    selected_energy.is_valid = is_valid;
    selected_energy.reason = reason;
    if is_valid
        [E_sel, alpha_sel] = deform_triangle_anisotropic(q1, q2, q3, opts.edgeLen, l1, l4, opt_beta(i_sel), opt_t(i_sel), nD, Nseg, false);
        selected_energy.E = E_sel;
        selected_energy.alpha = alpha_sel;
    else
        selected_energy.E = [];
        selected_energy.alpha = [];
    end
end

global_energy = struct();
if opts.compute_global_energy
    [E_total, E_unit, info_energy] = calculate_global_energy(params, opt_beta, opt_t, scale_facs, false);
    global_energy.E_total = E_total;
    global_energy.E_unit = E_unit;
    global_energy.info = info_energy;

    if opts.plot_flag
        alpha = linspace(0, 1, numel(E_total));
        Emax = max(E_total);
        E_norm = E_total ./ max(Emax, eps);
        figure('Name', 'Global Energy', 'Color', 'w');
        plot(alpha, E_norm, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
        xlabel('Deployment \xi'); ylabel('E/E_{max}');
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 16, 'LineWidth', 1.0);
        xlim([0 1]); ylim([0 1.05]);
        grid off;
    end
end

% ---------------- 11) Optional SVG export -----------------------------
svg_path = '';
svg_info = struct();
if opts.export_svg
    [svg_path, svg_info] = generate_svg(tessellation, opts.svg_name, false, [], 0.15);
    disp("SVG written to: " + svg_path);
end

% ---------------- Output bundle ----------------------------------------
out = struct();
out.opts = opts;
out.info_map = info_map;
out.v_initial_mesh = v_initial_mesh;
out.v_target_mesh = v_target_mesh;
out.f_mesh = f_mesh;
out.v_initial = v_initial;
out.v_target = v_target;
out.f_out = f_out;
out.i_out = i_out;
out.scale_facs = scale_facs;
out.opt_beta = opt_beta;
out.opt_t = opt_t;
out.assign_info = assign_info;
out.tessellation = tessellation;
out.tessellation_target = tessellation_target;
out.selected_energy = selected_energy;
out.global_energy = global_energy;
out.svg_path = svg_path;
out.svg_info = svg_info;
