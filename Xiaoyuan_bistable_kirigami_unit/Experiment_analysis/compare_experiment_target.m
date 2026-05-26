% compare_experiment_target.m
% Extract upper silhouette from experiment photo and rendered target image,
% normalise both to [-1,1] x [0,1], and plot the comparison.
%
% Experiment : white kirigami dome on black background
%   → binarise (Otsu), morphological closing to fill gaps between pieces,
%     upper boundary extracted directly (no convex hull)
% Target     : gray hemisphere on white background
%   → threshold near white to capture full dome, fill holes

clc; clear; close all;

%% ── 1. SETTINGS ─────────────────────────────────────────────────────────
model = 'doubledome';
exp_path    = strcat('Figures/',model,'_exp.jpg');
tar_path    = strcat('Figures/',model,'_tar.jpg');

n_compare = 300;   % uniform width samples for comparison

% ── Per-model tuning ──────────────────────────────────────────────────────
% close_radius : morphological closing disk (px) — bridge kirigami gaps
% sg_span      : Savitzky-Golay window for boundary smoothing (odd, > sg_order)
% reg_scales   : search range for horizontal scale registration
%                (scale < 1 stretches exp wider; use tight range when shape
%                 already matches, wider range when photo foreshortens dome)
% reg_shifts   : search range for lateral shift registration
sg_order = 3;
switch model
    case 'doubledome'
        close_radius = 60;
        sg_span      = 101;
        reg_scales   = linspace(0.95, 1.05, 11);   % shift-only effectively
        reg_shifts   = linspace(-0.20, 0.20, 81);
    case 'quarterdome'
        close_radius = 60;
        sg_span      = 101;
        reg_scales   = linspace(0.60, 1.05, 91);   % allow wide stretch
        reg_shifts   = linspace(-0.30, 0.30, 121);
    otherwise
        close_radius = 60;
        sg_span      = 101;
        reg_scales   = linspace(0.90, 1.05, 16);
        reg_shifts   = linspace(-0.20, 0.20, 81);
end

%% ── 2. EXPERIMENT IMAGE ──────────────────────────────────────────────────
% White kirigami structure on dark background.
% Strategy: binarise → close gaps → upper boundary.

img_exp  = imread(exp_path);
gray_exp = rgb2gray(img_exp);
bw_exp   = imbinarize(gray_exp);              % white=1, black=0

% Remove any white blobs that touch the image border (e.g. photo frame/mat).
% The kirigami dome sits away from the edges, so it is unaffected.
bw_exp = imclearborder(bw_exp);

% Close gaps between kirigami pieces to merge them into one solid shape
se          = strel('disk', close_radius);
bw_closed   = imclose(bw_exp, se);

% Keep only the largest blob (the dome)
bw_dome     = bwareafilt(bw_closed, 1);

% NOTE: bwconvhull is intentionally NOT used here. For non-convex shapes
% such as double-dome, the convex hull bridges the saddle and erases the dip.
% The morphological closing above is sufficient to fill inter-piece gaps.

[nrows_e, ~] = size(bw_dome);

% Upper boundary: topmost foreground pixel per column
cols_e = find(any(bw_dome, 1));
x_e    = zeros(1, numel(cols_e));
y_e    = zeros(1, numel(cols_e));
for k  = 1:numel(cols_e)
    rows   = find(bw_dome(:, cols_e(k)));
    x_e(k) = cols_e(k);
    y_e(k) = nrows_e - min(rows);         % flip: row 1 = top of image
end

% Smooth the pixel-grid staircase noise while preserving peaks and valley
y_e = sgolayfilt(y_e, sg_order, sg_span);

% Normalise
xc_e = (x_e(1) + x_e(end)) / 2;
xr_e = (x_e(end) - x_e(1)) / 2;
w_exp_raw = (x_e - xc_e) / xr_e;
h_exp_raw = (y_e - min(y_e)) / (max(y_e) - min(y_e));

%% ── 3. TARGET IMAGE ──────────────────────────────────────────────────────
% Gray hemisphere on white background with soft rendered edges.
% imclearborder fails because the dome base touches the image border.
% Strategy: Canny edge detection finds the dome outline regardless of
% absolute brightness, imclose bridges small gaps, imfill fills interior.

img_tar  = imread(tar_path);
gray_tar = rgb2gray(img_tar);

edges_tar = edge(gray_tar, 'Canny');          % detect dome boundary
edges_tar = imclose(edges_tar, strel('disk', 4));  % close small gaps
bw_tar    = imfill(edges_tar, 'holes');       % fill enclosed region → solid dome
bw_tar    = bwareafilt(bw_tar, 1);            % keep largest blob

[nrows_t, ~] = size(gray_tar);

% Upper boundary: topmost foreground pixel per column
cols_t = find(any(bw_tar, 1));
x_t    = zeros(1, numel(cols_t));
y_t    = zeros(1, numel(cols_t));
for k  = 1:numel(cols_t)
    rows   = find(bw_tar(:, cols_t(k)));
    x_t(k) = cols_t(k);
    y_t(k) = nrows_t - min(rows);
end

% Normalise
xc_t = (x_t(1) + x_t(end)) / 2;
xr_t = (x_t(end) - x_t(1)) / 2;
w_tar_raw = (x_t - xc_t) / xr_t;
h_tar_raw = (y_t - min(y_t)) / (max(y_t) - min(y_t));

%% ── 4. RESAMPLE ONTO COMMON GRID ─────────────────────────────────────────

w_grid   = linspace(-1, 1, n_compare);
h_exp_rs = interp1(w_exp_raw, h_exp_raw, w_grid, 'linear', NaN);
h_tar_rs = interp1(w_tar_raw, h_tar_raw, w_grid, 'linear', NaN);

% Final Gaussian smoothing pass on the resampled experiment curve to remove
% any residual staircase noise introduced during interpolation.
h_exp_rs = smoothdata(h_exp_rs, 'gaussian', 20);

% ── Optimal registration: shift + scale ───────────────────────────────────
% Corrects two independent biases:
%   b  (shift) : photo not perfectly centred laterally
%   a  (scale) : morphological closing expands blob base, compressing the
%                dome profile in normalised coordinates (a < 1 stretches it)
% Transform applied: h_reg(w) = h_exp( a*w + b )
% Registration objective: blend of RMSE and max error.
% alpha=1 → pure RMSE; alpha=0 → pure max error (Chebyshev).
% A value of 0.4 reduces the max error more aggressively than pure RMSE.
reg_alpha = 0.4;

scales    = reg_scales;
shifts    = reg_shifts;
best_cost = Inf;
best_a    = 1;
best_b    = 0;
for a = scales
    for b = shifts
        h_reg     = interp1(w_grid, h_exp_rs, a*w_grid + b, 'linear', NaN);
        valid_reg = ~isnan(h_reg) & ~isnan(h_tar_rs);
        if sum(valid_reg) < 50, continue; end
        d = h_reg(valid_reg) - h_tar_rs(valid_reg);
        cost = reg_alpha * sqrt(mean(d.^2)) + (1 - reg_alpha) * max(abs(d));
        if cost < best_cost
            best_cost = cost;
            best_a    = a;
            best_b    = b;
        end
    end
end
h_exp_rs = interp1(w_grid, h_exp_rs, best_a*w_grid + best_b, 'linear', NaN);
fprintf('Registration: scale a = %.4f, shift b = %.4f\n', best_a, best_b);

% ── Piecewise peak-landmark warp ──────────────────────────────────────────
% Corrects residual asymmetric peak shift (e.g. camera off-axis).
% For each target width position w_t, the lookup in exp space is:
%   w_exp_lookup = interp1([−1, lp_t, rp_t, 1], [−1, lp_e, rp_e, 1], w_t)
% i.e. w_src uses TARGET landmarks, w_dst uses EXP landmarks.
min_peak_h   = 0.80;
min_peak_sep = 0.20;

[~, locs_t] = findpeaks(h_tar_rs, w_grid, ...
    'MinPeakHeight', min_peak_h, 'MinPeakDistance', min_peak_sep);
[~, locs_e] = findpeaks(h_exp_rs, w_grid, ...
    'MinPeakHeight', min_peak_h, 'MinPeakDistance', min_peak_sep);

if numel(locs_t) >= 2 && numel(locs_e) >= 2
    locs_t = sort(locs_t);   locs_e = sort(locs_e);
    lp_t = locs_t(1);   rp_t = locs_t(end);
    lp_e = locs_e(1);   rp_e = locs_e(end);

    % Also find the valley (minimum between the two peaks) in both curves.
    % Adding it as a 3rd landmark gives 4 independent segments, correcting
    % asymmetric horizontal compression around the saddle region.
    mask_t  = w_grid > lp_t & w_grid < rp_t;
    wbt     = w_grid(mask_t);
    [~,ivt] = min(h_tar_rs(mask_t));
    vl_t    = wbt(ivt);

    mask_e  = w_grid > lp_e & w_grid < rp_e;
    wbe     = w_grid(mask_e);
    [~,ive] = min(h_exp_rs(mask_e));
    vl_e    = wbe(ive);

    % 5-point piecewise linear map (target space → exp lookup space):
    %   [-1, lp_t, vl_t, rp_t, 1]  →  [-1, lp_e, vl_e, rp_e, 1]
    w_src    = [-1,  lp_t,  vl_t,  rp_t,  1];
    w_dst    = [-1,  lp_e,  vl_e,  rp_e,  1];
    w_lookup = interp1(w_src, w_dst, w_grid, 'linear', 'extrap');
    h_exp_rs = interp1(w_grid, h_exp_rs, w_lookup, 'linear', NaN);
    fprintf('Peak warp  : L %.3f→%.3f  V %.3f→%.3f  R %.3f→%.3f\n', ...
        lp_e, lp_t, vl_e, vl_t, rp_e, rp_t);
else
    fprintf('Peak warp  : skipped (found %d tar peaks, %d exp peaks)\n', ...
        numel(locs_t), numel(locs_e));
end

%% ── 5. ERROR METRICS ─────────────────────────────────────────────────────

valid   = ~isnan(h_exp_rs) & ~isnan(h_tar_rs);
rmse    = sqrt(mean((h_exp_rs(valid) - h_tar_rs(valid)).^2));
max_err = max(abs(h_exp_rs(valid) - h_tar_rs(valid)));

fprintf('RMSE      = %.4f\n', rmse);
fprintf('Max error = %.4f\n', max_err);

%% ── 6. FIGURE 1: Normalised height comparison (line plot) ───────────────

figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 14 12]);
hold on; box on;

plot(w_grid, h_tar_rs, '-',  'Color', [0    0    0   ], 'LineWidth', 2.0, ...
    'DisplayName', 'Tar.');
plot(w_grid, h_exp_rs, ':',  'Color', [0.85 0    0   ], 'LineWidth', 2.0, ...
    'DisplayName', 'Exp.');

xlabel('Normalized Width',  'FontName', 'Helvetica', 'FontSize', 14);
ylabel('Normalized Height', 'FontName', 'Helvetica', 'FontSize', 14);

xlim([-1 1]);  ylim([0 1.05]);
legend('Location', 'northeast', 'Box', 'off', 'FontSize', 13);
set(gca, 'FontName', 'Helvetica', 'FontSize', 13, 'LineWidth', 1.2, ...
    'TickDir', 'in', 'XTick', -1:0.5:1, 'YTick', 0:0.2:1);
title(sprintf('RMSE = %.3f  |  Max = %.3f', rmse, max_err), ...
    'FontSize', 11, 'FontWeight', 'normal');
hold off;

%% ── 7. FIGURE 2: Experiment photo with smooth boundary overlay ───────────
% y_e is already SG-smoothed; convert back to image-row coordinates.
y_row_e = nrows_e - y_e;

figure('Color', 'k', 'Units', 'centimeters', 'Position', [18 2 22 14]);
imshow(img_exp);
hold on;
plot(x_e, y_row_e, '-', 'Color', [1 0.2 0.2], 'LineWidth', 2.5);
title('Experiment — extracted smooth boundary', ...
    'Color', 'w', 'FontName', 'Helvetica', 'FontSize', 13);
hold off;

%% ── 8. FIGURE 3: Silhouette shape overlap ────────────────────────────────
% valid_e = ~isnan(h_exp_rs);
% valid_t = ~isnan(h_tar_rs);
% 
% w_e    = w_grid(valid_e);   h_e_pl = h_exp_rs(valid_e);
% w_t    = w_grid(valid_t);   h_t_pl = h_tar_rs(valid_t);
% 
% figure('Color', 'w', 'Units', 'centimeters', 'Position', [42 2 14 12]);
% hold on; box on;
% 
% % Target: grey filled silhouette
% fill([w_t, fliplr(w_t)], [h_t_pl, zeros(1, numel(w_t))], ...
%     [0.75 0.75 0.75], ...
%     'EdgeColor', [0 0 0], 'LineWidth', 2.0, 'DisplayName', 'Tar.');
% 
% % Experiment: red dashed outline on top
% plot(w_e, h_e_pl, '--', 'Color', [0.85 0 0], 'LineWidth', 2.0, ...
%     'DisplayName', 'Exp.');
% 
% xlabel('Normalized Width',  'FontName', 'Helvetica', 'FontSize', 14);
% ylabel('Normalized Height', 'FontName', 'Helvetica', 'FontSize', 14);
% legend('Location', 'northeast', 'Box', 'off', 'FontSize', 13);
% set(gca, 'FontName', 'Helvetica', 'FontSize', 13, 'LineWidth', 1.2, ...
%     'TickDir', 'in', 'XTick', -1:0.5:1, 'YTick', 0:0.2:1);
% title(sprintf('RMSE = %.3f  |  Max = %.3f', rmse, max_err), ...
%     'FontSize', 11, 'FontWeight', 'normal');
% xlim([-1 1]);  ylim([0 1.05]);
% hold off;
