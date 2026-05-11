% compare_experiment_target.m
% Extract upper silhouette from experiment photo and rendered target image,
% normalise both to [-1,1] x [0,1], and plot the comparison.
%
% Experiment : white kirigami dome on black background
%   → binarise (Otsu), morphological closing to fill gaps between pieces,
%     convex hull to get clean outer silhouette
% Target     : gray hemisphere on white background
%   → threshold near white to capture full dome, fill holes

clc; clear; close all;

%% ── 1. SETTINGS ─────────────────────────────────────────────────────────
model = 'quarterdome';
exp_path    = strcat('Figures/',model,'_exp.jpg');
tar_path    = strcat('Figures/',model,'_tar.jpg');

n_compare   = 300;   % uniform width samples for comparison

% Morphological closing disk radius (pixels) for experiment image.
% Should be larger than the widest gap between kirigami pieces.
close_radius = 60;

%% ── 2. EXPERIMENT IMAGE ──────────────────────────────────────────────────
% White kirigami structure on dark background.
% Strategy: binarise → close gaps → convex hull → upper boundary.

img_exp  = imread(exp_path);
gray_exp = rgb2gray(img_exp);
bw_exp   = imbinarize(gray_exp);              % white=1, black=0

% Close gaps between kirigami pieces to merge them into one solid shape
se          = strel('disk', close_radius);
bw_closed   = imclose(bw_exp, se);

% Keep only the largest blob (the dome)
bw_dome     = bwareafilt(bw_closed, 1);

% Convex hull → smooth, gap-free silhouette
bw_hull     = bwconvhull(bw_dome);

[nrows_e, ~] = size(bw_hull);

% Upper boundary: topmost foreground pixel per column
cols_e = find(any(bw_hull, 1));
x_e    = zeros(1, numel(cols_e));
y_e    = zeros(1, numel(cols_e));
for k  = 1:numel(cols_e)
    rows   = find(bw_hull(:, cols_e(k)));
    x_e(k) = cols_e(k);
    y_e(k) = nrows_e - min(rows);         % flip: row 1 = top of image
end

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

%% ── 5. ERROR METRICS ─────────────────────────────────────────────────────

valid   = ~isnan(h_exp_rs) & ~isnan(h_tar_rs);
rmse    = sqrt(mean((h_exp_rs(valid) - h_tar_rs(valid)).^2));
max_err = max(abs(h_exp_rs(valid) - h_tar_rs(valid)));

fprintf('RMSE      = %.4f\n', rmse);
fprintf('Max error = %.4f\n', max_err);

%% ── 6. COMPARISON PLOT ───────────────────────────────────────────────────

figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 14 12]);
hold on; box on;

plot(w_grid, h_tar_rs, '-',  'Color', [0    0    0   ], 'LineWidth', 2.0, ...
    'DisplayName', 'Tar.');
plot(w_grid, h_exp_rs, ':',  'Color', [0.85 0    0   ], 'LineWidth', 2.0, ...
    'DisplayName', 'Exp.');

xlabel('Normalized Width',  'FontName', 'Helvetica', 'FontSize', 14);
ylabel('Normalized Height', 'FontName', 'Helvetica', 'FontSize', 14);

xlim([-1 1]);
ylim([ 0 1.05]);

legend('Location', 'northeast', 'Box', 'off', 'FontSize', 13);
set(gca, 'FontName', 'Helvetica', 'FontSize', 13, 'LineWidth', 1.2, ...
    'TickDir', 'in', 'XTick', -1:0.5:1, 'YTick', 0:0.2:1);
title(sprintf('RMSE = %.3f  |  Max = %.3f', rmse, max_err), ...
    'FontSize', 11, 'FontWeight', 'normal');

hold off;

