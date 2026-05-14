%% Visualize ternary study results
% Load a saved output folder from ternary_study.m and produce ternary plots.
%
% Usage:
%   1. Set out_root to the desired output folder (e.g. 'output/ternary_refine_...')
%   2. Run this script to generate 2-D ternary, 3-D ternary, and beta-contour figures.

clear; clc;

%% ---- Select output folder ----
out_root = fullfile('output', 'ternary_refine_20260430_123351');   % <-- change as needed

%% ---- Load saved data ----
S = load(fullfile(out_root, 'anisotropy_filled_refine.mat'));
anisotropy_filled = S.anisotropy_filled_refine;
anisotropy_clean  = S.anisotropy_clean;

%% ---- Filter for finite values ----
mask = isfinite(anisotropy_filled.eta_val) & isfinite(anisotropy_filled.eps_bist);
anisotropy_filter = anisotropy_filled(mask, :);

%% ---- Parameters for a single beta slice ----
etaThreshold = 0.03;
beta_plot    = 0;        % pick a beta value present in the data
beta_tol     = 1e-4;

%% ---- 2D ternary: plot_ternary ----
T0   = anisotropy_filter(abs(anisotropy_filter.beta - beta_plot) < beta_tol, :);
isB  = isfinite(T0.eps_bist) & isfinite(T0.eta_val) & (T0.eta_val > etaThreshold);
T_bist = T0(isB, :);

if ~isempty(T_bist)
    % eps_bist map
    plot_ternary(T_bist.a1, T_bist.a2, T_bist.a3, T_bist.eps_bist, 'scatter');
    title(sprintf('\\epsilon_{bist} ternary map  (\\beta = %.4f)', beta_plot), ...
        'Interpreter', 'tex', 'FontSize', 18);

    % eta map
    plot_ternary(T_bist.a1, T_bist.a2, T_bist.a3, T_bist.eta_val, 'scatter');
    title(sprintf('\\eta ternary map  (\\beta = %.4f)', beta_plot), ...
        'Interpreter', 'tex', 'FontSize', 18);
else
    warning('No bistable points at beta = %.4f (tol = %.1e). Check beta_plot.', beta_plot, beta_tol);
end

%% ---- 3D ternary: plot_ternary_3D ----
plot_ternary_3D(anisotropy_clean.a1, anisotropy_clean.a2, anisotropy_clean.a3, ...
    anisotropy_clean.eps_bist, anisotropy_clean.beta);

%% ---- Beta-contour stack: plot_ternary_beta_contour ----
opts_contour = struct('threshold', 0, 'nGrid', 100, 'faceAlpha', 0.4, 'cmap', 'parula');
plot_ternary_beta_contour(anisotropy_clean.a1, anisotropy_clean.a2, anisotropy_clean.a3, ...
    anisotropy_clean.eps_bist, anisotropy_clean.beta, opts_contour);
