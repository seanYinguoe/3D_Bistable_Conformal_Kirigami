%% Visualize ternary study results
% Load a saved output folder from ternary_study.m and produce ternary plots.
%
% Usage:
%   1. Set out_root to the desired output folder (e.g. 'output/ternary_refine_...')
%   2. Run this script to generate 2-D ternary, 3-D ternary, and beta-contour figures.

clear; clc;

%% ---- Select output folder ----
out_root = fullfile('output', 'ternary_refine_20260430_123351');%

% ---- Load saved data ----
S = load(fullfile(out_root, 'anisotropy_filter_refine.mat'));
data = S.anisotropy_filter_refine;
% data = anisotropy_filter;

%% ---- Parameters for a single beta slice ----
etaThreshold = 0.15;
beta_plot    = 0;        % pick a beta value present in the data
beta_tol     = 1e-4;

%% ---- 2D ternary: plot_ternary ----
T0   = data(abs(data.beta - beta_plot) < beta_tol, :);
isB  = isfinite(T0.eps_bist) & isfinite(T0.eta_val) & (T0.eta_val > etaThreshold);
T_bist = T0(isB, :);

alphaRange2D = [50, 80];   % degrees shown on 2D ternary axes

if ~isempty(T_bist)
    % eps_bist map
    plot_ternary(T_bist.a1, T_bist.a2, T_bist.a3, T_bist.eps_bist, 'scatter', [], alphaRange2D);
    title(sprintf('\\epsilon_{bist} ternary map  (\\beta = %.4f)', beta_plot), ...
        'Interpreter', 'tex', 'FontSize', 18);

    % eta map
    plot_ternary(T_bist.a1, T_bist.a2, T_bist.a3, T_bist.eta_val, 'scatter', [], alphaRange2D);
    title(sprintf('\\eta ternary map  (\\beta = %.4f)', beta_plot), ...
        'Interpreter', 'tex', 'FontSize', 18);
else
    warning('No bistable points at beta = %.4f (tol = %.1e). Check beta_plot.', beta_plot, beta_tol);
end

%% ---- 3D ternary: plot_ternary_3D ----
% plot_ternary_3D(data.a1, data.a2, data.a3, ...
%     data.eps_bist, data.beta);

%% ---- Beta-contour surface: plot_ternary_beta_contour ----
opts_contour = struct('threshold', 0, 'etaThreshold', 0.20, 'nGrid', 120, 'faceAlpha', 0.55, 'alphaRange', [50, 70]);
plot_ternary_beta_contour(data.a1, data.a2, data.a3, ...
    data.eps_bist, data.beta, data.eta_val, opts_contour);
