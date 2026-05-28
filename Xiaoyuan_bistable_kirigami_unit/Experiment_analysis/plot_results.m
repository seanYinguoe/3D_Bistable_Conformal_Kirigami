% Plot tensile test results for kirigami tessellation
% Data: tessellation_stretch_7.csv
%   Columns : Time (s) | Displacement (mm) | Force (kN)
%   All recorded values are negative (machine convention).
%   Processing: negate both channels → positive extension & tensile force.
%               Force converted kN → N.
%               Displacement zero-offset to start at 0.

clc; clear; close all;

%% ── Read CSV ──────────────────────────────────────────────────────────────
% Row 1 = variable names, Row 2 = unit labels  →  data starts at row 3.
% Values are stored as quoted strings; MATLAB strips quotes automatically.
data_dir = fileparts(mfilename('fullpath'));
filename  = fullfile(data_dir, 'tessellation_stretch_7.csv');

opts = delimitedTextImportOptions( ...
    'NumVariables', 3, ...
    'DataLines',    [3, Inf], ...
    'Delimiter',    ',', ...
    'VariableNames', {'Time_s', 'Displacement_mm', 'Force_kN'}, ...
    'VariableTypes', {'double',  'double',          'double'});

T = readtable(filename, opts);

%% ── Process ───────────────────────────────────────────────────────────────
% Negate: machine pulled in negative direction; flip to positive convention.
displacement = -T.Displacement_mm;          % mm  (positive extension)
force        = -T.Force_kN * 1e3;           % N   (positive tensile load)

% Zero-offset displacement so the curve starts at (0, F0).
displacement = displacement - displacement(1);

% Smooth force signal with a Savitzky-Golay filter.
% Window of 201 pts (~2 mm at the data density) removes high-freq noise
% while preserving the overall curve shape.
poly_order  = 3;
win_size    = 201;   % must be odd; increase to smooth more aggressively
force_smooth = sgolayfilt(force, poly_order, win_size);

%% ── Plot ──────────────────────────────────────────────────────────────────
figure('Color', 'w', 'Units', 'centimeters', 'Position', [5 5 14 10]);
hold on; box on;

plot(displacement, force_smooth, '-', ...
    'Color',     [0.85 0.33 0.10], ...
    'LineWidth', 1.5);

xlabel('Displacement (mm)', 'Interpreter', 'tex', 'FontSize', 20);
ylabel('Force (N)',          'Interpreter', 'tex', 'FontSize', 20);
title('Tessellation Stretch – Sample 7', ...
    'FontSize', 18, 'FontWeight', 'normal');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 18, 'LineWidth', 1.5);

xlim([0, 60]);
ylim([0, 20]);

grid off;
hold off;

fprintf('Displacement range : %.2f – %.2f mm\n', ...
    min(displacement), max(displacement));
fprintf('Force range        : %.2f – %.2f N\n', ...
    min(force), max(force));
