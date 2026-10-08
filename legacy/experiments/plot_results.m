% Plot tensile test results for kirigami tessellation
% Data: tessellation_stretch_7.csv, tessellation_stretch_11.csv
%   Columns : Time (s) | Displacement (mm) | Force (kN)
%   All recorded values are negative (machine convention).
%   Processing: negate both channels → positive extension & tensile force.
%               Force converted kN → N.
%               Displacement zero-offset to start at 0.

clc; clear; close all;

%% ── Helper: load and process one CSV ─────────────────────────────────────
    function [displacement, force_smooth] = load_sample(filepath)
        opts = delimitedTextImportOptions( ...
            'NumVariables', 3, ...
            'DataLines',    [3, Inf], ...
            'Delimiter',    ',', ...
            'VariableNames', {'Time_s', 'Displacement_mm', 'Force_kN'}, ...
            'VariableTypes', {'double',  'double',          'double'});
        T = readtable(filepath, opts);

        displacement = -T.Displacement_mm;
        force        = -T.Force_kN * 1e3;
        displacement = displacement - displacement(1);

        poly_order   = 3;
        win_size     = 501;   % larger window → stronger noise suppression
        force_smooth = sgolayfilt(force, poly_order, win_size);
    end

%% ── Read & process both samples ──────────────────────────────────────────
data_dir = fileparts(mfilename('fullpath'));

[disp7,  force7]  = load_sample(fullfile(data_dir, 'tessellation_stretch_8.csv'));
[disp11, force11] = load_sample(fullfile(data_dir, 'tessellation_stretch_11.csv'));

%% ── Plot ──────────────────────────────────────────────────────────────────
figure('Color', 'w', 'Units', 'centimeters', 'Position', [5 5 14 10]);
hold on; box on;

plot(disp7,  force7,  '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
plot(disp11, force11, '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5);

xlabel('Displacement (mm)', 'Interpreter', 'tex', 'FontSize', 20);
ylabel('Force (N)',          'Interpreter', 'tex', 'FontSize', 20);
title('Tessellation Stretch', 'FontSize', 18, 'FontWeight', 'normal');

legend('Sample 7', 'Sample 11', ...
    'Location', 'northwest', 'FontSize', 16, 'Box', 'off');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 18, 'LineWidth', 1.5);

xlim([0, 60]);
ylim([0, 20]);

grid off;
hold off;

fprintf('Sample 7  – Displacement: %.2f – %.2f mm | Force: %.2f – %.2f N\n', ...
    min(disp7),  max(disp7),  min(force7),  max(force7));
fprintf('Sample 11 – Displacement: %.2f – %.2f mm | Force: %.2f – %.2f N\n', ...
    min(disp11), max(disp11), min(force11), max(force11));
