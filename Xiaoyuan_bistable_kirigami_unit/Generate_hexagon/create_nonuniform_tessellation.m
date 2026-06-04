function tessellation = create_nonuniform_tessellation(num_x, num_y, edgeLen, l1, l4, t)
%CREATE_NONUNIFORM_TESSELLATION  Triangle tessellation with beta varying along x.
%
%   tessellation = create_nonuniform_tessellation(num_x, num_y, edgeLen, l1, l4, t)
%
% Beta profile (cosine, symmetric about the center column):
%   beta = 0       at the center column
%   beta = pi/20   at the left and right edges
%
% The tessellation is plotted and saved as an SVG cut pattern.

if nargin < 6
    error('create_nonuniform_tessellation requires num_x, num_y, edgeLen, l1, l4, t.');
end

% Add Conformal_mapping folder so generate_svg is accessible
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, '..', 'Conformal_mapping'));

% Build per-face beta map varying along x.
% Faces are ordered row-by-row; within each row the x-column index is
% mod(face_index-1, num_x), ranging 0 .. num_x-1.
beta_max = pi / 25;
n_faces  = num_x * num_y;
x_cols   = mod((0:n_faces-1)', num_x);  % [0 .. num_x-1] for each face

% |cos| profile: 1 at edges (x_col=0 and x_col=num_x-1), 0 at center
beta_map = beta_max * abs(cos(pi * x_cols / (num_x - 1)));

% Generate tessellation using the existing pipeline
tessellation = triangle_tessellation(num_x, num_y, beta_map, edgeLen, l1, l4, t);

% Colour scheme: void / flank / filament / inner triangle
colour = {'white', ...
          [0.9216  0.8863  0.4235], ...
          [0.7059  0.9608  0.4118], ...
          [0.9216  0.8863  0.4235]};

% Plot
figure; hold on;
for i = 1:numel(tessellation)
    plot_triangle(tessellation{i}, colour);
end
hold off;
axis equal; axis off;
title(sprintf('Non-uniform \\beta: 0 (center) \\to \\pi/20 (edges),  %d\\times%d', num_x, num_y));

% Save SVG
svg_name = sprintf('nonuniform_beta_%dx%d.svg', num_x, num_y);
svg_path = generate_svg(tessellation, svg_name, false, [], 0.20);
fprintf('SVG saved to: %s\n', svg_path);

end
