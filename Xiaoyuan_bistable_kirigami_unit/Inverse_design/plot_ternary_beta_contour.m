function plot_ternary_beta_contour(a1, a2, a3, val, beta, opts)
% PLOT_TERNARY_BETA_CONTOUR  3-D ternary + beta figure with stacked contours
%
% XY base  : ternary plot of (alpha1, alpha2, alpha3)
% Z axis   : tilting angle beta
% Contour  : filled patch + boundary line at each beta slice,
%            where val >= opts.threshold  defines the bistable region.
%
% Usage:
%   plot_ternary_beta_contour(a1, a2, a3, eps_bist, beta)
%   plot_ternary_beta_contour(a1, a2, a3, eps_bist, beta, opts)
%
% opts fields (all optional):
%   .threshold  - bistability cut-off (default 0)
%   .nGrid      - interpolation grid size inside triangle (default 100)
%   .faceAlpha  - patch transparency 0-1 (default 0.4)
%   .cmap       - colormap name for beta levels (default 'parula')
%   .lineWidth  - contour boundary line width (default 1.5)
%   .symmetrize - true = apply all 6 angle permutations (default true)

%% ---- 0. Default options -----------------------------------------------
if nargin < 6 || isempty(opts), opts = struct(); end
threshold  = getf(opts, 'threshold',  0);
nGrid      = getf(opts, 'nGrid',      100);
faceAlpha  = getf(opts, 'faceAlpha',  0.4);
cmapName   = getf(opts, 'cmap',       'parula');
lw         = getf(opts, 'lineWidth',  1.5);
do_sym     = getf(opts, 'symmetrize', true);

%% ---- 1. (Optional) Symmetrize all 6 permutations ----------------------
if do_sym
    [a1, a2, a3, val, beta] = symmetrize6(a1(:), a2(:), a3(:), val(:), beta(:));
end

%% ---- 2. Ternary -> Cartesian ------------------------------------------
S  = a1 + a2 + a3;
A  = a1 ./ S;   % -> left-bottom vertex (alpha1)
B  = a3 ./ S;   % -> right-bottom vertex (alpha3)
C  = a2 ./ S;   % -> top vertex (alpha2)
[Xp, Yp] = ternary2cart(A, B, C);

ok = isfinite(val) & isfinite(beta) & isfinite(Xp) & isfinite(Yp);
Xp = Xp(ok); Yp = Yp(ok); beta = beta(ok); val = val(ok);

%% ---- 3. Regular rectangular grid masked to ternary triangle -----------
xlin = linspace(0, 1,           nGrid);
ylin = linspace(0, sqrt(3)/2,  nGrid);
[Xg, Yg] = meshgrid(xlin, ylin);
in_tri    = inEquilateralTriangle(Xg, Yg);

%% ---- 4. Figure & ternary frame ----------------------------------------
figure; hold on; axis equal;
drawTernaryFrame3D();

%% ---- 5. Unique beta levels + colormap ---------------------------------
beta_unique = unique(beta);
nB          = numel(beta_unique);
cmap        = eval([cmapName '(nB)']);     % nB x 3

for ib = 1:nB
    b   = beta_unique(ib);
    sel = abs(beta - b) < 1e-12;
    xi  = Xp(sel);  yi = Yp(sel);  vi = val(sel);

    if numel(xi) < 6
        continue
    end

    % Interpolate scattered data onto the regular grid
    Vg         = griddata(xi, yi, vi, Xg, Yg, 'linear');
    % Points outside the data hull, or outside triangle, -> NaN
    Vg(~in_tri) = NaN;

    % Replace NaN with (threshold - 1) only for contourc (NaN breaks it)
    Vg_c          = Vg;
    Vg_c(isnan(Vg_c)) = threshold - 1;

    col = cmap(ib, :);

    % -- a) filled bistable patch at z = b --
    drawFilledPatch(xlin, ylin, Vg_c, b, threshold, col, faceAlpha);

    % -- b) boundary contour line at z = b --
    drawContourLine(xlin, ylin, Vg_c, b, threshold, col, lw);
end

%% ---- 6. Axes, colorbar ------------------------------------------------
xlim([0 1]); ylim([0 sqrt(3)/2]);
zlim([min(beta_unique) max(beta_unique)]);

% Beta colorbar using the actual beta range
colormap(cmap);
caxis([min(beta_unique) max(beta_unique)]);
cb = colorbar('eastoutside');
cb.Label.String    = '\beta  (rad)';
cb.Label.Interpreter = 'tex';
cb.Label.FontSize  = 14;
cb.FontSize        = 12;

zlabel('\beta  (rad)', 'FontSize', 14, 'Interpreter', 'tex');
title('Bistable region in ternary space vs \beta', ...
      'Interpreter', 'tex', 'FontSize', 15);

grid on;
view(40, 25);
end


%% ========================================================================
%  LOCAL HELPERS
%% ========================================================================

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

% ---- Symmetrize across all 6 angle permutations ------------------------
function [a1o, a2o, a3o, vo, bo] = symmetrize6(a1, a2, a3, v, b)
P = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];
A  = [a1, a2, a3];
a1o = []; a2o = []; a3o = []; vo = []; bo = [];
for k = 1:6
    a1o = [a1o; A(:, P(k,1))]; %#ok<AGROW>
    a2o = [a2o; A(:, P(k,2))]; %#ok<AGROW>
    a3o = [a3o; A(:, P(k,3))]; %#ok<AGROW>
    vo  = [vo;  v]; %#ok<AGROW>
    bo  = [bo;  b]; %#ok<AGROW>
end
end

% ---- Ternary barycentric -> Cartesian ----------------------------------
function [x, y] = ternary2cart(A, B, C)
x = 0.5 .* (2.*B + C);
y = (sqrt(3)/2) .* C;
end

% ---- Point-in-equilateral-triangle test --------------------------------
function inside = inEquilateralTriangle(X, Y)
% Vertices: (0,0), (1,0), (0.5, sqrt(3)/2)
% Uses the signed-area / half-plane method
h = sqrt(3)/2;
% Edge 1: (0,0) -> (1,0)     normal points up
d1 = Y >= 0;
% Edge 2: (1,0) -> (0.5,h)   normal points left-down
d2 = (-h).*(X - 1) - (-0.5).*(Y - 0) >= 0;
% Edge 3: (0.5,h) -> (0,0)   normal points right-down
d3 = (-h).*(X - 0.5) - (0.5).*(Y - h) <= 0;
inside = d1 & d2 & ~d3;   % careful with sign

% Simpler: use barycentric coords
% Convert (X,Y) -> (A,B,C) and check all >= 0
denom = h;   % height of triangle
C_bary = Y ./ denom;
B_bary = X - Y ./ sqrt(3);
A_bary = 1 - B_bary - C_bary;
inside = (A_bary >= -1e-9) & (B_bary >= -1e-9) & (C_bary >= -1e-9);
end

% ---- Draw filled bistable patch at a given beta height -----------------
function drawFilledPatch(xlin, ylin, Vg_c, b, thresh, col, alpha)
try
    C = contourc(xlin, ylin, Vg_c, [thresh thresh]);
catch
    return
end
i = 1;
while i < size(C, 2)
    n  = C(2, i);
    cx = C(1, i+1 : i+n);
    cy = C(2, i+1 : i+n);
    patch('XData', cx, 'YData', cy, 'ZData', b*ones(1, n), ...
          'FaceColor', col, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    i = i + n + 1;
end
end

% ---- Draw contour boundary line at a given beta height -----------------
function drawContourLine(xlin, ylin, Vg_c, b, thresh, col, lw)
try
    C = contourc(xlin, ylin, Vg_c, [thresh thresh]);
catch
    return
end
i = 1;
while i < size(C, 2)
    n  = C(2, i);
    cx = C(1, i+1 : i+n);
    cy = C(2, i+1 : i+n);
    plot3(cx, cy, b*ones(1, n), '-', 'Color', col, 'LineWidth', lw);
    i = i + n + 1;
end
end

% ---- Draw the 3-D ternary frame (outer triangle + grid + tick labels) --
function drawTernaryFrame3D()
h = sqrt(3)/2;

% Outer triangle at z = 0
plot3([0 1 0.5 0], [0 0 h 0], [0 0 0 0], 'k-', 'LineWidth', 2);

g        = linspace(0, 1, 7);
symTicks = {'0','\pi/6','\pi/3','\pi/2','2\pi/3','5\pi/6','\pi'};

% Grid lines (constant alpha1, alpha2, alpha3)
for i = 2:6
    v = g(i);
    col_g = [0.85 0.85 0.85];

    [x1,y1] = ternary2cart(v,   0,   1-v);
    [x2,y2] = ternary2cart(v, 1-v,   0  );
    plot3([x1 x2],[y1 y2],[0 0], 'Color', col_g);

    [x1,y1] = ternary2cart(0,   v,   1-v);
    [x2,y2] = ternary2cart(1-v, v,   0  );
    plot3([x1 x2],[y1 y2],[0 0], 'Color', col_g);

    [x1,y1] = ternary2cart(0,   1-v, v  );
    [x2,y2] = ternary2cart(1-v, 0,   v  );
    plot3([x1 x2],[y1 y2],[0 0], 'Color', col_g);
end

% Tick labels
for i = 1:numel(g)
    r  = g(i);
    lb = symTicks{i};

    % alpha3 bottom edge (left -> right)
    [xb, yb] = ternary2cart(1-r, r, 0);
    text(xb, yb - 0.05, 0, lb, 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'Interpreter', 'tex');

    % alpha1 left edge (top -> bottom)
    [x1, y1] = ternary2cart(r, 0, 1-r);
    text(x1 - 0.05, y1, 0, lb, 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'Rotation', 60, 'Interpreter', 'tex');

    % alpha2 right edge (bottom -> top)
    [x2, y2] = ternary2cart(0, 1-r, r);
    text(x2 + 0.05, y2, 0, lb, 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'Rotation', -60, 'Interpreter', 'tex');
end

% Axis labels
text(-0.10, h/2, 0, '\alpha_1', 'FontSize', 14, ...
    'Rotation', 60, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
text(1.10,  h/2, 0, '\alpha_2', 'FontSize', 14, ...
    'Rotation', -60, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
text(0.50,  -0.14, 0, '\alpha_3', 'FontSize', 14, ...
    'HorizontalAlignment', 'center', 'Interpreter', 'tex');
end
