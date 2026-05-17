function plot_ternary(alpha1, alpha2, alpha3, index, mode, clim, alphaRange)
% PLOT_TERNARY  Standalone symmetric ternary plot (no toolbox required)
%
% alpha1, alpha2, alpha3 : internal angles (rad)
% index                  : scalar field (e.g. eps_bist or eta)
%                          Points with index = NaN will not be plotted.
% mode       (optional)  : 'scatter' (default) | 'field'
% clim       (optional)  : color limits [cmin cmax]
% alphaRange (optional)  : [alpha_lo alpha_hi] in degrees, default [0 180]
%                          Zooms the ternary plot to this sub-range.
%
% Axis convention:
%   - alpha3: bottom edge
%   - alpha2: right  edge
%   - alpha1: left   edge

if nargin < 5 || isempty(mode),       mode = 'scatter'; end
if nargin < 6,                         clim = [];        end
if nargin < 7 || isempty(alphaRange), alphaRange = [0 180]; end

alpha_lo_deg  = alphaRange(1);
alpha_hi_deg  = alphaRange(2);
alpha_lo_rad  = alpha_lo_deg * pi / 180;
alpha_denom_rad = (180 - 3*alpha_lo_deg) * pi / 180;   % pi - 3*alpha_lo_rad

%% ---- 1) Symmetrize points across all 6 permutations (full triangle) ----
Aang = [alpha1(:), alpha2(:), alpha3(:)];
V    = index(:);

P = [
    1 2 3;
    1 3 2;
    2 1 3;
    2 3 1;
    3 1 2;
    3 2 1
];

Afull = [];
Vfull = [];

for i = 1:6
    Afull = [Afull; Aang(:,P(i,:))]; %#ok<AGROW>
    Vfull = [Vfull; V];              %#ok<AGROW>
end

alpha1 = Afull(:,1);
alpha2 = Afull(:,2);
alpha3 = Afull(:,3);
index  = Vfull;

%% ---- 2) Rescale to zoomed barycentric coordinates ----
% A' = (alpha_i - alpha_lo) / (pi - 3*alpha_lo)  →  sum = 1 exactly
Ap = (alpha1 - alpha_lo_rad) ./ alpha_denom_rad;
Bp = (alpha3 - alpha_lo_rad) ./ alpha_denom_rad;
Cp = (alpha2 - alpha_lo_rad) ./ alpha_denom_rad;

[X, Y] = ternaryToCartesian(Ap, Bp, Cp);

%% ---- 3) Keep only finite index values and points inside zoomed range ----
inRange = (Ap >= -1e-6) & (Bp >= -1e-6) & (Cp >= -1e-6) & ...
          (Ap <=  1+1e-6) & (Bp <= 1+1e-6) & (Cp <= 1+1e-6);
mask = isfinite(index) & isfinite(X) & isfinite(Y) & inRange;
X   = X(mask);
Y   = Y(mask);
val = index(mask);

%% ---- 4) Plot frame, grid, axis ticks ----
figure; hold on; axis equal; axis off;
drawTernaryAxesWithZoomedTicks(alpha_lo_deg, alpha_hi_deg, alpha_denom_rad);

%% ---- 5) Draw data ----
switch lower(mode)
    case 'field'
        if numel(X) >= 3
            tri = delaunay(X, Y);
            trisurf(tri, X, Y, val, 'EdgeColor', 'none');
            view(2);
        end
    otherwise
        scatter(X, Y, 55, val, 'filled');
end

colormap(jet);
if ~isempty(clim) && numel(clim) == 2 && all(isfinite(clim))
    caxis(clim);
end

cb = colorbar('southoutside');
cb.Label.Interpreter = 'none';
cb.Label.FontSize = 18;
cb.FontSize = 18;
pos = cb.Position;
pos(1) = 0.15;
pos(3) = 0.45;
pos(2) = pos(2) - 0.05;
cb.Position = pos;

%% ---- 6) Axis labels ----
text(0.05,  sqrt(3)/4, '\alpha_1', 'FontSize', 18, ...
    'Rotation', 60, 'HorizontalAlignment','center');
text(0.95,  sqrt(3)/4, '\alpha_2', 'FontSize', 18, ...
    'Rotation', -60, 'HorizontalAlignment','center');
text(0.50, -0.12,       '\alpha_3', 'FontSize', 18, ...
    'HorizontalAlignment','center');

title('Ternary map of (\alpha_1,\alpha_2,\alpha_3)', ...
    'Interpreter','tex','FontSize',18);

end


%% ------------------------------------------------------------------------
% BARYCENTRIC → CARTESIAN
%% ------------------------------------------------------------------------
function [x, y] = ternaryToCartesian(A, B, C)
x = 0.5*(2*B + C);
y = (sqrt(3)/2)*C;
end


%% ------------------------------------------------------------------------
% DRAW TERNARY AXES + GRID + DEGREE TICKS (zoomed range)
%% ------------------------------------------------------------------------
function drawTernaryAxesWithZoomedTicks(alpha_lo_deg, alpha_hi_deg, alpha_denom_rad)

% Outer triangle
plot([0 1 0.5 0], [0 0 sqrt(3)/2 0], 'k-', 'LineWidth', 2);
hold on;

% Tick positions: every 5 degrees within [alpha_lo, alpha_hi]
tick_vals_deg = alpha_lo_deg : 5 : alpha_hi_deg;
alpha_denom_deg = (180 - 3*alpha_lo_deg);   % same denominator in degrees
g = (tick_vals_deg - alpha_lo_deg) / alpha_denom_deg;  % normalised [0,1]

% Grid lines at interior ticks (skip endpoints 0 and 1)
for i = 2 : numel(g)-1
    v = g(i);

    % constant A (alpha1)
    [x1,y1] = ternaryToCartesian(v,    0, 1-v);
    [x2,y2] = ternaryToCartesian(v, 1-v,   0);
    plot([x1 x2],[y1 y2],'Color',[0.9 0.9 0.9]);

    % constant B (alpha3)
    [x1,y1] = ternaryToCartesian(0,   v, 1-v);
    [x2,y2] = ternaryToCartesian(1-v, v,   0);
    plot([x1 x2],[y1 y2],'Color',[0.9 0.9 0.9]);

    % constant C (alpha2)
    [x1,y1] = ternaryToCartesian(0, 1-v,   v);
    [x2,y2] = ternaryToCartesian(1-v, 0,   v);
    plot([x1 x2],[y1 y2],'Color',[0.9 0.9 0.9]);
end

% Outward tick normals (unit vectors pointing away from triangle centre)
n_bottom = [0, -1];
n_left   = [-sqrt(3)/2,  0.5];
n_right  = [ sqrt(3)/2,  0.5];
tick_len = 0.025;
label_gap = 0.06;

% Tick marks and labels
for i = 1:numel(g)
    r     = g(i);
    label = sprintf('%g\\circ', tick_vals_deg(i));

    % --- α3 bottom edge: left → right  (B=r, C=0, A=1-r) ---
    [xb,yb] = ternaryToCartesian(1-r, r, 0);
    plot([xb, xb + tick_len*n_bottom(1)], [yb, yb + tick_len*n_bottom(2)], 'k-', 'LineWidth', 1);
    text(xb + label_gap*n_bottom(1), yb + label_gap*n_bottom(2), label, ...
        'FontSize', 13, 'HorizontalAlignment','center','Interpreter','tex');

    % --- α1 left edge: top → bottom  (A=r, C=1-r, B=0) ---
    [x1,y1] = ternaryToCartesian(r, 0, 1-r);
    plot([x1, x1 + tick_len*n_left(1)], [y1, y1 + tick_len*n_left(2)], 'k-', 'LineWidth', 1);
    text(x1 + label_gap*n_left(1), y1 + label_gap*n_left(2), label, ...
        'FontSize', 13, 'HorizontalAlignment','center','Rotation', 60, 'Interpreter','tex');

    % --- α2 right edge: bottom → top  (C=r, B=1-r, A=0) ---
    [x2,y2] = ternaryToCartesian(0, 1-r, r);
    plot([x2, x2 + tick_len*n_right(1)], [y2, y2 + tick_len*n_right(2)], 'k-', 'LineWidth', 1);
    text(x2 + label_gap*n_right(1), y2 + label_gap*n_right(2), label, ...
        'FontSize', 13, 'HorizontalAlignment','center','Rotation', -60, 'Interpreter','tex');
end

end
