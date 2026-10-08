function plot_ternary(alpha1, alpha2, alpha3, index, mode, clim, alphaRange)
% PLOT_TERNARY  Standalone symmetric ternary plot (no toolbox required)
%
% alpha1, alpha2, alpha3 : internal angles (rad)
% index                  : scalar field (e.g. eps_bist or eta)
% mode       (optional)  : 'scatter' (default) | 'field'
% clim       (optional)  : color limits [cmin cmax]
% alphaRange (optional)  : [alpha_lo alpha_hi] in degrees, default [0 180]
%
% Font conventions:
%   Greek symbols : LaTeX math italic (Times New Roman style)
%   All other text: Arial

if nargin < 5 || isempty(mode),       mode = 'scatter'; end
if nargin < 6,                         clim = [];        end
if nargin < 7 || isempty(alphaRange), alphaRange = [0 180]; end

%% ---- Shared style constants ----
fsize   = 16;
fnArial = 'Arial';

alpha_lo_deg    = alphaRange(1);
alpha_hi_deg    = alphaRange(2);
alpha_lo_rad    = alpha_lo_deg * pi / 180;
alpha_denom_rad = (180 - 3*alpha_lo_deg) * pi / 180;

%% ---- 1) Symmetrize across all 6 permutations ----
Aang = [alpha1(:), alpha2(:), alpha3(:)];
V    = index(:);
P = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];
Afull = []; Vfull = [];
for i = 1:6
    Afull = [Afull; Aang(:,P(i,:))]; %#ok<AGROW>
    Vfull = [Vfull; V];              %#ok<AGROW>
end
alpha1 = Afull(:,1); alpha2 = Afull(:,2); alpha3 = Afull(:,3);
index  = Vfull;

%% ---- 2) Rescale to zoomed barycentric coordinates ----
Ap = (alpha1 - alpha_lo_rad) ./ alpha_denom_rad;
Bp = (alpha3 - alpha_lo_rad) ./ alpha_denom_rad;
Cp = (alpha2 - alpha_lo_rad) ./ alpha_denom_rad;
[X, Y] = ternaryToCartesian(Ap, Bp, Cp);

%% ---- 3) Filter to valid points inside range ----
inRange = (Ap >= -1e-6) & (Bp >= -1e-6) & (Cp >= -1e-6) & ...
          (Ap <= 1+1e-6) & (Bp <= 1+1e-6) & (Cp <= 1+1e-6);
mask = isfinite(index) & isfinite(X) & isfinite(Y) & inRange;
X = X(mask); Y = Y(mask); val = index(mask);

%% ---- 4) Figure layout ----
% Fixed figure size with explicit margins so title, triangle, and
% colorbar never overlap.
fig = figure('Units', 'centimeters', 'Position', [2 2 13 15]);
set(fig, 'DefaultTextFontName', fnArial, 'DefaultTextFontSize', fsize);

% Axes occupies the middle band: 18 % top margin (title), 25 % bottom (colorbar + α3 label)
ax = axes('Parent', fig, 'Position', [0.10 0.24 0.80 0.60]);
hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');

drawTernaryAxesWithZoomedTicks(alpha_lo_deg, alpha_hi_deg, fsize, fnArial, ax);

%% ---- 5) Data ----
switch lower(mode)
    case 'field'
        if numel(X) >= 3
            tri = delaunay(X, Y);
            trisurf(tri, X, Y, val, 'EdgeColor', 'none', 'Parent', ax);
            view(ax, 2);
        end
    otherwise
        scatter(ax, X, Y, 60, val, 'filled');
end

colormap(ax, jet);
if ~isempty(clim) && numel(clim) == 2 && all(isfinite(clim))
    caxis(ax, clim);
end

%% ---- 6) Colorbar — fixed position, 5 ticks, 2 sig-figs ----
cb = colorbar(ax, 'southoutside');
cb.FontName  = fnArial;
cb.FontSize  = fsize;
cb.Label.FontName    = fnArial;
cb.Label.FontSize    = fsize;
cb.Label.Interpreter = 'none';

clim_cur      = caxis(ax);
cb.Ticks      = linspace(clim_cur(1), clim_cur(2), 5);
cb.TickLabels = arrayfun(@(v) sprintf('%.2g', v), cb.Ticks, 'UniformOutput', false);

% Place colorbar explicitly at the bottom of the figure
drawnow;
ax_pos = ax.Position;          % [left bottom width height] in normalised units
cb_w   = ax_pos(3) * 0.72;
cb_l   = ax_pos(1) + (ax_pos(3) - cb_w) / 2;
cb.Position = [cb_l, 0.07, cb_w, 0.030];

%% ---- 7) Axis labels (LaTeX italic for Greek) ----
text(ax, 0.05, sqrt(3)/4, '$\alpha_1$', ...
    'Interpreter', 'latex', 'FontSize', fsize+2, ...
    'Rotation', 60, 'HorizontalAlignment', 'center');
text(ax, 0.95, sqrt(3)/4, '$\alpha_2$', ...
    'Interpreter', 'latex', 'FontSize', fsize+2, ...
    'Rotation', -60, 'HorizontalAlignment', 'center');
text(ax, 0.50, -0.20, '$\alpha_3$', ...
    'Interpreter', 'latex', 'FontSize', fsize+2, ...
    'HorizontalAlignment', 'center');

end


%% ------------------------------------------------------------------------
function [x, y] = ternaryToCartesian(A, B, C)
x = 0.5*(2*B + C);
y = (sqrt(3)/2)*C;
end


%% ------------------------------------------------------------------------
function drawTernaryAxesWithZoomedTicks(alpha_lo_deg, alpha_hi_deg, fsize, fnArial, ax)

% Outer triangle
plot(ax, [0 1 0.5 0], [0 0 sqrt(3)/2 0], 'k-', 'LineWidth', 2);

tick_vals_deg   = alpha_lo_deg : 10 : alpha_hi_deg;
alpha_denom_deg = 180 - 3*alpha_lo_deg;
g = (tick_vals_deg - alpha_lo_deg) / alpha_denom_deg;

% Grid lines (interior ticks only)
for i = 2 : numel(g)-1
    v = g(i);
    [x1,y1] = ternaryToCartesian(v,    0, 1-v);
    [x2,y2] = ternaryToCartesian(v, 1-v,   0);
    plot(ax, [x1 x2],[y1 y2],'Color',[0.88 0.88 0.88]);

    [x1,y1] = ternaryToCartesian(0,   v, 1-v);
    [x2,y2] = ternaryToCartesian(1-v, v,   0);
    plot(ax, [x1 x2],[y1 y2],'Color',[0.88 0.88 0.88]);

    [x1,y1] = ternaryToCartesian(0, 1-v,   v);
    [x2,y2] = ternaryToCartesian(1-v, 0,   v);
    plot(ax, [x1 x2],[y1 y2],'Color',[0.88 0.88 0.88]);
end

% Outward tick normals
n_bottom = [0, -1];
n_left   = [-sqrt(3)/2,  0.5];
n_right  = [ sqrt(3)/2,  0.5];
tick_len  = 0.025;
label_gap = 0.085;   % larger gap to separate labels from triangle edges

for i = 1:numel(g)
    r     = g(i);
    label = sprintf('%g^\\circ', tick_vals_deg(i));

    % α3 bottom edge (skip endpoints shared with left/right edge labels)
    [xb,yb] = ternaryToCartesian(1-r, r, 0);
    plot(ax, [xb, xb + tick_len*n_bottom(1)], [yb, yb + tick_len*n_bottom(2)], 'k-', 'LineWidth', 1);
    text(ax, xb + label_gap*n_bottom(1), yb + label_gap*n_bottom(2), label, ...
        'FontSize', fsize, 'FontName', fnArial, ...
        'HorizontalAlignment', 'center', 'Interpreter', 'tex');

    % α1 left edge
    [x1,y1] = ternaryToCartesian(r, 0, 1-r);
    plot(ax, [x1, x1 + tick_len*n_left(1)], [y1, y1 + tick_len*n_left(2)], 'k-', 'LineWidth', 1);
    text(ax, x1 + label_gap*n_left(1), y1 + label_gap*n_left(2), label, ...
        'FontSize', fsize, 'FontName', fnArial, ...
        'HorizontalAlignment', 'center', 'Rotation', 60, 'Interpreter', 'tex');

    % α2 right edge
    [x2,y2] = ternaryToCartesian(0, 1-r, r);
    plot(ax, [x2, x2 + tick_len*n_right(1)], [y2, y2 + tick_len*n_right(2)], 'k-', 'LineWidth', 1);
    text(ax, x2 + label_gap*n_right(1), y2 + label_gap*n_right(2), label, ...
        'FontSize', fsize, 'FontName', fnArial, ...
        'HorizontalAlignment', 'center', 'Rotation', -60, 'Interpreter', 'tex');
end

end
