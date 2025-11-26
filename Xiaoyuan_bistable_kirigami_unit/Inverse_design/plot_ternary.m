function plot_ternary(alpha1, alpha2, alpha3, index)
% PLOT_TERNARY  Standalone symmetric ternary plot (no toolbox required)
% Now includes symbolic radian ticks: 0, π/6, π/3, π/2, 2π/3, 5π/6, π
%
% alpha1, alpha2, alpha3 : internal angles (rad)
% index                  : scalar field (e.g. eps_bist or eta)
%                          Points with index = NaN will not be plotted.
%
% Axis convention (matching your figure):
%   - alpha3: bottom edge, 0 → π from left to right
%   - alpha2: right  edge, 0 → π from bottom to top
%   - alpha1: left   edge, 0 → π from top to bottom

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
    Afull = [Afull; Aang(:,P(i,:))];
    Vfull = [Vfull; V];
end

alpha1 = Afull(:,1);
alpha2 = Afull(:,2);
alpha3 = Afull(:,3);
index  = Vfull;

%% ---- 2) Convert to barycentric coordinates ----
% Mapping to vertices:
%   A ↔ alpha1  (left-bottom vertex)
%   B ↔ alpha3  (right-bottom vertex)
%   C ↔ alpha2  (top vertex)
S = alpha1 + alpha2 + alpha3;   % ≈ π
A = alpha1 ./ S;
B = alpha3 ./ S;
C = alpha2 ./ S;

[X, Y] = ternaryToCartesian(A,B,C);

%% ---- 3) Keep only finite index values (外面已经决定单/双稳) ----
mask = isfinite(index);
X = X(mask);
Y = Y(mask);
val = index(mask);

%% ---- 4) Plot frame, grid, axis ticks ----
figure; hold on; axis equal; axis off;
drawTernaryAxesWithSymbolicTicks();

%% ---- 5) Draw points ----
scatter(X, Y, 55, val, 'filled');

colormap(jet);

cb = colorbar('southoutside');   % put colorbar UNDER the triangle
cb.Label.Interpreter = 'none';
cb.Label.FontSize = 18;           % make label big
cb.FontSize = 18;                 % make tick labels bigger
pos = cb.Position;

pos(1) = 0.15;   % left edge
pos(3) = 0.45;   % width (0.15 → 0.60)
pos(2) = pos(2) - 0.05;   % move downward a bit

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
% A,B,C >= 0, A+B+C = 1
x = 0.5*(2*B + C);
y = (sqrt(3)/2)*C;
end


%% ------------------------------------------------------------------------
% DRAW TERNARY AXES + GRID + SYMBOLIC RADIAN TICKS
%% ------------------------------------------------------------------------
function drawTernaryAxesWithSymbolicTicks()

% Outer triangle
plot([0 1 0.5 0], [0 0 sqrt(3)/2 0], 'k-', 'LineWidth', 2);
hold on;

g   = linspace(0,1,7);   % 0, 1/6, ..., 1
symTicks = { ...
    '0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi'};

%% ---- Grid lines: constant A, B, C ----
for i = 2:6   % 跳过 0 和 1，避免画在边框上
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

%% ---- Symbolic tick labels: 0, π/6, ..., π ----
for i = 1:numel(g)
    r     = g(i);          % = angle / π
    label = symTicks{i};

    % --- α3 bottom edge: left → right ---
    % B = r, C = 0, A = 1-r
    [xb,yb] = ternaryToCartesian(1-r, r, 0);
    text(xb, yb - 0.05, label, 'FontSize', 18, ...
        'HorizontalAlignment','center','Interpreter','tex');

    % --- α1 left edge: top → bottom ---
    % A = r, C = 1-r, B = 0
    [x1,y1] = ternaryToCartesian(r, 0, 1-r);
    text(x1 - 0.05, y1, label, 'FontSize', 18, ...
        'HorizontalAlignment','center','Rotation', 60,...
        'Interpreter','tex');

    % --- α2 right edge: bottom → top ---
    % C = r, B = 1-r, A = 0
    [x2,y2] = ternaryToCartesian(0, 1-r, r);
    text(x2 + 0.05, y2, label, 'FontSize', 18, ...
        'HorizontalAlignment','center','Rotation', -60,...
        'Interpreter','tex');
end

end
