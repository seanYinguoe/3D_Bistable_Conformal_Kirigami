function plot_ternary_3D(alpha1, alpha2, alpha3, index, beta)
% PLOT_TERNARY_3D  3D ternary plot with beta as Z-axis
%
% alpha1, alpha2, alpha3 : internal angles (rad)
% index                  : scalar field for colour (e.g. eps_bist)
% beta                   : tilting angle, used as Z coordinate
%
% Axis convention same as plot_ternary:
%   - alpha3: bottom edge, 0 → π from left to right
%   - alpha2: right  edge, 0 → π from bottom to top
%   - alpha1: left   edge, 0 → π from top to bottom

%% ---- 1) Symmetrize points across all 6 permutations (full triangle) ----
Aang = [alpha1(:), alpha2(:), alpha3(:)];
V    = index(:);
Bz   = beta(:);

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
Bfull = [];

for i = 1:6
    Afull = [Afull; Aang(:,P(i,:))];
    Vfull = [Vfull; V];
    Bfull = [Bfull; Bz];
end

alpha1 = Afull(:,1);
alpha2 = Afull(:,2);
alpha3 = Afull(:,3);
index  = Vfull;
beta   = Bfull;

%% ---- 2) Convert to barycentric → Cartesian (same as your 2D plot) ----
S = alpha1 + alpha2 + alpha3;   % ≈ π
A = alpha1 ./ S;
B = alpha3 ./ S;
C = alpha2 ./ S;

[X, Y] = ternaryToCartesian(A,B,C);

%% ---- 3) Keep only finite values ----
mask = isfinite(index) & isfinite(beta);
X    = X(mask);
Y    = Y(mask);
Z    = beta(mask);
val  = index(mask);

%% ---- 4) Draw ternary frame on z = 0 ----
figure; hold on; axis equal;
drawTernaryAxesWithSymbolicTicks3D();  % 3D version

%% ---- 5) Draw 3D points ----
scatter3(X, Y, Z, 55, val, 'filled');

colormap(jet);
cb = colorbar('eastoutside');
cb.Label.String = 'index (e.g. \epsilon_{bist})';
cb.Label.Interpreter = 'tex';
cb.Label.FontSize = 14;
cb.FontSize = 12;

xlabel('x (ternary)');
ylabel('y (ternary)');
zlabel('\beta');

title('3D Ternary map with \beta as Z-axis', ...
    'Interpreter','tex','FontSize',16);

grid on;
view(40, 25);
end


%% ------------------------------------------------------------------------
% BARYCENTRIC → CARTESIAN (copy from your existing file)
%% ------------------------------------------------------------------------
function [x, y] = ternaryToCartesian(A, B, C)
% A,B,C >= 0, A+B+C = 1
x = 0.5*(2*B + C);
y = (sqrt(3)/2)*C;
end


%% ------------------------------------------------------------------------
% 3D version of DRAW TERNARY AXES + GRID + SYMBOLIC RADIAN TICKS
%% ------------------------------------------------------------------------
function drawTernaryAxesWithSymbolicTicks3D()

% Outer triangle on z = 0
plot3([0 1 0.5 0], [0 0 sqrt(3)/2 0], [0 0 0 0], 'k-', 'LineWidth', 2);
hold on;

g   = linspace(0,1,7);   % 0, 1/6, ..., 1
symTicks = { ...
    '0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi'};

%% ---- Grid lines: constant A, B, C ----
for i = 2:6   % skip 0 and 1, avoid drawing on outer edges
    v = g(i);

    % constant A (alpha1)
    [x1,y1] = ternaryToCartesian(v,    0, 1-v);
    [x2,y2] = ternaryToCartesian(v, 1-v,   0);
    plot3([x1 x2],[y1 y2],[0 0],'Color',[0.9 0.9 0.9]);

    % constant B (alpha3)
    [x1,y1] = ternaryToCartesian(0,   v, 1-v);
    [x2,y2] = ternaryToCartesian(1-v, v,   0);
    plot3([x1 x2],[y1 y2],[0 0],'Color',[0.9 0.9 0.9]);

    % constant C (alpha2)
    [x1,y1] = ternaryToCartesian(0, 1-v,   v);
    [x2,y2] = ternaryToCartesian(1-v, 0,   v);
    plot3([x1 x2],[y1 y2],[0 0],'Color',[0.9 0.9 0.9]);
end

%% ---- Symbolic tick labels: 0, π/6, ..., π on z = 0 ----
for i = 1:numel(g)
    r     = g(i);          % = angle / π
    label = symTicks{i};

    % --- α3 bottom edge: left → right ---
    [xb,yb] = ternaryToCartesian(1-r, r, 0);
    text(xb, yb - 0.05, 0, label, 'FontSize', 14, ...
        'HorizontalAlignment','center','Interpreter','tex');

    % --- α1 left edge: top → bottom ---
    [x1,y1] = ternaryToCartesian(r, 0, 1-r);
    text(x1 - 0.05, y1, 0, label, 'FontSize', 14, ...
        'HorizontalAlignment','center','Rotation', 60,...
        'Interpreter','tex');

    % --- α2 right edge: bottom → top ---
    [x2,y2] = ternaryToCartesian(0, 1-r, r);
    text(x2 + 0.05, y2, 0, label, 'FontSize', 14, ...
        'HorizontalAlignment','center','Rotation', -60,...
        'Interpreter','tex');
end

end
