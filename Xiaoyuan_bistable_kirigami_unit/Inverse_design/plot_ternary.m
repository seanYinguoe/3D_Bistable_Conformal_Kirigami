function plot_ternary(alpha1, alpha2, alpha3, bist_strain, bistability)
% PLOT_TERNARY  Standalone symmetric ternary plot (no toolbox required)
% Now includes symbolic radian ticks: 0, π/6, π/3, π/2, 2π/3, 5π/6, π

%% ---- 1) Symmetrize points across all 6 permutations (full triangle) ----
A  = [alpha1(:), alpha2(:), alpha3(:)];
E1 = bist_strain(:);
E2 = bistability(:);

P = [
    1 2 3;
    1 3 2;
    2 1 3;
    2 3 1;
    3 1 2;
    3 2 1
];

Afull  = [];
E1full = [];
E2full = [];

for i = 1:6
    Afull  = [Afull;  A(:,P(i,:))];
    E1full = [E1full; E1];
    E2full = [E2full; E2];
end

alpha1 = Afull(:,1);
alpha2 = Afull(:,2);
alpha3 = Afull(:,3);
bist_strain = E1full;
bistability = E2full;

%% ---- 2) Convert to barycentric coordinates ----
S = alpha1 + alpha2 + alpha3;
A = alpha1 ./ S;
B = alpha2 ./ S;
C = alpha3 ./ S;

[X, Y] = ternaryToCartesian(A,B,C);

%% ---- 3) Threshold bistability ----
isBist = bistability >= 0.1 & isfinite(bist_strain);
isMono = ~isBist;

%% ---- 4) Plot frame, grid, axis ticks ----
figure; hold on; axis equal; axis off;

drawTernaryAxesWithSymbolicTicks();

%% ---- 5) Draw points ----
scatter(X(isMono), Y(isMono), 35, [0.8 0.8 0.8], 'filled');
scatter(X(isBist), Y(isBist), 55, bist_strain(isBist), 'filled');

colormap(jet);
cb = colorbar;
ylabel(cb, '\epsilon_{bist}', 'Interpreter','tex');

%% ---- 6) Axis labels ----
text(0.05, sqrt(3)/4, '\alpha_1', 'FontSize', 18, ...
    'Rotation', 60, 'HorizontalAlignment','center');
text(0.95, sqrt(3)/4, '\alpha_2', 'FontSize', 18, ...
    'Rotation', -60, 'HorizontalAlignment','center');
text(0.50, -0.12, '\alpha_3', 'FontSize', 18, ...
    'HorizontalAlignment','center');

title('\epsilon_{bist} ternary map of (\alpha_1,\alpha_2,\alpha_3)',...
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
% DRAW TERNARY AXES + GRID + SYMBOLIC RADIAN TICKS
%% ------------------------------------------------------------------------
function drawTernaryAxesWithSymbolicTicks()

% Outer triangle
plot([0 1 0.5 0], [0 0 sqrt(3)/2 0], 'k-', 'LineWidth', 2);
hold on;

%% ---- Grid lines ----
g = linspace(0,1,7);   % 6 intervals (π/6)
for t = g
    [x1,y1] = ternaryToCartesian(1-t,t,0);
    [x2,y2] = ternaryToCartesian(1-t,0,t);
    plot([x1 x2],[y1 y2], 'Color',[0.9 0.9 0.9]);

    [x1,y1] = ternaryToCartesian(t,1-t,0);
    [x2,y2] = ternaryToCartesian(0,1-t,t);
    plot([x1 x2],[y1 y2], 'Color',[0.9 0.9 0.9]);

    [x1,y1] = ternaryToCartesian(t,0,1-t);
    [x2,y2] = ternaryToCartesian(0,t,1-t);
    plot([x1 x2],[y1 y2], 'Color',[0.9 0.9 0.9]);
end

%% ---- Symbolic tick labels: 0, π/6, ..., π ----
symTicks = { ...
    '0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi'};

rt3 = sqrt(3);

for i = 1:numel(g)
    r = g(i);
    label = symTicks{i};

    %% --- α3 bottom edge: left → right ---
    xb = r;
    yb = 0;
    text(xb, yb - 0.05, label, 'FontSize', 18, ...
        'HorizontalAlignment','center','Interpreter','tex');

    %% --- α1 left edge: bottom-left → top ---
    x1 = 0.5*r;
    y1 = (rt3/2)*r;
    text(x1 - 0.05, y1, symTicks{end+1-i}, 'FontSize', 18, ...
        'HorizontalAlignment','center','Rotation', 60,...
        'Interpreter','tex');

    %% --- α2 right edge: bottom-right → top ---
    x2 = 1 - 0.5*r;
    y2 = (rt3/2)*r;
    text(x2 + 0.05, y2, label, 'FontSize', 18, ...
        'HorizontalAlignment','center','Rotation', -60,...
        'Interpreter','tex');
end

end
