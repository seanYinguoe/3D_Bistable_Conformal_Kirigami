function plot_ternary(alpha1, alpha2, alpha3, bist_strain, bistability)
% PLOT_TERNARY  Standalone ternary plot (no toolbox required)
%
% Inputs:
%   alpha1, alpha2, alpha3 : angles (must satisfy alpha1+alpha2+alpha3 = const)
%   bist_strain            : bistable strain (colour)
%   bistability            : bistability index (threshold < 0.1 → grey)
%
% This function draws:
%   • True ternary axes (three axes, ticks, gridlines)
%   • Grey monostable points (eta < 0.1)
%   • Colour bistable points using bistable strain

%% --- Normalize barycentric coordinates ---
S = alpha1 + alpha2 + alpha3;     % ~= pi
A = alpha1 ./ S;
B = alpha2 ./ S;
C = alpha3 ./ S;

% Convert barycentric → Cartesian for plotting
[X, Y] = ternaryToCartesian(A, B, C);

%% --- Threshold bistability ---
isBist = bistability >= 0.1 & isfinite(bist_strain);
isMono = ~isBist;

%% --- Start plotting ---
figure; hold on; axis equal; axis off;

% Draw ternary axes and grid
drawTernaryAxes();

% Plot monostable points (grey)
scatter(X(isMono), Y(isMono), 35, [0.8 0.8 0.8], 'filled');

% Plot bistable points (colored)
scatter(X(isBist), Y(isBist), 50, bist_strain(isBist), 'filled');

colormap(jet);
cb = colorbar;
ylabel(cb, '\epsilon_{bist}', 'Interpreter','tex');

% Labels
text(-0.05, -0.03, '\alpha_1', 'FontSize', 14, 'HorizontalAlignment','right');
text(1.05, -0.03, '\alpha_2', 'FontSize', 14, 'HorizontalAlignment','left');
text(0.5, sqrt(3)/2 + 0.05, '\alpha_3', 'FontSize', 14, 'HorizontalAlignment','center');

title('\epsilon_{bist} ternary map of (\alpha_1,\alpha_2,\alpha_3)', ...
      'Interpreter','tex','FontSize',14);

end


%% --- BARYCENTRIC → CARTESIAN --------------------------------------------
function [x, y] = ternaryToCartesian(A, B, C)
% Convert (A,B,C) barycentric coordinates into 2D triangle coordinates.

x = 0.5*(2*B + C);
y = (sqrt(3)/2)*C;

end


%% --- DRAW AXES ----------------------------------------------------------
function drawTernaryAxes()
% Draw a complete ternary axis system with ticks and gridlines.

% Outer triangle
plot([0 1 0.5 0], [0 0 sqrt(3)/2 0], 'k-', 'LineWidth', 2);
hold on;

% Tick settings
N = 5;   % number of tick/grid intervals
tickVals = linspace(0,1,N+1);

% Draw grid lines
for t = tickVals
    % Lines parallel to α1 axis
    [x1,y1] = ternaryToCartesian(1-t, t, 0);
    [x2,y2] = ternaryToCartesian(1-t, 0, t);
    plot([x1,x2], [y1,y2], 'Color', [0.85 0.85 0.85]);

    % Lines parallel to α2 axis
    [x1,y1] = ternaryToCartesian(t, 1-t, 0);
    [x2,y2] = ternaryToCartesian(0, 1-t, t);
    plot([x1,x2], [y1,y2], 'Color', [0.85 0.85 0.85]);

    % Lines parallel to α3 axis
    [x1,y1] = ternaryToCartesian(t, 0, 1-t);
    [x2,y2] = ternaryToCartesian(0, t, 1-t);
    plot([x1,x2], [y1,y2], 'Color', [0.85 0.85 0.85]);
end

end
