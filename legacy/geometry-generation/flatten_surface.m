function flatten_surface()
close all; clc;

%% --- Base "flat" parameter domain in z = x + i y
nu = 220; nv = 170;
xvec = linspace(-1.25, 1.25, nu);
yvec = linspace(-1.00, 1.00, nv);
[X,Y] = meshgrid(xvec, yvec);
Z = X + 1i*Y;   % complex grid

%% --- Conformal map w = f(z)
% Choose epsilon small so f'(z) != 0 everywhere in the domain.
eps = 0.08;

% Example 1: cubic perturbation (nice smooth curvy boundary)
W = -Z + 0.10*sin(Z);

% (Optional) add another harmonic for more organic but still conformal
% eps2 = 0.02;
% W = W + eps2*Z.^5;

Uf = real(W);
Vf = imag(W);

%% --- Pick 4x4 grid lines => 5 lines each direction (including boundary)
nDivX = 5; nDivY = 5;
colLines = unique(round(linspace(1, nu, nDivX+1)));
rowLines = unique(round(linspace(1, nv, nDivY+1)));

%% --- Plot style (lavender + black)
fig = figure('Color','w','Position',[220 220 560 450]);
ax  = axes(fig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');

% ---- boundary from mapped outer ring
bx = [Uf(1,:) Uf(2:end,end)' fliplr(Uf(end,1:end-1)) fliplr(Uf(2:end-1,1)')];
by = [Vf(1,:) Vf(2:end,end)' fliplr(Vf(end,1:end-1)) fliplr(Vf(2:end-1,1)')];

patch(ax, bx, by, [0.82 0.83 0.92], ...
    'FaceAlpha',0.70, 'EdgeColor','k', 'LineWidth',1.6);

gridCol = [0 0 0 0.85];
gridLW  = 1.2;

% ---- mapped grid lines (images of straight lines -> orthogonal curves)
for c = colLines
    plot(ax, Uf(:,c), Vf(:,c), 'Color',gridCol, 'LineWidth',gridLW);
end
for r = rowLines
    plot(ax, Uf(r,:), Vf(r,:), 'Color',gridCol, 'LineWidth',gridLW);
end

end