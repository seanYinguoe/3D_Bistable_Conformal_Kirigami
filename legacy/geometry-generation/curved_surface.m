function curved_surface()
close all; clc;

%% ===================== 1) Define the surface =====================
nu = 220; nv = 170;
uvec = linspace(-1.25, 1.25, nu);
vvec = linspace(-1.00, 1.00, nv);
[U,V] = meshgrid(uvec, vvec);

a = 1.6;
b = 1.6;

X = U;
Y = V;

% Hyperbolic paraboloid (analytical)
Z = (X.^2)/(a^2) - (Y.^2)/(b^2);

%% ===================== 2) Plot the surface =====================
fig = figure('Color','w','Position',[200 200 620 460]);
ax  = axes(fig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');

% main surface (lavender + translucent)
hSurf = surf(ax, X, Y, Z, ...
    'EdgeColor','none', ...
    'FaceColor',[0.82 0.83 0.92], ...
    'FaceAlpha',0.70);

%% ===================== 3) Grid lines on surface (EVEN spacing on surface) =====================
nU = 6; nV = 6;
gridRGBA = [0 0 0 0.85];
gridLW   = 1.2;

plot_surface_grid_even(ax, X, Y, Z, uvec, vvec, nU, nV, gridRGBA, gridLW, '-');

%% ===================== 4) Outline =====================
outlineLW = 1.2;
outlineColor = [0 0 0];

plot3(ax, X(1,:),   Y(1,:),   Z(1,:),   '-', 'Color',outlineColor,'LineWidth',outlineLW); % v=min
plot3(ax, X(end,:), Y(end,:), Z(end,:), '-', 'Color',outlineColor,'LineWidth',outlineLW); % v=max
plot3(ax, X(:,1),   Y(:,1),   Z(:,1),   '-', 'Color',outlineColor,'LineWidth',outlineLW); % u=min
plot3(ax, X(:,end), Y(:,end), Z(:,end), '-', 'Color',outlineColor,'LineWidth',outlineLW); % u=max

%% ===================== 5) Camera + lighting (paper style) =====================
view(ax, 315, 18);
camproj(ax,'orthographic');

material(ax,'dull');
lighting(ax,'gouraud');

camlight(ax,'headlight');
L = light(ax,'Style','infinite');
L.Position = [-1 1 1];

hSurf.AmbientStrength  = 0.55;
hSurf.DiffuseStrength  = 0.60;
hSurf.SpecularStrength = 0.06;
hSurf.SpecularExponent = 25;

xlim(ax, [min(X(:)) max(X(:))]);
ylim(ax, [min(Y(:)) max(Y(:))]);
zlim(ax, [min(Z(:)) max(Z(:))]);

%% ===================== 6) Export surface as OBJ =====================
% Build triangulation from grid
%% EXPORT OBJ
[nr,nc] = size(X);
Vtx = [X(:), Y(:), Z(:)];
F = [];

for i = 1:nr-1
    for j = 1:nc-1
        v1 = sub2ind([nr,nc], i,   j);
        v2 = sub2ind([nr,nc], i+1, j);
        v3 = sub2ind([nr,nc], i+1, j+1);
        v4 = sub2ind([nr,nc], i,   j+1);
        F = [F; v1 v2 v3; v1 v3 v4];
    end
end

write_obj('curved_surface.obj', Vtx, F);
disp('OBJ exported: curved_surface.obj');
end

%% ===================== helper: EVEN grid spacing ON THE SURFACE =====================
function plot_surface_grid_even(ax, X, Y, Z, uvec, vvec, nU, nV, col, lw, ls)
% Evenly divides the surface *by arc-length* along a representative mid-slice
% in each direction, then uses those u/v positions to draw constant-u / constant-v curves.
%
% This is much more "even" visually on a curved surface than uniform u,v spacing.

% ---- U-direction: choose u positions by equal arc-length along mid-v curve
midV = round(numel(vvec)/2);

Xu = X(midV,:);
Yu = Y(midV,:);
Zu = Z(midV,:);

dsU = sqrt(diff(Xu).^2 + diff(Yu).^2 + diff(Zu).^2);
sU  = [0, cumsum(dsU)];
if sU(end) < eps
    uEven = linspace(min(uvec), max(uvec), nU);
else
    sU = sU / sU(end); % normalize [0,1]
    uEven = interp1(sU, uvec, linspace(0,1,nU), 'linear', 'extrap');
end

for k = 1:numel(uEven)
    [~, iu] = min(abs(uvec - uEven(k)));
    plot3(ax, X(:,iu), Y(:,iu), Z(:,iu), ...
        'LineStyle', ls, 'Color', col, 'LineWidth', lw);
end

% ---- V-direction: choose v positions by equal arc-length along mid-u curve
midU = round(numel(uvec)/2);

Xv = X(:,midU);
Yv = Y(:,midU);
Zv = Z(:,midU);

dsV = sqrt(diff(Xv).^2 + diff(Yv).^2 + diff(Zv).^2);
sV  = [0; cumsum(dsV)];
if sV(end) < eps
    vEven = linspace(min(vvec), max(vvec), nV);
else
    sV = sV / sV(end);
    vEven = interp1(sV, vvec, linspace(0,1,nV), 'linear', 'extrap');
end

for k = 1:numel(vEven)
    [~, iv] = min(abs(vvec - vEven(k)));
    plot3(ax, X(iv,:), Y(iv,:), Z(iv,:), ...
        'LineStyle', ls, 'Color', col, 'LineWidth', lw);
end

end

function write_obj(filename, V, F)
fid = fopen(filename,'w');
if fid==-1
    error('Cannot open file for writing.');
end

for i = 1:size(V,1)
    fprintf(fid,'v %.8f %.8f %.8f\n',V(i,1),V(i,2),V(i,3));
end

for i = 1:size(F,1)
    fprintf(fid,'f %d %d %d\n',F(i,1),F(i,2),F(i,3));
end

fclose(fid);
end