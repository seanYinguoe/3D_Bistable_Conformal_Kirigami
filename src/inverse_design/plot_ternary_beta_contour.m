function plot_ternary_beta_contour(a1, a2, a3, val, beta, eta, opts)
% PLOT_TERNARY_BETA_CONTOUR  3-D tube surface of the bistable region vs beta
%
% At every beta level the boundary of the bistable region is extracted as a
% closed ring in the ternary XY plane.  Rings at successive beta values are
% stitched together with surf() to form a smooth prism (柱体).  A Gaussian
% smooth is applied along the beta axis to remove the staircase between
% slices.
%
% Usage:
%   plot_ternary_beta_contour(a1, a2, a3, eps_bist, beta, eta_val)
%   plot_ternary_beta_contour(a1, a2, a3, eps_bist, beta, eta_val, opts)
%
% opts fields (all optional):
%   .threshold    - eps_bist cut-off for bistability   (default 0)
%   .etaThreshold - eta_val  cut-off                   (default 0.15)
%   .nGrid        - 2-D interpolation grid resolution  (default 150)
%   .nPts         - azimuthal points per ring          (default 200)
%   .faceAlpha    - surface transparency 0-1           (default 0.55)
%   .surfColor    - [r g b] surface colour             (default sky-blue)
%   .smoothSigma  - Gaussian smooth along beta axis    (default 1.5)
%   .symmetrize   - apply all 6 angle permutations     (default true)

%% ---- 0. Options -----------------------------------------------------------
if nargin < 7 || isempty(opts), opts = struct(); end
threshold    = getf(opts, 'threshold',    0);
etaThreshold = getf(opts, 'etaThreshold', 0.15);
nGrid        = getf(opts, 'nGrid',        150);
nPts         = getf(opts, 'nPts',         200);
faceAlpha    = getf(opts, 'faceAlpha',    0.55);
surfColor    = getf(opts, 'surfColor',    [0.35 0.72 0.42]);
smoothSigma  = getf(opts, 'smoothSigma',  1.5);
do_sym       = getf(opts, 'symmetrize',   true);
zScale       = getf(opts, 'zScale',       2);   % multiply Z data to stretch axis visually
alphaRange   = getf(opts, 'alphaRange',   [50, 70]);  % degrees [lo, hi]

% Rescaled ternary: the triangle vertex corresponds to alpha_vertex = 180-2*lo.
% Barycentric A' = (alpha - lo) / (alpha_vertex - lo) = (alpha_deg - lo) / (180-3*lo)
alpha_lo_deg  = alphaRange(1);
alpha_hi_deg  = alphaRange(2);
alpha_vtx_deg = 180 - 2*alpha_lo_deg;          % angle at each vertex (e.g. 80° for lo=50°)
alpha_denom   = alpha_vtx_deg - alpha_lo_deg;  % range in degrees (e.g. 30°)
alpha_lo_rad  = alpha_lo_deg  * pi/180;
alpha_denom_rad = alpha_denom * pi/180;

%% ---- 1. Eta filter --------------------------------------------------------
a1=a1(:); a2=a2(:); a3=a3(:); val=val(:); beta=beta(:); eta=eta(:);
keep = isfinite(eta) & (eta > etaThreshold);
a1=a1(keep); a2=a2(keep); a3=a3(keep); val=val(keep); beta=beta(keep);

%% ---- 2. Symmetrize --------------------------------------------------------
if do_sym
    [a1,a2,a3,val,beta] = symmetrize6(a1,a2,a3,val,beta);
end

%% ---- 3. Ternary -> Cartesian (rescaled to alphaRange) --------------------
% A' = (alpha1 - lo) / denom,  B' = (alpha3 - lo) / denom,  C' = (alpha2 - lo) / denom
% A'+B'+C' = (alpha1+alpha2+alpha3 - 3*lo) / denom = (pi - 3*lo) / denom = 1  (exact)
Ap = (a1 - alpha_lo_rad) ./ alpha_denom_rad;
Bp = (a3 - alpha_lo_rad) ./ alpha_denom_rad;
Cp = (a2 - alpha_lo_rad) ./ alpha_denom_rad;
[Xp,Yp] = ternary2cart(Ap, Bp, Cp);
ok = isfinite(val) & isfinite(beta) & isfinite(Xp) & isfinite(Yp);
Xp=Xp(ok); Yp=Yp(ok); beta=beta(ok); val=val(ok);

%% ---- 4. 2-D ternary grid --------------------------------------------------
xlin = linspace(0, 1,          nGrid);
ylin = linspace(0, sqrt(3)/2,  nGrid);
[Xg, Yg] = meshgrid(xlin, ylin);
in_tri   = inEquilateralTriangle(Xg, Yg);

%% ---- 5. Extract a closed boundary ring at every beta level ----------------
beta_unique = sort(unique(beta));
nB          = numel(beta_unique);

ring_x = NaN(nB, nPts);
ring_y = NaN(nB, nPts);
ring_z = NaN(nB, 1);

for ib = 1:nB
    b   = beta_unique(ib);
    sel = abs(beta - b) < 1e-12;
    xi  = Xp(sel); yi = Yp(sel); vi = val(sel);
    if numel(xi) < 3, continue; end

    Vg = griddata(xi, yi, vi, Xg, Yg, 'natural');
    Vg(~in_tri)  = threshold - 1;
    Vg(isnan(Vg)) = threshold - 1;

    C = contourc(xlin, ylin, Vg, [threshold threshold]);
    if isempty(C), continue; end

    [cx, cy] = getLargestContour(C);
    if numel(cx) < 4, continue; end

    [cx_r, cy_r]  = resampleContour(cx, cy, nPts);
    ring_x(ib,:)  = cx_r;
    ring_y(ib,:)  = cy_r;
    ring_z(ib)    = b;
end

%% ---- 6. Drop beta levels with no valid contour ----------------------------
valid  = ~isnan(ring_z);
ring_x = ring_x(valid, :);
ring_y = ring_y(valid, :);
ring_z = ring_z(valid);
nBv    = size(ring_x, 1);

if nBv < 2
    warning('plot_ternary_beta_contour: fewer than 2 valid beta levels — cannot build surface.');
    return;
end

%% ---- 7. Align rings to prevent surface twisting ---------------------------
for ib = 2:nBv
    [ring_x(ib,:), ring_y(ib,:)] = alignRing( ...
        ring_x(ib-1,:), ring_y(ib-1,:), ring_x(ib,:), ring_y(ib,:));
end

%% ---- 8. Smooth along beta axis (Gaussian) ---------------------------------
if smoothSigma > 0
    for j = 1:nPts
        ring_x(:,j) = gaussSmooth1D(ring_x(:,j), smoothSigma);
        ring_y(:,j) = gaussSmooth1D(ring_y(:,j), smoothSigma);
    end
end

%% ---- 9. Build surf matrices (close the ring by appending column 1) --------
CX = [ring_x, ring_x(:,1)];
CY = [ring_y, ring_y(:,1)];
CZ = repmat(ring_z, 1, nPts+1);

%% ---- 10. Draw -------------------------------------------------------------
z_max_data = max(ring_z) * 1.08;   % true beta ceiling (for labels)
z_max_plot = z_max_data * zScale;  % scaled Z used for all 3-D coordinates

% Scale Z coordinates — this physically stretches the axis
CZ_plot      = CZ       * zScale;
ring_z_plot  = ring_z   * zScale;

fig = figure;
ax  = axes(fig);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, 'off');
set(fig, 'Color', 'w');
set(ax, 'FontName', 'Arial');

% Triangular prism frame — pass both scaled height (for geometry) and
% original z_max (for tick labels so they still read in rad)
drawTriPrismFrame(z_max_plot, ax, alpha_lo_deg, alpha_hi_deg, alpha_vtx_deg, z_max_data);

% Tube walls
surf(CX, CY, CZ_plot, ...
    'FaceColor', surfColor, ...
    'FaceAlpha', faceAlpha, ...
    'EdgeColor', surfColor * 0.65, ...
    'EdgeAlpha', 0.20, ...
    'LineWidth', 0.3, ...
    'Parent', ax);

% Top cap
fill3(ring_x(end,:), ring_y(end,:), ring_z_plot(end)*ones(1,nPts), surfColor, ...
    'FaceAlpha', faceAlpha * 0.8, 'EdgeColor', 'none', 'Parent', ax);


%% ---- 11. View / limits ----------------------------------------------------
xlim(ax, [-0.15 1.15]);
ylim(ax, [-0.18 sqrt(3)/2 + 0.05]);
zlim(ax, [0 z_max_plot]);
view(ax, 40, 28);
end


%% =========================================================================
%  LOCAL HELPERS
%% =========================================================================

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

% ---- 6 angle permutations ------------------------------------------------
function [a1o,a2o,a3o,vo,bo] = symmetrize6(a1,a2,a3,v,b)
P = [1 2 3; 1 3 2; 2 1 3; 2 3 1; 3 1 2; 3 2 1];
A = [a1,a2,a3];
a1o=[]; a2o=[]; a3o=[]; vo=[]; bo=[];
for k = 1:6
    a1o=[a1o; A(:,P(k,1))]; %#ok<AGROW>
    a2o=[a2o; A(:,P(k,2))]; %#ok<AGROW>
    a3o=[a3o; A(:,P(k,3))]; %#ok<AGROW>
    vo=[vo; v];              %#ok<AGROW>
    bo=[bo; b];              %#ok<AGROW>
end
end

% ---- Ternary barycentric -> Cartesian ------------------------------------
function [x,y] = ternary2cart(A,B,C)
x = 0.5*(2*B+C);
y = (sqrt(3)/2)*C;
end

% ---- Point-in-equilateral-triangle (barycentric test) --------------------
function inside = inEquilateralTriangle(X,Y)
h = sqrt(3)/2;
C_b = Y./h;
B_b = X - Y/sqrt(3);
A_b = 1 - B_b - C_b;
inside = (A_b >= -1e-9) & (B_b >= -1e-9) & (C_b >= -1e-9);
end

% ---- Extract longest contour segment from contourc output ----------------
function [cx,cy] = getLargestContour(C)
cx=[]; cy=[]; max_n=0;
i=1;
while i < size(C,2)
    n = C(2,i);
    if n > max_n
        max_n = n;
        cx = C(1, i+1:i+n)';
        cy = C(2, i+1:i+n)';
    end
    i = i+n+1;
end
end

% ---- Resample contour to nPts equally spaced points (arc length) ---------
function [cx_r, cy_r] = resampleContour(cx, cy, nPts)
cx=cx(:); cy=cy(:);
% Remove duplicate endpoint if contourc closed it
if (cx(end)-cx(1))^2+(cy(end)-cy(1))^2 < 1e-14
    cx=cx(1:end-1); cy=cy(1:end-1);
end
% Arc length over closed loop
ds  = sqrt(diff([cx;cx(1)]).^2 + diff([cy;cy(1)]).^2);
arc = [0; cumsum(ds)];
t   = linspace(0, arc(end), nPts+1);
t   = t(1:end-1);   % nPts points, last = first omitted (periodic)
cx_r = interp1(arc, [cx;cx(1)], t, 'linear')';
cy_r = interp1(arc, [cy;cy(1)], t, 'linear')';
end

% ---- Align ring2 to ring1 by finding best cyclic shift + orientation -----
function [cx2a, cy2a] = alignRing(cx1,cy1,cx2,cy2)
% Try both orientations then pick the shift that minimises squared distance
[s1,d1] = bestShift(cx1,cy1,cx2,      cy2      );
[s2,d2] = bestShift(cx1,cy1,cx2(end:-1:1), cy2(end:-1:1));
if d1 <= d2
    cx2a = circshift(cx2,      s1);
    cy2a = circshift(cy2,      s1);
else
    cx2a = circshift(cx2(end:-1:1), s2);
    cy2a = circshift(cy2(end:-1:1), s2);
end
end

function [s, dist] = bestShift(cx1,cy1,cx2,cy2)
% FFT cross-correlation to find the cyclic shift maximising inner product
cc = real(ifft(fft(cx1(:).').*conj(fft(cx2(:).')))) + ...
     real(ifft(fft(cy1(:).').*conj(fft(cy2(:).'))));
[~,idx] = max(cc);
s  = idx-1;
cx2s = circshift(cx2, s);
cy2s = circshift(cy2, s);
dist = sum((cx1(:)-cx2s(:)).^2+(cy1(:)-cy2s(:)).^2);
end

% ---- 1-D Gaussian smooth with edge-padded convolution --------------------
function ys = gaussSmooth1D(y, sigma)
y=y(:);
if sigma <= 0 || numel(y) < 3, ys=y; return; end
r = ceil(3*sigma);
k = exp(-((-r:r).^2)/(2*sigma^2));
k = k/sum(k);
% Pad with edge values to reduce boundary error
yp = [repmat(y(1),r,1); y; repmat(y(end),r,1)];
ys = conv(yp', k, 'valid')';
end

% ---- Triangular prism frame: base triangle + vertical edges + beta axis --
function drawTriPrismFrame(z_max, ax, alpha_lo, alpha_hi, alpha_vtx, z_max_data)
% alpha_lo, alpha_hi, alpha_vtx : degrees defining the displayed range.
% The triangle vertex corresponds to alpha_vtx (= 180 - 2*alpha_lo).
% Ticks are shown from alpha_lo to alpha_hi in 5° steps.
h   = sqrt(3)/2;
grd = [0.80 0.80 0.80];
fsz = 13;
fsz_lbl = 16;
fnt = 'Arial';

% Triangle vertex positions (Cartesian, unchanged)
Vx = [0,   1,   0.5];
Vy = [0,   0,   h  ];

% Tick positions in rescaled barycentric [0,1]
% g=0 -> alpha_lo, g=1 -> alpha_vtx (= 180-2*alpha_lo)
tick_deg = alpha_lo : 5 : alpha_vtx;                       % e.g. 50,55,...,80
g_ticks  = (tick_deg - alpha_lo) / (alpha_vtx - alpha_lo); % 0 ... 1
tick_lbl = arrayfun(@(d) sprintf('%d\\circ', d), tick_deg, 'UniformOutput', false);

% Grid positions: same grid lines as ticks, skip 0 (drawn by triangle edge)
g_grid    = g_ticks(g_ticks > 0.01);
tick_len  = 0.030;   % outward tick length in Cartesian ternary units

%% -- Bottom triangle (thick) --
plot3(ax, [Vx, Vx(1)], [Vy, Vy(1)], [0,0,0,0], 'k-', 'LineWidth', 2);

%% -- Ternary grid at z = 0 --
for v = g_grid
    [x1,y1]=ternary2cart(v,  0,1-v); [x2,y2]=ternary2cart(v,1-v,0  ); plot3(ax,[x1 x2],[y1 y2],[0 0],'Color',grd);
    [x1,y1]=ternary2cart(0,  v,1-v); [x2,y2]=ternary2cart(1-v,v,0  ); plot3(ax,[x1 x2],[y1 y2],[0 0],'Color',grd);
    [x1,y1]=ternary2cart(0,1-v,v  ); [x2,y2]=ternary2cart(1-v,0,v  ); plot3(ax,[x1 x2],[y1 y2],[0 0],'Color',grd);
end

%% -- Outward tick marks + labels at z = 0 --
% Outward unit normals for each edge (away from triangle interior):
%   bottom edge  : (0, -1)
%   left edge    : (-√3/2,  1/2)
%   right edge   : ( √3/2,  1/2)
n_bot = [0,          -1      ];
n_lft = [-sqrt(3)/2,  0.5   ];
n_rgt = [ sqrt(3)/2,  0.5   ];

for i = 1:numel(g_ticks)
    r  = g_ticks(i);
    lb = tick_lbl{i};

    % --- α3 bottom edge ---
    [xb,yb] = ternary2cart(1-r, r, 0);
    plot3(ax, [xb, xb + n_bot(1)*tick_len], [yb, yb + n_bot(2)*tick_len], [0,0], 'k-', 'LineWidth', 1);
    text(ax, xb + n_bot(1)*0.055, yb + n_bot(2)*0.055, 0, lb, ...
        'FontSize',fsz,'FontName',fnt,'HorizontalAlignment','center','Interpreter','tex');

    % --- α1 left edge ---
    [x1,y1] = ternary2cart(r, 0, 1-r);
    plot3(ax, [x1, x1 + n_lft(1)*tick_len], [y1, y1 + n_lft(2)*tick_len], [0,0], 'k-', 'LineWidth', 1);
    text(ax, x1 + n_lft(1)*0.055, y1 + n_lft(2)*0.055, 0, lb, ...
        'FontSize',fsz,'FontName',fnt,'HorizontalAlignment','center','Rotation',60,'Interpreter','tex');

    % --- α2 right edge ---
    [x2,y2] = ternary2cart(0, 1-r, r);
    plot3(ax, [x2, x2 + n_rgt(1)*tick_len], [y2, y2 + n_rgt(2)*tick_len], [0,0], 'k-', 'LineWidth', 1);
    text(ax, x2 + n_rgt(1)*0.055, y2 + n_rgt(2)*0.055, 0, lb, ...
        'FontSize',fsz,'FontName',fnt,'HorizontalAlignment','center','Rotation',-60,'Interpreter','tex');
end

% Axis labels
text(ax,-0.07, h/2, 0, '\alpha_1','FontSize',fsz_lbl,'FontName',fnt,'Rotation', 60,'HorizontalAlignment','center','Interpreter','tex');
text(ax, 1.07, h/2, 0, '\alpha_2','FontSize',fsz_lbl,'FontName',fnt,'Rotation',-60,'HorizontalAlignment','center','Interpreter','tex');
text(ax, 0.50,-0.10, 0, '\alpha_3','FontSize',fsz_lbl,'FontName',fnt,'HorizontalAlignment','center','Interpreter','tex');

%% -- Three vertical edges --
for k = 1:3
    plot3(ax, [Vx(k) Vx(k)], [Vy(k) Vy(k)], [0 z_max], 'k-', 'LineWidth', 1.5);
end

%% -- Beta axis ticks + labels on vertex 1 edge (left-bottom) --
% Tick positions are in scaled Z; labels display original beta values.
n_btick          = 6;
beta_ticks_plot  = linspace(0, z_max,      n_btick);   % scaled positions
beta_ticks_label = linspace(0, z_max_data, n_btick);   % original values for labels
for i = 1:n_btick
    bt_p = beta_ticks_plot(i);
    bt_l = beta_ticks_label(i);
    plot3(ax, [Vx(1)-0.018, Vx(1)], [Vy(1), Vy(1)], [bt_p, bt_p], 'k-', 'LineWidth', 1);
    text(ax, Vx(1)-0.025, Vy(1), bt_p, sprintf('%.2g', bt_l), ...
        'FontSize', fsz, 'FontName', fnt, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end
text(ax, Vx(1)-0.10, Vy(1), z_max*0.5, '\beta  (rad)', ...
    'FontSize', fsz_lbl, 'FontName', fnt, 'Interpreter', 'tex', ...
    'HorizontalAlignment', 'center', 'Rotation', 90);

%% -- Top triangle (closes the prism) --
plot3(ax, [Vx, Vx(1)], [Vy, Vy(1)], z_max*[1,1,1,1], 'k-', 'LineWidth', 1.5);
end
