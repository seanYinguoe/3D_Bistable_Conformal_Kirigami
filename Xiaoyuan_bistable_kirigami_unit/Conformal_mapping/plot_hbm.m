function plot_hbm(sol, varargin)
% Visualize Hencky bar-chain:
%   - undeformed & deformed centerlines
%   - internal rotational springs at joints (proportional to |kappa|)
%   - optional end springs (angle mismatch)
%
% Inputs
%   sol : struct returned by hbm_energy (fields: XY0,XY1,phi0,phi1,kappa0,kappa1,a,EI, etc.)
%
% Name-Value options
%   'ShowEndSprings'   : true/false (default: true)
%   'Turns'            : number of turns for each spring coil (default: 4)
%   'SpringScale'      : scalar for spring size (default: 0.6)
%   'ColorMap'         : 'signed' (blue/red by sign of kappa) or 'mono' (single color), default 'signed'

    p = inputParser;
    addParameter(p, 'ShowEndSprings', true, @(x)islogical(x));
    addParameter(p, 'Turns', 4, @(x)isnumeric(x) && x>=1);
    addParameter(p, 'SpringScale', 0.6, @(x)isnumeric(x) && x>0);
    addParameter(p, 'ColorMap', 'signed', @(s)ischar(s) || isstring(s));
    parse(p, varargin{:});
    opt = p.Results;

    % figure
    figure('Color','w'); tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % ===== Left panel: shapes =====
    ax1 = nexttile; hold(ax1,'on'); axis(ax1,'equal'); grid(ax1,'on');
    title(ax1, sprintf('Bar-chain shape  |  Eb(deformed)=%.3e', sol.Ebend1));
    xlabel(ax1,'x'); ylabel(ax1,'y');

    % undeformed centerline
    plot(ax1, sol.XY0(:,1), sol.XY0(:,2), '-', 'LineWidth',2, 'Color',[0.25 0.60 1.00], ...
         'DisplayName','undeformed centerline');
    draw_links(ax1, sol.XY0, sol.phi0, [0.25 0.60 1.00]*0.85);

    % deformed centerline
    plot(ax1, sol.XY1(:,1), sol.XY1(:,2), '-', 'LineWidth',2, 'Color',[1.00 0.40 0.25], ...
         'DisplayName','deformed centerline');
    draw_links(ax1, sol.XY1, sol.phi1, [1.00 0.40 0.25]*0.85);

    % end markers
    plot(ax1, sol.A0(1),sol.A0(2),'ko','MarkerFaceColor','k','DisplayName','A0');
    plot(ax1, sol.B0(1),sol.B0(2),'ks','MarkerFaceColor','k','DisplayName','B0');
    plot(ax1, sol.A1(1),sol.A1(2),'mo','MarkerFaceColor','m','DisplayName','A1');
    plot(ax1, sol.B1(1),sol.B1(2),'ms','MarkerFaceColor','m','DisplayName','B1');

    legend(ax1,'Location','bestoutside');

    % internal springs on deformed shape
    draw_internal_springs(ax1, sol.XY1, sol.phi1, sol.kappa1, sol.a, ...
                          'Turns', opt.Turns, 'SpringScale', opt.SpringScale, 'ColorMap', opt.ColorMap);

    % optional end springs
    if opt.ShowEndSprings
        draw_end_springs(ax1, sol, 'SpringScale', opt.SpringScale);
    end

    % ===== Right panel: angles & curvature =====
    ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
    title(ax2, '\phi (segment angle) and \kappa (discrete curvature)');
    xlabel(ax2,'normalized arc-length s/L'); ylabel(ax2,'rad');

    N  = numel(sol.phi1);
    a  = sol.a;
    s  = (0:N)'*a; sN = s/s(end);
    sc = ((0.5:1:N-0.5)')*a; scN = sc/s(end);     % centers for phi
    sk = (1:N-1)'*a;           skN = sk/s(end);   % joints for kappa

    plot(ax2, scN, sol.phi0(:), '-o', 'LineWidth',1.6, 'MarkerSize',4, 'Color',[0.25 0.60 1.00], ...
         'DisplayName','\phi (undeformed)');
    plot(ax2, scN, sol.phi1(:), '-o', 'LineWidth',1.6, 'MarkerSize',4, 'Color',[1.00 0.40 0.25], ...
         'DisplayName','\phi (deformed)');

    plot(ax2, skN, sol.kappa1(:), '-s', 'LineWidth',1.6, 'MarkerSize',4, 'Color',[0.20 0.20 0.20], ...
         'DisplayName','\kappa (deformed)');

    legend(ax2,'Location','best');
end

% ---------- helpers ----------

function draw_links(ax, XY, phi, color)
    N = numel(phi);
    for i=1:N
        plot(ax, XY(i:i+1,1), XY(i:i+1,2), '-', 'LineWidth',3, 'Color',color, 'HandleVisibility','off');
        % tiny orientation tick
        mid = 0.5*(XY(i,:)+XY(i+1,:));
        q = 0.25*norm(XY(end,:)-XY(1,:))/N;
        quiver(ax, mid(1),mid(2), q*cos(phi(i)), q*sin(phi(i)), 0, ...
               'MaxHeadSize',0.8, 'LineWidth',0.8, 'Color',color, 'HandleVisibility','off');
    end
end

function draw_internal_springs(ax, XY, phi, kappa, a, varargin)
% Draw a coil at each interior joint (i = 1..N-1) centered at joint node,
% oriented along the local normal to indicate a torsional spring at the hinge.
    p = inputParser;
    addParameter(p,'Turns',4,@(x)isnumeric(x) && x>=1);
    addParameter(p,'SpringScale',0.6,@(x)isnumeric(x) && x>0);
    addParameter(p,'ColorMap','signed');
    parse(p,varargin{:});
    Turns = p.Results.Turns;
    sscale = p.Results.SpringScale;
    cmap   = p.Results.ColorMap;

    N  = numel(phi);
    if N <= 1, return; end

    % color function by sign/magnitude of kappa
    function col = kappa_color(val)
        switch lower(cmap)
            case 'signed'
                % blue for negative, red for positive, intensity by magnitude
                m = min(1, 0.5 + 0.5*min(1, abs(val)/max(1e-12, max(abs(kappa)))));
                if val >= 0, col = [1.0 m*0.6 m*0.6];
                else,        col = [m*0.6 m*0.6 1.0];
                end
            otherwise
                col = [0.1 0.1 0.1];
        end
    end

    Ltot = norm(XY(end,:)-XY(1,:));
    base_len = sscale * (Ltot / N);    % coil length scale

    for i = 1:N-1
        % joint at node i+1 (between segment i and i+1)
        P  = XY(i+1,:);

        % average tangent at the joint
        t1 = [cos(phi(i)),     sin(phi(i))];
        t2 = [cos(phi(i+1)),   sin(phi(i+1))];
        t  = (t1 + t2);  if norm(t)<1e-12, t = t1; end
        t  = t / norm(t);

        % local normal to place a "torsional coil" sticking out
        n  = [-t(2), t(1)];

        % coil geometry (centered at P, along n)
        k  = kappa(i);
        len = base_len * (0.6 + 0.4*min(1, abs(k)/(max(abs(kappa))+1e-12))); % scale with |kappa|
        amp = 0.35*len;  % amplitude
        coil = coil_poly(P, n, len, amp, Turns);

        plot(ax, coil(:,1), coil(:,2), '-', 'LineWidth',1.8, 'Color',kappa_color(k), ...
             'DisplayName', ternary(i==1, 'internal springs', ''));
    end
end

function XY = coil_poly(P, dir, L, A, turns)
% Create a simple 2D coil polyline starting at P, extending along "dir".
% L: coil length; A: amplitude; turns: number of oscillations.
    dir = dir / norm(dir);
    % choose a perpendicular vector for the coil oscillation plane
    per = [-dir(2), dir(1)];
    npts = max(20, 8*turns);
    s = linspace(-L/2, L/2, npts);
    w = 2*pi*turns / L;
    XY = P + s.' .* dir + (A*sin(w*s)).' .* per;
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function draw_end_springs(ax, sol, varargin)
% Visualize end torsional springs as small coils at A1 and B1,
% sized by the mismatch between end link angle and prescribed clamp angle.
    p = inputParser;
    addParameter(p,'SpringScale',0.6,@(x)isnumeric(x) && x>0);
    parse(p,varargin{:});
    sscale = p.Results.SpringScale;

    % left end
    phi1 = sol.phi1(1);
    dthL = wrapToPi(phi1 - sol.thetaA1);
    tL = [cos(phi1), sin(phi1)]; nL = [-tL(2), tL(1)];
    len = sscale * (sol.a * sol.N)/sol.N; amp = 0.35*len;
    cL = coil_poly(sol.XY1(1,:), nL, len*(0.6+0.4*min(1,abs(dthL)/0.5)), amp, 3);
    plot(ax, cL(:,1), cL(:,2), '-', 'LineWidth',1.8, 'Color',[0 0 0], 'DisplayName','end springs');

    % right end
    phiN = sol.phi1(end);
    dthR = wrapToPi(phiN - sol.thetaB1);
    tR = [cos(phiN), sin(phiN)]; nR = [-tR(2), tR(1)];
    cR = coil_poly(sol.XY1(end,:), nR, len*(0.6+0.4*min(1,abs(dthR)/0.5)), amp, 3);
    plot(ax, cR(:,1), cR(:,2), '-', 'LineWidth',1.8, 'Color',[0 0 0], 'HandleVisibility','off');
end
