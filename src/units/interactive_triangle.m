% interactive_triangle
% Mode A (legacy):
%   interactive_triangle(edgeLen,l1,l4,beta,t)
% Mode B (anisotropic interactive deploy):
%   interactive_triangle(q1,q2,q3,edgeLen,l1,l4,beta,t,nD,Nseg,plot_flag)
function interactive_triangle(varargin)

if nargin == 5
    run_legacy_mode(varargin{:});
    return;
end

if nargin >= 10
    run_anisotropic_mode(varargin{:});
    return;
end

error('interactive_triangle: invalid input signature.');

end

function run_legacy_mode(edgeLen,l1,l4,beta,t)
% Original behavior preserved.
fig = figure('Position', [100, 100, 800, 800]);
set(fig, 'CloseRequestFcn', @safeCloseFigure);
xlim([-1.6*edgeLen, 0.4*edgeLen]);
ylim([-edgeLen, edgeLen]);

% Initial delta value
delta = 0;
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;

% Set the range of phi value
strain_min = 0;
strain_max = 0.8;

% Define colors
colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255};

% Create the triangle unit
[triangle,alpha_1_optimal,alpha_2_optimal] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);

% Add slider in the bottom half
uicontrol(fig, 'Style', 'slider', 'Min', strain_min*edgeLen, 'Max', strain_max*edgeLen, 'Value', delta, ...
    'Position', [150, 30, 500, 30], 'Callback', {@update_plot, fig});

uicontrol(fig, 'Style', 'text', 'String', ['Strain = ', num2str(delta/edgeLen)], ...
    'Position', [650, 25, 120, 40], 'FontSize', 12, 'Tag', 'strainText');

plot_triangle(triangle,colour);

    function update_plot(hObject, ~, fig_handle)
        delta_loc = get(hObject, 'Value');
        triangle_loc = triangle_unit(alpha_1_optimal, alpha_2_optimal, delta_loc, beta, edgeLen, l1,l4,t);

        clf(fig_handle);
        axis equal;
        hold on;
        xlim([-1.6*edgeLen, 0.4*edgeLen]);
        ylim([-edgeLen, edgeLen]);
        plot_triangle(triangle_loc,colour);

        uicontrol(fig_handle, 'Style', 'slider', 'Min', strain_min*edgeLen, 'Max', strain_max*edgeLen, 'Value', delta_loc, ...
            'Position', [150, 30, 500, 30], 'Callback', {@update_plot, fig_handle});

        uicontrol(fig_handle, 'Style', 'text', 'String', ['Strain = ', num2str(delta_loc/edgeLen)], ...
            'Position', [650, 25, 120, 40], 'FontSize', 12, 'Tag', 'strainText');
    end
end

function run_anisotropic_mode(q1,q2,q3,edgeLen,l1,l4,beta,t,nD,Nseg,plot_flag)
if nargin < 11 || isempty(plot_flag)
    plot_flag = true;
end

% Ensure row vectors
q1 = q1(:).';
q2 = q2(:).';
q3 = q3(:).';

% Helpers
wrap = @(th) atan2(sin(th), cos(th));
pack_x = @(phiB,eB,phiA,eA,phiC,eC,theta,xG,yG) [phiB; eB; phiA; eA; phiC; eC; theta; xG; yG];

% Geometry constants (same as deform_triangle_anisotropic)
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(pi/3 - beta);
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) - cos(pi/3 - beta) .* l4;
l5 = ((sqrt(3)/2) .* l4 + sin(beta) .* l6) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 - l2 .* ((sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)));
r_vertex = (sqrt(3)/3) * l3;

Emod = 4.3e11;
b = 1.0;

% Original/deformed boundary in same convention
p1_orig = [0, 0];
p2_orig = [0, -edgeLen];
p3_orig = [-sqrt(3)/2*edgeLen, -1/2*edgeLen];

edge1 = norm(q2 - q3);
edge2 = norm(q1 - q3);
edge3 = norm(q1 - q2);

p1_def = [0, 0];
p2_def = [0, -edge3];
p3_y = (edge1^2 - edge2^2 - edge3^2) / (2 * edge3);
p3_x = -sqrt(max(edge2^2 - p3_y^2, 0));
p3_def = [p3_x, p3_y];

% Reference unit for flank transforms
[tri0, ~] = deform_triangle_isotropic(0, edgeLen, l1, l4, beta, t, Nseg);

% Optimization warm-start state (same as deform_triangle_anisotropic)
theta_prev = -pi + beta;
K = Nseg - 1;
DeltaTot = wrap(theta_prev + pi - beta);
phi0 = (DeltaTot / K) * ones(K,1);
e0 = zeros(Nseg,1);
G0 = [-sqrt(3)/6*edgeLen, -1/2*edgeLen];
xG0 = G0(1);
yG0 = G0(2);
x_prev = pack_x(phi0,e0,phi0,e0,phi0,e0,theta_prev,xG0,yG0);

% Histories
E_curve = nan(nD,1);
alpha_curve = nan(nD,1);
pack_hist = cell(nD,1);
flankB_hist = cell(nD,1);
flankA_hist = cell(nD,1);
flankC_hist = cell(nD,1);

for k = 1:nD
    alpha_curve(k) = (k - 1) / (nD - 1);

    % Interpolate boundary
    p1 = (1 - alpha_curve(k)) * p1_orig + alpha_curve(k) * p1_def;
    p2 = (1 - alpha_curve(k)) * p2_orig + alpha_curve(k) * p2_def;
    p3 = (1 - alpha_curve(k)) * p3_orig + alpha_curve(k) * p3_def;

    [B_flank, A_flank, C_flank, alphaL_B, alphaL_A, alphaL_C, flank_B, flank_A, flank_C] = get_flank(p1, p2, p3);

    params = struct('t',t,'E',Emod,'b',b,'beta',beta,'B_flank',B_flank, ...
        'A_flank',A_flank,'C_flank',C_flank,'r_vertex',r_vertex,'l2',l2, ...
        'l3',l3,'alphaL_B',alphaL_B,'alphaL_A',alphaL_A,'alphaL_C',alphaL_C, ...
        'xG0',xG0,'yG0',yG0);

    [energy, pack] = energy_lig_anisotropic(x_prev, Nseg, params);

    E_curve(k) = energy;
    pack_hist{k} = pack;
    flankB_hist{k} = flank_B;
    flankA_hist{k} = flank_A;
    flankC_hist{k} = flank_C;

    theta_prev = pack.theta;
    xG0 = pack.xG;
    yG0 = pack.yG;
    x_prev = pack_x(pack.phiB,pack.eB,pack.phiA,pack.eA,pack.phiC,pack.eC,theta_prev,xG0,yG0);
end

% Figure
fig = figure('Color','w', 'Position',[120 120 1200 560]);
set(fig, 'CloseRequestFcn', @safeCloseFigure);
setappdata(fig, 'isClosing', false);
ax1 = subplot(1,2,1, 'Parent', fig);
ax2 = subplot(1,2,2, 'Parent', fig);

k0 = 1;
plot_unit_step(ax1, k0);
plot_energy(ax2, alpha_curve, E_curve, k0);

uicontrol(fig, 'Style','slider', 'Min',1, 'Max',nD, 'Value',k0, ...
    'SliderStep',[1/max(1,nD-1), 10/max(1,nD-1)], ...
    'Position',[120 20 560 28], 'Callback', @onSlide, ...
    'Interruptible','off', 'BusyAction','cancel');

txt = uicontrol(fig, 'Style','text', 'String', sprintf('k=%d, alpha=%.4f', k0, alpha_curve(k0)), ...
    'Position',[700 16 220 30], 'FontSize',11, 'HorizontalAlignment','left');

if plot_flag
    uicontrol(fig, 'Style','pushbutton', 'String','Open detailed state figure', ...
        'Position',[920 16 220 30], 'Callback', @showDetailed);
end

    function onSlide(hObj, ~)
        if ~isvalid(fig) || getappdata(fig, 'isClosing')
            return;
        end
        k = max(1, min(nD, round(get(hObj, 'Value'))));
        plot_unit_step(ax1, k);
        plot_energy(ax2, alpha_curve, E_curve, k);
        set(txt, 'String', sprintf('k=%d, alpha=%.4f', k, alpha_curve(k)));
        drawnow limitrate;
    end

    function showDetailed(~, ~)
        if ~isvalid(fig) || getappdata(fig, 'isClosing')
            return;
        end
        k = get(ax2, 'UserData');
        if isempty(k) || ~isscalar(k)
            k = 1;
        end
        a = alpha_curve(k);
        p1 = (1 - a) * p1_orig + a * p1_def;
        p2 = (1 - a) * p2_orig + a * p2_def;
        p3 = (1 - a) * p3_orig + a * p3_def;
        q1k = [p1, 0]; q2k = [p2, 0]; q3k = [p3, 0];
        deform_triangle_anisotropic(q1k, q2k, q3k, edgeLen, l1, l4, beta, t, nD, Nseg, true);
    end

    function plot_unit_step(ax, k)
        cla(ax);
        axes(ax); %#ok<LAXES>
        hold(ax, 'on');
        axis(ax, 'equal');

        colour = {'white', [0.9216 0.8863 0.4235], [0.7059 0.9608 0.4118], [0.9216 0.8863 0.4235]};
        pack = pack_hist{k};
        flank_B = flankB_hist{k};
        flank_A = flankA_hist{k};
        flank_C = flankC_hist{k};

        update_triangle_local(pack.B, pack.A, pack.C, flank_B, flank_A, flank_C, ...
            pack.XYB, pack.XYA, pack.XYC, colour, ax, t, edgeLen);

        title(ax, sprintf('Unit configuration at alpha = %.4f', alpha_curve(k)));
        axis(ax, 'off');
    end

    function [B_flank_def, A_flank_def, C_flank_def, alphaL_B_def, alphaL_A_def, alphaL_C_def, flank1, flank2, flank3] = get_flank(p1, p2, p3)
        f_flank = [19 20 21 22;
                   27 28 29 30;
                   23 24 25 26];

        p = [p1; p2; p3];
        base_idx = [22 30 26];
        edge_target = [p3 - p1; p1 - p2; p2 - p3];
        flanks = cell(3,1);

        for ii = 1:3
            flk = tri0(f_flank(ii,:),:) + (p(ii,:) - tri0(base_idx(ii),:));
            v1 = flk(3,:) - flk(4,:);
            v2 = edge_target(ii,:);
            rot_ang = atan2(v1(1)*v2(2) - v1(2)*v2(1), dot(v1, v2));
            flanks{ii} = real((flk - p(ii,:)) * rotation(-rot_ang) + p(ii,:));
        end

        flank1 = flanks{1};
        flank2 = flanks{2};
        flank3 = flanks{3};

        B_flank_def = flank1(2,:);
        A_flank_def = flank2(2,:);
        C_flank_def = flank3(2,:);

        vB = flank1(2,:) - flank1(3,:);
        alphaL_B_def = atan2(vB(2), vB(1));
        vA = flank2(2,:) - flank2(3,:);
        alphaL_A_def = atan2(vA(2), vA(1));
        vC = flank3(2,:) - flank3(3,:);
        alphaL_C_def = atan2(vC(2), vC(1));
    end
end

function safeCloseFigure(src, ~)
% Robust close handler to avoid UI freeze when callbacks are active.
try
    if isvalid(src)
        setappdata(src, 'isClosing', true);
    end
catch
end
delete(src);
end

function triangle = update_triangle_local(B,A,C,flankB,flankA,flankC,XYB,XYA,XYC,colour,ax,t,edgeLen)
triangle = zeros(45,2);
triangle(43,:) = A;
triangle(44,:) = B;
triangle(45,:) = C;
triangle([19,20,21,22],:) = flankB;
triangle([27,28,29,30],:) = flankA;
triangle([23,24,25,26],:) = flankC;

axes(ax); %#ok<LAXES>
plot_triangle(triangle, colour);
build_ligament_local(XYB, flankB, A, B, t, colour{3}, ax);
build_ligament_local(XYA, flankA, C, A, t, colour{3}, ax);
build_ligament_local(XYC, flankC, B, C, t, colour{3}, ax);

xlim(ax, [-1.4*edgeLen, 0.5*edgeLen]);
ylim(ax, [-1.1*edgeLen, 0.25*edgeLen]);
end

function [XY_inner, XY_outer] = build_ligament_local(XY_outer, flank, A, B, t, colour, ax)
if nargin < 6
    colour = [0.5 0.5 0.5];
end

n = size(XY_outer,1);

dir_start = flank(1,:) - flank(2,:);
dir_start = dir_start / norm(dir_start);

dir_end = A - B;
dir_end = dir_end / norm(dir_end);

dir_interp = zeros(n,2);
for i = 1:n
    s = (i - 1) / (n - 1);
    dir = (1 - s) * dir_start + s * dir_end;
    dir_interp(i,:) = dir / norm(dir);
end

XY_inner_1 = XY_outer + t * dir_interp;
XY_inner_2 = XY_outer - t * dir_interp;

target = 0.5 * (A + B);
d1 = norm(mean(XY_inner_1,1) - target);
d2 = norm(mean(XY_inner_2,1) - target);

if d1 < d2
    XY_inner = XY_inner_1;
else
    XY_inner = XY_inner_2;
end

verts = [XY_outer; flipud(XY_inner)];
faces = 1:(2*n);

patch('Parent',ax, 'Vertices',verts, 'Faces',faces, ...
    'FaceColor',colour, 'FaceAlpha',0.5, 'EdgeColor','k', 'LineWidth',0.5);
end

function plot_energy(ax, alpha, E, k)
cla(ax);
plot(ax, alpha, E, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
hold(ax, 'on');
plot(ax, alpha(k), E(k), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
set(ax, 'UserData', k);
xlabel(ax, 'Deployment alpha');
ylabel(ax, 'Strain Energy');
grid(ax, 'on');
axis(ax, 'tight');
title(ax, 'Energy curve (marker = current step)');
end
