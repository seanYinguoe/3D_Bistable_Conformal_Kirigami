function [E_total,alpha] = deform_triangle_anisotropic(q1,q2,q3,edgeLen,l1,l4,beta,t,nD,N,plot_flag)
% Inputs:
%   q1, q2, q3 - Unit node coordinates (1x3 vectors)
%   edgeLen    - Original triangle edge length
%   l1, l4     - Geometric parameters of flanks
%   beta       - Tilting angle
%   t          - Thickness of filaments
%   nD         - Number of deployment steps
%   N          - Number of nodes per ligament
%   plot_flag  - Optional, true to plot energy curve and deployed unit
%
% Output:
%   E_total    - Deployed energy
%   alpha      - Deployment fraction (0 to 1)

if nargin < 11 || isempty(plot_flag)
    plot_flag = false;
end

%% Helper
wrap  = @(th) atan2(sin(th), cos(th));
pack_x = @(phiB,eB,phiA,eA,phiC,eC,theta,xG,yG) [phiB; eB; phiA; eA; phiC; eC; theta; xG; yG];

%% Geometry constants (independent of delta)
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(pi/3 - beta);
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) - cos(pi/3 - beta) .* l4;
l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 ...
    - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) );
r_vertex  = (sqrt(3)/3) * l3;     % centroid–vertex radius (r = l3/sqrt(3))

Emod = 4.3e11;  b = 1.0;          % material constants

%% Outer-edge lengths
% Original unit
p1_orig = [0,0];
p2_orig = [0,-edgeLen];
p3_orig = [-sqrt(3)/2*edgeLen, -1/2*edgeLen];
% current unit
edge1 = norm(q2 - q3);  % p2p3
edge2 = norm(q1 - q3);  % p1p3
edge3 = norm(q1 - q2);  % p1p2
p1_def = [0,0];
p2_def = [0,-edge3];
p3_y = (edge1^2 - edge2^2 - edge3^2) / (2 * edge3);
p3_x = -sqrt(max(edge2^2 - p3_y^2,0));
p3_def = [p3_x,p3_y];

%% create reference unit(delta = 0)
[tri0, ~] = deform_triangle_isotropic(0, edgeLen, l1, l4, beta, t, N);
E_total  = nan(1, nD);
theta_all = nan(1, nD);
theta_prev = -pi + beta;             % start with undefomred one

%% input initial guess for optimisation
K = N-1;  % number of torsional spring
DeltaTot = wrap(theta_prev + pi - beta); % total angle difference
phi0 = (DeltaTot / K) * ones(K,1); % initial rotational angle of torsional springs
e0   = zeros(N,1);                 % initial length change in ligaments
G0 = [-sqrt(3)/6*edgeLen, -1/2*edgeLen];
xG0 = G0(1);
yG0 = G0(2);
x_prev = pack_x(phi0,e0,phi0,e0,phi0,e0,theta_prev,xG0,yG0); % initial guess
alpha = zeros(nD,1);
%% Loop over deltas
for k = 1:nD
    % deploy unit using interpolation between original and target one
    alpha(k) = (k - 1) / (nD - 1); % interpolation fraction 0 → 1

    % Interpolate node positions
    p1 = (1 - alpha(k)) * p1_orig + alpha(k) * p1_def;
    p2 = (1 - alpha(k)) * p2_orig + alpha(k) * p2_def;
    p3 = (1 - alpha(k)) * p3_orig + alpha(k) * p3_def;
    % Get
    [B_flank, A_flank, C_flank,...
        alphaL_B, alphaL_A, alphaL_C,...
        flank_B, flank_A, flank_C] = get_flank(p1, p2, p3);

    % update params for optimisation
    params = struct('t',t,'E',Emod,'b',b,'beta',beta,'B_flank',B_flank,...
        'A_flank',A_flank,'C_flank',C_flank,'r_vertex',r_vertex,'l2',l2,...
        'l3',l3,'alphaL_B',alphaL_B,'alphaL_A',alphaL_A,'alphaL_C',alphaL_C,...
        'xG0',xG0,'yG0',yG0);

    [energy, pack] = energy_lig_anisotropic(x_prev, N, params);
    assert(pack.diagnostics.exitflag>0&&pack.diagnostics.max_constraint<1e-6, ...
        'kirigami:EnergySolveFailed','Unit energy step %d failed its solver checks.',k); % N=10 segments
    % Save the reuslt
    theta_all(k) = pack.theta;
    E_total(k)   = energy;               % 3 ligaments total
    theta_prev = pack.theta;
    xG0 = pack.xG;
    yG0 = pack.yG;
    x_prev = pack_x(pack.phiB,pack.eB,pack.phiA,pack.eA,pack.phiC,pack.eC,theta_prev,xG0,yG0);
    last_pack = pack;
end

%% Plot energy curve
if plot_flag
    figure('Color','w');
    hold on; box on;
    plot(alpha, E_total, '-', 'Color',[0.85 0.33 0.10],'LineWidth', 1, ...
        'DisplayName', 'Unit energy');
    xlabel('Deployment', 'Interpreter','tex', 'FontSize',20);
    ylabel('Strain Energy(N/mm^2)', 'Interpreter','tex', 'FontSize',20);
    set(gca, 'FontName','Times New Roman','FontSize',20);
    legend('Location','northwest','Box','off', 'Fontsize',18);
    grid off;
    axis square;
end

%% Plot configuration
if plot_flag
    figure()
    colour = {'white', [0.9216 0.8863 0.4235], [0.7059 0.9608 0.4118], [0.9216 0.8863 0.4235]};
    triangle = update_triangle(last_pack.B,last_pack.A,last_pack.C,flank_B,flank_A,flank_C, ...
                               last_pack.XYB,last_pack.XYA,last_pack.XYC,colour); %#ok<NASGU>
end

% Define function that can get deployed unit
    function [B_flank_def, A_flank_def, C_flank_def,...
            alphaL_B_def, alphaL_A_def, alphaL_C_def,...
            flank1, flank2, flank3] = get_flank(p1, p2, p3)
        % Indecies of inner triangle
        f_flank = [19 20 21 22;
            27 28 29 30;
            23 24 25 26
            ];

        % Move and rotate each flank to match the current boundary
        p = [p1; p2; p3];
        base_idx = [22 30 26];
        edge_target = [p3-p1; p1-p2; p2-p3];
        flanks = cell(3,1);
        for ii = 1:3
            flk = tri0(f_flank(ii,:),:) + (p(ii,:) - tri0(base_idx(ii),:));
            v1 = flk(3,:) - flk(4,:);
            v2 = edge_target(ii,:);
            rot_ang = atan2(v1(1)*v2(2) - v1(2)*v2(1), dot(v1, v2));
            flanks{ii} = real((flk - p(ii,:)) * rotation(-rot_ang) + p(ii,:));
        end
        flank1 = flanks{1}; flank2 = flanks{2}; flank3 = flanks{3};

        B_flank_def = flank1(2,:);
        A_flank_def = flank2(2,:);
        C_flank_def = flank3(2,:);

        vB = flank1(2,:) - flank1(3,:);
        alphaL_B_def = atan2(vB(2),vB(1));
        vA = flank2(2,:) - flank2(3,:);
        alphaL_A_def = atan2(vA(2),vA(1));
        vC = flank3(2,:) - flank3(3,:);
        alphaL_C_def = atan2(vC(2),vC(1));
    end

%% Build final configuration
    function triangle = update_triangle(B,A,C,flankB,flankA,flankC,XYB,XYA,XYC,colour)

        % Inner triangle
        triangle = zeros(45,2);
        triangle(43,:) = A;
        triangle(44,:) = B;
        triangle(45,:) = C;

        % flanks
        triangle([19,20,21,22],:) = flankB;
        triangle([27 28 29 30],:) = flankA;
        triangle([23 24 25 26],:) = flankC;

        % plot unit
        hold on
        plot_triangle(triangle,colour);

        % filaments
        build_ligament(XYB, flankB, A, B, t, colour{3});
        build_ligament(XYA, flankA, C, A, t, colour{3});
        build_ligament(XYC, flankC, B, C, t, colour{3});
        hold off
        axis off
    end

    function [XY_inner, XY_outer] = build_ligament(XY_outer, flank, A, B, t, colour)
        % XY_outer : outer boundary polyline of ligament
        % flank    : flank points (use to determine direction on the flank side)
        % A, B     : attachment points on the inner triangle
        % t        : ligament thickness
        % colour   : fill color for patch

        if nargin < 6
            colour = [0.5 0.5 0.5];
        end

        n = size(XY_outer, 1);

        % start direction comes from flank
        dir_start = flank(1,:) - flank(2,:);
        dir_start = dir_start / norm(dir_start);

        % end direction comes from triangle attachment A-B
        dir_end = A - B;
        dir_end = dir_end / norm(dir_end);

        % create interpolated direction along the ligament
        dir_interp = zeros(n,2);
        for i = 1:n
            s = (i-1)/(n-1);
            dir = (1-s)*dir_start + s*dir_end;
            dir_interp(i,:) = dir / norm(dir);
        end

        % inner offset candidates
        XY_inner_1 = XY_outer + t * dir_interp;
        XY_inner_2 = XY_outer - t * dir_interp;

        % decide which side is actually the inner side (closer to triangle)
        target = 0.5*(A + B);
        d1 = norm(mean(XY_inner_1,1) - target);
        d2 = norm(mean(XY_inner_2,1) - target);

        if d1 < d2
            XY_inner = XY_inner_1;
        else
            XY_inner = XY_inner_2;
        end

        % create closed polygon
        verts = [XY_outer; flipud(XY_inner)];
        faces = 1:(2*n);

        % draw the solid ligament
        patch('Vertices', verts, 'Faces', faces, ...
            'FaceColor', colour, ...
            'FaceAlpha', 0.5, ...
            'EdgeColor', 'k', ...
            'LineWidth', 0.5,...
            'MarkerFaceColor', 'cyan');
    end

end
