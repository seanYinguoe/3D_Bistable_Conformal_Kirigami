    function [E_total,alpha] = deform_triangle_anisotropic(q1,q2,q3,edgeLen,l1,l4,beta,t,nD,N,colour)
    % Inputs:
    %   q1, q2, q3   - Unit node coordinates (1x3 vectors)
    %   edgeLen      - Original triangle edge length
    %   l1, l4, beta   - l1: length of flanks l4:thickness of flanks
    %   beta:tilting angle
    %   t            - Thickness of filaments
    %   i_out        - Orientation flag (0 for upwards, 1 for downwards)
    %   nD           - Number of steps
    %
    % Output:
    %   triangle_new - Deformed triangle coordinates
    %   E_total      - Deployed energy
    
    %% Helper
    wrap  = @(th) atan2(sin(th), cos(th));
    
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
    % plot_triangle(tri0);
    E_total  = nan(1, nD);
    theta_all = nan(1, nD);
    theta_prev = -pi + beta;             % start with undefomred one
    
    %% input initial guess for optimisation
    K = N-1;  % number of torsional spring
    DeltaTot = wrap(theta_prev + pi - beta); % total angle difference
    phi_prevB = (DeltaTot / K) * ones(K,1); % initial rotational angle of torsional springs
    e_prevB   = zeros(N,1); % initial length change in ligaments
    phi_prevA = (DeltaTot / K) * ones(K,1);
    e_prevA   = zeros(N,1);
    phi_prevC = (DeltaTot / K) * ones(K,1);
    e_prevC   = zeros(N,1);
    G0 = [-sqrt(3)/6*edgeLen, -1/2*edgeLen];
    xG0 = G0(1);
    yG0 = G0(2);
    x_prev = [phi_prevB; e_prevB; phi_prevA; e_prevA; phi_prevC; e_prevC; theta_prev; xG0; yG0]; % initial guess
    %x_prev = [phi_prevA; e_prevA; phi_prevB; e_prevB;theta_prev]; % initial guess
    alpha = zeros(nD,1);
    %% Loop over deltas
    for k = 1:nD*0.75
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
    
        %G0 = 1/3*(p1 + p2 + p3);
        %xG0 = G0(1);
        %yG0 = G0(2);
    
        % update params for optimisation
        params = struct('t',t,'E',Emod,'b',b,'beta',beta,'B_flank',B_flank,...
            'A_flank',A_flank,'C_flank',C_flank,'r_vertex',r_vertex,'l2',l2,...
            'l3',l3,'alphaL_B',alphaL_B,'alphaL_A',alphaL_A,'alphaL_C',alphaL_C,...
            'xG0',xG0,'yG0',yG0);
    
        [energy, pack] = energy_lig_anisotropic(x_prev, N, params); % N=10 segments
        % Save the reuslt
        theta_all(k) = pack.theta;
        E_total(k)   = energy;               % 3 ligaments total
        theta_prev   = pack.theta;
        phi_prevB = pack.phiB;
        e_prevB = pack.eB;
        phi_prevA = pack.phiA;
        e_prevA = pack.eA;
        phi_prevC = pack.phiC;
        e_prevC = pack.eC;
        xG_prev = pack.xG;
        yG_prev = pack.yG;
        x_prev = [phi_prevB; e_prevB; phi_prevA; e_prevA; phi_prevC; e_prevC; theta_prev; xG_prev; yG_prev];
        xG0 = xG_prev;
        yG0 = yG_prev;
        %x_prev = [phi_prevA; e_prevA;phi_prevB; e_prevB; theta_prev];
        last_pack = pack;
    end
    
    
    %% Plot energy curve
    figure('Color','w');
    hold on; box on;
    plot(alpha, E_total, '-', 'Color',[0.85 0.33 0.10],'LineWidth', 1);
    xlabel('Deployment', 'Interpreter','tex', ...
           'FontSize',20);
    ylabel('Strain Energy(N/mm^2)', 'Interpreter','tex', ...
           'FontSize',20);
    set(gca, 'FontName','Times New Roman','FontSize',20);
    legend('Location','northwest','Box','off', 'Fontsize',18);
    grid off;
    axis square;

    %% Plot configuration
    %figure()
    %triangle = update_triangle(last_pack.B,last_pack.A,last_pack.C,flank_B,flank_A,flank_C, last_pack.XYB,last_pack.XYA, last_pack.XYC,colour);
    
    % Define function that can get deployed unit
        function [B_flank_def, A_flank_def, C_flank_def,...
                alphaL_B_def, alphaL_A_def, alphaL_C_def,...
                flank1, flank2, flank3] = get_flank(p1, p2, p3)
            % Indecies of inner triangle
            f_flank = [19 20 21 22;
                27 28 29 30;
                23 24 25 26
                ];
    
            % Move the flanks to fit the outer boundary(rigid conditions)
            flank1_t = p1 - tri0(22,:);
            flank1 = tri0(f_flank(1,:),:) + flank1_t;
            flank2_t = p2 - tri0(30,:);
            flank2 = tri0(f_flank(2,:),:) + flank2_t;
            flank3_t = p3 - tri0(26,:);
            flank3 = tri0(f_flank(3,:),:) + flank3_t;
    
            % Rotate the flanks to fit the outer boundary
            vector1 = flank1(3,:) - flank1(4,:);
            vector2 = p3 - p1;
            flank1_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));
            vector1 = flank2(3,:) - flank2(4,:);
            vector2 = p1 - p2;
            flank2_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));
            vector1 = flank3(3,:) - flank3(4,:);
            vector2 = p2 - p3;
            flank3_r = atan2(vector1(1)*vector2(2) - vector1(2)*vector2(1), dot(vector1, vector2));
    
            flank1 = real((flank1-p1)*rotation(-flank1_r) + p1);
            flank2 = real((flank2-p2)*rotation(-flank2_r) + p2);
            flank3 = real((flank3-p3)*rotation(-flank3_r) + p3);
    
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
            plot_triangle(triangle,colour)
    
            % filaments
            build_ligament(XYB, flankB, A, B, t, colour{3})
            build_ligament(XYA, flankA, C, A, t, colour{3})
            build_ligament(XYC, flankC, B, C, t, colour{3})
            
            % plot(XYB(:,1),  XYB(:,2),'LineWidth', 1.2, 'MarkerSize', 4, ...
            %     'DisplayName', 'deformed');
            % plot(XYA(:,1),  XYA(:,2),'LineWidth', 1.2, 'MarkerSize', 4, ...
            %     'DisplayName', 'deformed');
            % plot(XYC(:,1),  XYC(:,2),'LineWidth', 1.2, 'MarkerSize', 4, ...
            %     'DisplayName', 'deformed');
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
    
         % compute tangent direction along XY_outer
        dXY = diff(XY_outer, 1, 1);
        segLen = vecnorm(dXY, 2, 2);
        T = dXY ./ segLen;
    
        T_pts = zeros(n,2);
        T_pts(1,:) = T(1,:);
        T_pts(end,:) = T(end,:);
        if n > 2
            T_pts(2:end-1,:) = 0.5 * (T(1:end-1,:) + T(2:end,:));
        end
    
        % curve normal
        N_curve = [-T_pts(:,2), T_pts(:,1)];
        N_curve = N_curve ./ vecnorm(N_curve,2,2);
    
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

