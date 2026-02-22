function [strain_bist, bistability, info] = bistability_analysis(l1, l4, beta, t, edgeLen, q1, q2, q3, do_plot, L_confirm, peakRiseMin, dropMin, bistabilityMin)
% Compute bistability metrics using incremental warm-start deployment.
%
% Inputs:
%   l1, l4, beta : Geometric parameters
%   t            : Filament width
%   edgeLen      : Original outer-triangle edge length
%   q1, q2, q3   : Target deformed boundary vertices (1x3)
%   do_plot      : true/false, plot energy-deployment curve (optional)
%   L_confirm    : slope persistence look-ahead length (optional, default 3)
%   peakRiseMin  : minimum normalized rise to accept peak (optional, default 0.02)
%   dropMin      : minimum normalized drop to accept valley (optional, default 0.02)
%   bistabilityMin : minimum bistability threshold (optional, default 0.1)
%
% Outputs:
%   strain_bist  : Strain at first meaningful local minimum after first meaningful local maximum
%   bistability  : (E_max - E_min) / E_max
%   info         : Diagnostics struct with fields:
%                  status, k_stop, alpha_stop,
%                  E_history, alpha_history, strain_history,
%                  max_idx, min_idx, bistability_raw

%% Settings
N = 8;
dalpha = 1/200;
alpha_max = 1.0;

if nargin < 9 || isempty(do_plot)
    do_plot = false;
end
if nargin < 10 || isempty(L_confirm)
    L_confirm = 3;
end
if nargin < 11 || isempty(peakRiseMin)
    peakRiseMin = 0.02;
end
if nargin < 12 || isempty(dropMin)
    dropMin = 0.02;
end
if nargin < 13 || isempty(bistabilityMin)
    bistabilityMin = 0.04;
end

do_plot = logical(do_plot);
L_confirm = max(1, round(L_confirm));

% Detection settings
slope_eps_rel = 1e-4;           % relative slope deadband
persist_frac = 0.67;            % "mostly" sign persistence threshold

% Post-peak reliability gate (single metric)
maxFlipFrac = 0.25;

% Peak neighborhood validity gate (cheap local smoothness test)
% Suggested defaults:
%   mPeak=4       -> 4 slopes on each side of candidate peak
%   asymFrac=0.15 -> reject if right-side descent is too flat vs left ascent
%   jumpMult=8    -> reject severe one-step slope kinks at peak
%   useJumpGate   -> enable/disable jump gate
mPeak = 4;
asymFrac = 0.15;
jumpMult = 8;
useJumpGate = true;

%% Geometry constants (same model as deform_triangle_anisotropic)
wrap  = @(th) atan2(sin(th), cos(th));
pack_x = @(phiB,eB,phiA,eA,phiC,eC,theta,xG,yG) [phiB; eB; phiA; eA; phiC; eC; theta; xG; yG];

l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(pi/3 - beta);
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) - cos(pi/3 - beta) .* l4;
l5 = ((sqrt(3)/2) .* l4 + sin(beta) .* l6) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 ...
    - l2 .* ((sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)));
r_vertex = (sqrt(3)/3) * l3;

Emod = 4.3e11;
b = 1.0;

%% Outer-edge interpolation endpoints
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

%% Reference unit and warm-start initial state
[tri0, ~] = deform_triangle_isotropic(0, edgeLen, l1, l4, beta, t, N);

theta_prev = -pi + beta;
K = N - 1;
DeltaTot = wrap(theta_prev + pi - beta);
phi0 = (DeltaTot / K) * ones(K,1);
e0 = zeros(N,1);
G0 = [-sqrt(3)/6*edgeLen, -1/2*edgeLen];
xG0 = G0(1);
yG0 = G0(2);
x_prev = pack_x(phi0, e0, phi0, e0, phi0, e0, theta_prev, xG0, yG0);

%% Params template (update dynamic fields only each step)
params = struct('t',t,'E',Emod,'b',b,'beta',beta, ...
    'B_flank',[0 0],'A_flank',[0 0],'C_flank',[0 0], ...
    'r_vertex',r_vertex,'l2',l2,'l3',l3, ...
    'alphaL_B',0,'alphaL_A',0,'alphaL_C',0, ...
    'xG0',xG0,'yG0',yG0);

%% Histories and outputs
maxSteps = floor(alpha_max / dalpha) + 1;
alpha_hist = nan(maxSteps,1);
E_hist = nan(maxSteps,1);
strain_hist = nan(maxSteps,1);

% === NEW: store per-step configuration for bistable-shape plotting ===
pack_hist = cell(maxSteps,1);
flankB_hist = cell(maxSteps,1);
flankA_hist = cell(maxSteps,1);
flankC_hist = cell(maxSteps,1);

strain_scale = norm(q1 - q2) / edgeLen - 1;

strain_bist = NaN;
bistability = NaN;
bistability_raw = NaN;
status = 'no_transition';

found_max = false;
max_idx = NaN;
min_idx = NaN;
last_checked_i = 3; % candidate index i checked once k >= i+L_confirm

% === NEW: store bistable configuration ===
pack_bist = [];
flankB_bist = [];
flankA_bist = [];
flankC_bist = [];
alpha_bist = NaN;
stop_now = false;

% Diagnostics for post-peak quality gate
flipFrac_postpeak = NaN;
postpeak_ok = false;

% Diagnostics for peak-neighborhood gate
peak_magL = NaN;
peak_magR = NaN;
peak_ratioRtoL = NaN;
peak_jumpRel = NaN;
peak_ok = false;

%% Incremental deployment loop
for k = 1:maxSteps
    alpha_k = min((k - 1) * dalpha, alpha_max);

    [energy, pack, flankB_step, flankA_step, flankC_step] = step_deploy(alpha_k);

    alpha_hist(k) = alpha_k;
    E_hist(k) = energy;
    strain_hist(k) = alpha_k * strain_scale;

    % Cache step state for potential later retrieval (e.g., min_idx < current k)
    pack_hist{k} = pack;
    flankB_hist{k} = flankB_step;
    flankA_hist{k} = flankA_step;
    flankC_hist{k} = flankC_step;

    % Warm-start continuity for next step
    xG0 = pack.xG;
    yG0 = pack.yG;
    x_prev = pack_x(pack.phiB, pack.eB, pack.phiA, pack.eA, pack.phiC, pack.eC, ...
        pack.theta, xG0, yG0);

    % Persistent extrema detection with look-ahead confirmation
    if k >= (L_confirm + 4)
        i_max_check = k - L_confirm - 1; % ensures post slopes i:i+L_confirm exist
        processed_upto = i_max_check;
        for i = last_checked_i:i_max_check
            if i < 3
                continue;
            end

            [pre_slopes, post_slopes, slope_eps] = slope_windows(alpha_hist(1:k), E_hist(1:k), i);

            if ~found_max
                pre_ok = mostly_sign(pre_slopes, +1, slope_eps, persist_frac);
                post_ok = mostly_sign(post_slopes, -1, slope_eps, persist_frac);
                if pre_ok && post_ok
                    Ei = E_hist(i);
                    Emin_pre = min(E_hist(1:i));
                    rise_ratio = (Ei - Emin_pre) / max(abs(Ei), eps);
                    if rise_ratio >= peakRiseMin
                        [peak_ok, peak_magL, peak_magR, peak_ratioRtoL, peak_jumpRel] = ...
                            peak_neighborhood_gate(i, k);
                        if ~peak_ok
                            status = 'invalid_peak';
                            strain_bist = NaN;
                            bistability = NaN;
                            max_idx = NaN;
                            min_idx = NaN;
                            bistability_raw = NaN;
                            stop_now = true;
                            processed_upto = i;
                            break;
                        end
                        found_max = true;
                        max_idx = i;
                        processed_upto = i;
                        break;
                    end
                end
            else
                if i <= max_idx
                    continue;
                end
                pre_ok = mostly_sign(pre_slopes, -1, slope_eps, persist_frac);
                post_ok = mostly_sign(post_slopes, +1, slope_eps, persist_frac);
                if pre_ok && post_ok
                    Ei = E_hist(i);
                    Emax = E_hist(max_idx);
                    drop_ratio = (Emax - Ei) / max(abs(Emax), eps);
                    if drop_ratio >= dropMin
                        % Candidate valley found: validate post-peak smoothness/robustness
                        cand_min_idx = i;
                        cand_bistability_raw = drop_ratio;
                        processed_upto = i;

                        [flipFrac_postpeak, postpeak_ok] = postpeak_flip_gate(max_idx, cand_min_idx);

                        if ~postpeak_ok
                            status = 'monostable_osc';
                            strain_bist = NaN;
                            bistability = NaN;
                            min_idx = NaN;
                            bistability_raw = NaN;
                            stop_now = true; % break early to save time
                        elseif cand_bistability_raw < bistabilityMin
                            min_idx = cand_min_idx;
                            bistability_raw = cand_bistability_raw;
                            strain_bist = NaN;
                            bistability = NaN;
                            status = 'monostable_weak';
                        else
                            min_idx = cand_min_idx;
                            bistability_raw = cand_bistability_raw;
                            strain_bist = strain_hist(min_idx);
                            bistability = bistability_raw;
                            status = 'bistable';

                            % === NEW: store bistable configuration ===
                            pack_bist = pack_hist{min_idx};
                            flankB_bist = flankB_hist{min_idx};
                            flankA_bist = flankA_hist{min_idx};
                            flankC_bist = flankC_hist{min_idx};
                            alpha_bist = alpha_hist(min_idx);
                        end
                        break;
                    end
                end
            end
        end
        last_checked_i = processed_upto + 1;

        if stop_now || ~isnan(min_idx)
            break;
        end
    end

    if alpha_k >= alpha_max
        break;
    end
end

% Trim histories to evaluated samples
k_stop = find(~isnan(alpha_hist), 1, 'last');
alpha_hist = alpha_hist(1:k_stop);
E_hist = E_hist(1:k_stop);
strain_hist = strain_hist(1:k_stop);

if strcmp(status, 'no_transition')
    strain_bist = NaN;
    bistability = NaN;
end

info = struct();
info.status = status;
info.k_stop = k_stop;
info.alpha_stop = alpha_hist(end);
info.E_history = E_hist;
info.alpha_history = alpha_hist;
info.strain_history = strain_hist;
info.max_idx = max_idx;
info.min_idx = min_idx;
info.bistability_raw = bistability_raw;
info.flipFrac_postpeak = flipFrac_postpeak;
info.postpeak_ok = postpeak_ok;
info.peak_magL = peak_magL;
info.peak_magR = peak_magR;
info.peak_ratioRtoL = peak_ratioRtoL;
info.peak_jumpRel = peak_jumpRel;
info.peak_ok = peak_ok;
if strcmp(status, 'bistable')
    info.alpha_bist = alpha_bist;
else
    info.alpha_bist = NaN;
end
info.has_shape_plot = false;

if do_plot
    % Same plotting style as deform_triangle_anisotropic
    figure('Color','w');
    hold on; box on;
    plot(alpha_hist, E_hist, '-', 'Color',[0.85 0.33 0.10], 'LineWidth', 1);

    % Mark only final bistable pair (not weak/monostable classifications)
    if strcmp(status, 'bistable') && ~isnan(max_idx) && max_idx <= k_stop
        plot(alpha_hist(max_idx), E_hist(max_idx), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        text(alpha_hist(max_idx), E_hist(max_idx), ' max energy', ...
            'VerticalAlignment', 'bottom', 'FontName', 'Times New Roman', 'FontSize', 14);
    end
    if strcmp(status, 'bistable') && ~isnan(min_idx) && min_idx <= k_stop
        plot(alpha_hist(min_idx), E_hist(min_idx), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        text(alpha_hist(min_idx), E_hist(min_idx), ' bistable energy', ...
            'VerticalAlignment', 'top', 'FontName', 'Times New Roman', 'FontSize', 14);
    else
        % Optional status annotation for non-bistable cases
        text(alpha_hist(end), E_hist(end), [' status: ' status], ...
            'VerticalAlignment', 'bottom', 'FontName', 'Times New Roman', 'FontSize', 12);
    end

    xlabel('Deployment', 'Interpreter','tex', 'FontSize',20);
    ylabel('Strain Energy(N/mm^2)', 'Interpreter','tex', 'FontSize',20);
    set(gca, 'FontName','Times New Roman','FontSize',20);
    legend('Location','northwest','Box','off', 'Fontsize',18);
    grid off;
    axis square;

    % === NEW: plot bistable configuration ===
    if strcmp(status, 'bistable') && ~isempty(pack_bist)
        figure('Color','w');
        colour = {'white', [0.9216 0.8863 0.4235], [0.7059 0.9608 0.4118], [0.9216 0.8863 0.4235]};
        update_triangle(pack_bist.B, pack_bist.A, pack_bist.C, ...
            flankB_bist, flankA_bist, flankC_bist, ...
            pack_bist.XYB, pack_bist.XYA, pack_bist.XYC, colour);
        title(sprintf('Bistable configuration at \\alpha = %.3f', alpha_bist), ...
            'FontName', 'Times New Roman', 'FontSize', 16);
        info.has_shape_plot = true;
    end
end

    function [energy, pack, flankB, flankA, flankC] = step_deploy(alpha_k)
        % Interpolate boundary nodes
        p1 = (1 - alpha_k) * p1_orig + alpha_k * p1_def;
        p2 = (1 - alpha_k) * p2_orig + alpha_k * p2_def;
        p3 = (1 - alpha_k) * p3_orig + alpha_k * p3_def;

        [B_flank, A_flank, C_flank, alphaL_B, alphaL_A, alphaL_C, flankB, flankA, flankC] = get_flank(p1, p2, p3);

        % Update dynamic params only
        params.B_flank = B_flank;
        params.A_flank = A_flank;
        params.C_flank = C_flank;
        params.alphaL_B = alphaL_B;
        params.alphaL_A = alphaL_A;
        params.alphaL_C = alphaL_C;
        params.xG0 = xG0;
        params.yG0 = yG0;

        [energy, pack] = energy_lig_anisotropic(x_prev, N, params);
    end

    function [B_flank_def, A_flank_def, C_flank_def, alphaL_B_def, alphaL_A_def, alphaL_C_def, flank1, flank2, flank3] = get_flank(p1, p2, p3)
        f_flank = [19 20 21 22;
                   27 28 29 30;
                   23 24 25 26];

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

    function [pre_slopes, post_slopes, slope_eps] = slope_windows(alpha_v, E_v, i)
        % slope j corresponds to segment (j -> j+1)
        dE = diff(E_v);
        da = diff(alpha_v);
        s = dE ./ da;

        pre_slopes = s(i-2:i-1);
        post_slopes = s(i:min(i+L_confirm, numel(s)));

        Escale = max(abs(E_v));
        slope_eps = max(1e-12, slope_eps_rel * max(Escale, 1));
    end

    function ok = mostly_sign(vals, sign_target, tol, frac)
        if isempty(vals)
            ok = false;
            return;
        end
        if sign_target > 0
            n_good = sum(vals > tol);
        else
            n_good = sum(vals < -tol);
        end
        ok = (n_good >= ceil(frac * numel(vals)));
    end

    function [flipFrac, isOk] = postpeak_flip_gate(iMax, iMin)
        % Post-peak gate using only derivative sign flips between peak and valley.
        if iMin <= iMax + 2
            flipFrac = inf;
            isOk = false;
            return;
        end

        idx = iMax:iMin;
        aSeg = alpha_hist(idx);
        eSeg = E_hist(idx);
        da = diff(aSeg);
        de = diff(eSeg);
        s = de ./ da;

        Escale = max(abs(eSeg));
        s_eps = max(1e-12, slope_eps_rel * max(Escale, 1));
        sig = zeros(size(s));
        sig(s > s_eps) = 1;
        sig(s < -s_eps) = -1;
        sigNZ = sig(sig ~= 0);

        if numel(sigNZ) < 2
            flipFrac = 0;
        else
            flipFrac = sum(diff(sigNZ) ~= 0) / (numel(sigNZ) - 1);
        end

        isOk = (flipFrac <= maxFlipFrac);
    end

    function [ok, magL, magR, ratioRtoL, jumpRel] = peak_neighborhood_gate(iPeak, kNow)
        % Validate candidate peak with local slopes around iPeak only.
        magL = NaN;
        magR = NaN;
        ratioRtoL = NaN;
        jumpRel = NaN;
        ok = false;

        % Need mPeak slopes on left and right of iPeak.
        if iPeak < (mPeak + 1)
            return;
        end
        if (iPeak + mPeak) > kNow
            return;
        end

        aSeg = alpha_hist(1:kNow);
        eSeg = E_hist(1:kNow);
        s = diff(eSeg) ./ diff(aSeg);

        sL = s(iPeak-mPeak:iPeak-1);
        sR = s(iPeak:iPeak+mPeak-1);

        Escale = max(abs(eSeg));
        s_eps = max(1e-12, slope_eps_rel * max(Escale, 1));

        medL = median(sL);
        medR = median(sR);
        magL = median(abs(sL));
        magR = median(abs(sR));
        ratioRtoL = magR / max(magL, eps);

        trend_ok = (medL > s_eps) && (medR < -s_eps);
        asym_ok = (magR >= asymFrac * magL);

        if useJumpGate
            baseMag = median(abs([sL; sR]));
            jumpRel = abs(s(iPeak) - s(iPeak-1)) / max(baseMag, eps);
            jump_ok = (jumpRel <= jumpMult);
        else
            jump_ok = true;
        end

        ok = trend_ok && asym_ok && jump_ok;
    end

    function triangle = update_triangle(B,A,C,flankB,flankA,flankC,XYB,XYA,XYC,colour)
        triangle = zeros(45,2);
        triangle(43,:) = A;
        triangle(44,:) = B;
        triangle(45,:) = C;

        triangle([19,20,21,22],:) = flankB;
        triangle([27 28 29 30],:) = flankA;
        triangle([23 24 25 26],:) = flankC;

        hold on
        plot_triangle(triangle,colour)

        build_ligament(XYB, flankB, A, B, t, colour{3})
        build_ligament(XYA, flankA, C, A, t, colour{3})
        build_ligament(XYC, flankC, B, C, t, colour{3})
        hold off
        axis off
    end

    function [XY_inner, XY_outer] = build_ligament(XY_outer, flank, A, B, t_lig, colour)
        if nargin < 6
            colour = [0.5 0.5 0.5];
        end

        n = size(XY_outer, 1);

        dir_start = flank(1,:) - flank(2,:);
        dir_start = dir_start / norm(dir_start);

        dir_end = A - B;
        dir_end = dir_end / norm(dir_end);

        dir_interp = zeros(n,2);
        for i_pt = 1:n
            s_loc = (i_pt-1)/(n-1);
            dir = (1-s_loc)*dir_start + s_loc*dir_end;
            dir_interp(i_pt,:) = dir / norm(dir);
        end

        XY_inner_1 = XY_outer + t_lig * dir_interp;
        XY_inner_2 = XY_outer - t_lig * dir_interp;

        target = 0.5*(A + B);
        d1 = norm(mean(XY_inner_1,1) - target);
        d2 = norm(mean(XY_inner_2,1) - target);

        if d1 < d2
            XY_inner = XY_inner_1;
        else
            XY_inner = XY_inner_2;
        end

        verts = [XY_outer; flipud(XY_inner)];
        faces = 1:(2*n);

        patch('Vertices', verts, 'Faces', faces, ...
            'FaceColor', colour, ...
            'FaceAlpha', 0.5, ...
            'EdgeColor', 'k', ...
            'LineWidth', 0.5,...
            'MarkerFaceColor', 'cyan');
    end

end
