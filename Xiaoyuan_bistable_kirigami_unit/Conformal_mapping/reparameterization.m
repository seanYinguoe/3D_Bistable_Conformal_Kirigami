function [v_initial, v_target, info] = reparameterization(v_out, f_out, vt_mesh, v_mesh, f_mesh, anisotropy_filter, opts)
%REPARAMETERIZATION Surface-constrained reparameterization with global bounds.
%   [v_initial, v_target, info] = reparameterization(v_out, f_out, ...
%       vt_mesh, v_mesh, f_mesh, anisotropy_filter, opts)
%
% Uses global per-triangle scale factors:
%   scale_facs = calculate_scale_facs(v_out, v_target, f_out)   % Mx3
%   lam_all = scale_facs(:)
%   ratio = max(lam_all) / min(lam_all)
% with monotone acceptance on ratio and angle violations.

if nargin < 6
    error('Not enough input arguments.');
end

if nargin < 7
    opts = struct();
end
if isempty(opts)
    opts = struct();
end

[ratio_max, alphaMax, opts, admissible] = parse_constraint_inputs(anisotropy_filter, opts);
opts = fill_default_opts(opts);

if size(v_out,2) == 2
    v_out = [v_out, zeros(size(v_out,1),1, 'like', v_out)];
elseif size(v_out,2) ~= 3
    error('v_out must be Nx2 or Nx3.');
end
if size(f_out,2) ~= 3
    error('f_out must be Mx3.');
end

if size(v_mesh,2) ~= 3
    error('v_mesh must be Ns x 3.');
end
if size(vt_mesh,2) > 2
    vt_mesh = vt_mesh(:,1:2);
end
if size(vt_mesh,1) ~= size(v_mesh,1)
    error('vt_mesh and v_mesh must have same number of vertices.');
end
if size(f_mesh,2) ~= 3
    error('f_mesh must be Mx3.');
end
if any(f_mesh(:) < 1) || any(f_mesh(:) > size(v_mesh,1))
    error('f_mesh contains invalid vertex indices.');
end
vt_mesh3 = [vt_mesh, zeros(size(vt_mesh,1),1, 'like', vt_mesh)];

% Register v_out(:,1:2) into vt_mesh(:,1:2) frame using centroid + uniform scale.
q_src = v_out(:,1:2);
uv = vt_mesh3(:,1:2);
c_src = mean(q_src, 1);
c_uv = mean(uv, 1);
r_src = q_src - c_src;
r_uv = uv - c_uv;
s_src = sqrt(mean(sum(r_src.^2, 2)));
s_uv = sqrt(mean(sum(r_uv.^2, 2)));
if s_src > eps(class(s_src))
    sim_s = s_uv / s_src;
else
    sim_s = 1;
end
q_reg = (q_src - c_src) * sim_s + c_uv;
vt_mesh = uv;

surfV = v_mesh;
surfF = f_mesh;

[v_initial, v_target] = initialize_pair_from_surface(v_out, q_reg, vt_mesh, surfV, surfF);

nV = size(v_target,1);
[~, vertNbrs] = build_edge_list_and_vertex_neighbors(f_out, nV);
[edgeI, edgeJ] = tri_edge_pairs(f_out);
L0_tri_edges = tri_edge_lengths(v_out, f_out);
L0_safe = max(L0_tri_edges, eps(class(L0_tri_edges)));
triData = precompute_surface_data(surfV, surfF);
admModel = build_admissible_model(anisotropy_filter, opts);

[v_target, anchorFace, faceNormalsAtPts] = project_points_to_surface( ...
    v_target, surfV, surfF, triData, zeros(size(v_target,1),1));

nIter = opts.nIter;
histRatio = zeros(nIter,1);
histMinLam = zeros(nIter,1);
histMaxLam = zeros(nIter,1);
histMaxAngle = zeros(nIter,1);
histObj = zeros(nIter,1);
histStepAccepted = zeros(nIter,1);
histRmsUpdate = zeros(nIter,1);
histAdmDistMean = zeros(nIter,1);
histAdmViolFrac = zeros(nIter,1);
histBoundViolSum = zeros(nIter,1);
histBoundViolCount = zeros(nIter,1);

status = "maxIter";

for it = 1:nIter
    state_old = evaluate_state(v_out, f_out, v_target, ratio_max, alphaMax, vertNbrs, admModel, opts);
    merit_old = constraint_merit(state_old.ratioViol, state_old.angleViol, state_old.phi_adm, state_old.boundViolSum);

    histRatio(it) = state_old.ratio;
    histMinLam(it) = state_old.minLam;
    histMaxLam(it) = state_old.maxLam;
    histMaxAngle(it) = state_old.maxAngle;
    histObj(it) = state_old.J;
    histAdmDistMean(it) = state_old.admDistMean;
    histAdmViolFrac(it) = state_old.admViolFrac;
    histBoundViolSum(it) = state_old.boundViolSum;
    histBoundViolCount(it) = state_old.boundViolCount;

    g = zeros(size(v_target), 'like', v_target);

    if state_old.phi_adm > 0
        g = g + admissible_edge_correction(v_target, edgeI, edgeJ, L0_safe, state_old.targetLam, state_old.distViol, opts.w_adm);
        g = g + local_angle_correction(v_target, f_out, state_old.triAngles, state_old.targetAlpha, state_old.distViol, opts.w_adm * opts.w_ang_local);
    end

    if state_old.angleViol > 0
        g = g + global_angle_correction(v_target, f_out, state_old.triAngles, alphaMax, opts.w_ang);
    end

    if state_old.ratioViol > 0
        lam_ref = exp(mean(log(max(state_old.lam_all, eps(class(state_old.lam_all))))));
        g = g + global_ratio_tail_correction(v_target, edgeI, edgeJ, L0_safe, state_old.scale_facs, lam_ref, opts.w_ratio, opts.qFrac);
    end

    if opts.w_smooth > 0
        g = g + opts.w_smooth * laplacian_smoothing(v_target, vertNbrs);
    end

    n = faceNormalsAtPts;
    dotgn = sum(g .* n, 2);
    g_t = g - dotgn .* n;
    % Remove rigid translation mode to avoid global drift.
    g_t = g_t - mean(g_t, 1);
    g_t = cap_row_norms(g_t, opts.maxDisp);

    accepted = false;
    stepCur = opts.step;
    stepUsed = 0;
    v_prev = v_target;

    for bt = 1:opts.maxBacktrack
        dv = stepCur * g_t;
        v_try = v_target + dv;
        [v_new, anchorFaceNew, faceNormalsNew] = project_points_to_surface(v_try, surfV, surfF, triData, anchorFace);

        state_new = evaluate_state(v_out, f_out, v_new, ratio_max, alphaMax, vertNbrs, admModel, opts);

        merit_new = constraint_merit(state_new.ratioViol, state_new.angleViol, state_new.phi_adm, state_new.boundViolSum);

        acceptStep = false;
        if merit_new < merit_old - opts.acceptTol
            acceptStep = true;
        elseif merit_new <= merit_old + opts.acceptTol && state_new.J <= state_old.J + opts.acceptTol ...
                && state_new.phi_adm <= state_old.phi_adm + opts.acceptTol ...
                && state_new.ratioViol <= state_old.ratioViol + opts.acceptTol ...
                && state_new.angleViol <= state_old.angleViol + opts.acceptTol ...
                && state_new.boundViolSum <= state_old.boundViolSum + opts.acceptTol ...
                && state_new.boundViolCount <= state_old.boundViolCount
            acceptStep = true;
        end

        if acceptStep
            accepted = true;
            v_target = v_new;
            anchorFace = anchorFaceNew;
            faceNormalsAtPts = faceNormalsNew;
            stepUsed = stepCur;
            histRatio(it) = state_new.ratio;
            histMinLam(it) = state_new.minLam;
            histMaxLam(it) = state_new.maxLam;
            histMaxAngle(it) = state_new.maxAngle;
            histObj(it) = state_new.J;
            histAdmDistMean(it) = state_new.admDistMean;
            histAdmViolFrac(it) = state_new.admViolFrac;
            histBoundViolSum(it) = state_new.boundViolSum;
            histBoundViolCount(it) = state_new.boundViolCount;
            histRmsUpdate(it) = sqrt(mean(sum((v_new - v_prev).^2, 2)));
            break;
        end

        stepCur = 0.5 * stepCur;
    end

    histStepAccepted(it) = stepUsed;

    if ~accepted
        status = "stalled";
        histRatio = histRatio(1:it);
        histMinLam = histMinLam(1:it);
        histMaxLam = histMaxLam(1:it);
        histMaxAngle = histMaxAngle(1:it);
        histObj = histObj(1:it);
        histStepAccepted = histStepAccepted(1:it);
        histRmsUpdate = histRmsUpdate(1:it);
        histAdmDistMean = histAdmDistMean(1:it);
        histAdmViolFrac = histAdmViolFrac(1:it);
        histBoundViolSum = histBoundViolSum(1:it);
        histBoundViolCount = histBoundViolCount(1:it);
        break;
    end

    if opts.verbose
        fprintf(['[reparameterization] iter %d/%d | ratio=%.6f | minLam=%.6f | ' ...
                 'maxLam=%.6f | maxAngleDeg=%.4f | admDist=%.4e | violFrac=%.4f | ' ...
                 'boundViol=%d(%.4e) | stepUsed=%.4e | J=%.6e\n'], ...
                 it, nIter, histRatio(it), histMinLam(it), histMaxLam(it), ...
                 rad2deg(histMaxAngle(it)), histAdmDistMean(it), histAdmViolFrac(it), ...
                 round(histBoundViolCount(it)), histBoundViolSum(it), stepUsed, histObj(it));
    end

    if state_new.ratioViol < opts.tol_ratio && state_new.angleViol < opts.tol_ang ...
            && state_new.admDistMean < opts.tol_adm ...
            && state_new.boundViolSum < opts.tol_bound ...
            && state_new.boundViolCount == 0 ...
            && histRmsUpdate(it) < opts.tol_rms
        status = "converged";
        histRatio = histRatio(1:it);
        histMinLam = histMinLam(1:it);
        histMaxLam = histMaxLam(1:it);
        histMaxAngle = histMaxAngle(1:it);
        histObj = histObj(1:it);
        histStepAccepted = histStepAccepted(1:it);
        histRmsUpdate = histRmsUpdate(1:it);
        histAdmDistMean = histAdmDistMean(1:it);
        histAdmViolFrac = histAdmViolFrac(1:it);
        histBoundViolSum = histBoundViolSum(1:it);
        histBoundViolCount = histBoundViolCount(1:it);
        break;
    end
end

scale_pre = calculate_scale_facs(v_initial, v_target, f_out);
lam_pre = scale_pre(:);
lam_pre = lam_pre(isfinite(lam_pre));
if isempty(lam_pre)
    rescale_facs = 1.0;
    rescale_feasible = true;
else
    if ~isempty(opts.lamHardMin) && ~isempty(opts.lamHardMax)
        minPre = min(lam_pre);
        maxPre = max(lam_pre);
        s_low = opts.lamHardMin / max(minPre, eps(class(minPre)));
        s_high = opts.lamHardMax / max(maxPre, eps(class(maxPre)));
        if s_low <= s_high
            s_pref = opts.lamMinTarget / max(minPre, eps(class(minPre)));
            rescale_facs = min(max(s_pref, s_low), s_high);
            rescale_feasible = true;
        else
            rescale_facs = 1.0;
            rescale_feasible = false;
        end
    else
        rescale_facs = opts.lamMinTarget / max(min(lam_pre), eps(class(lam_pre)));
        rescale_feasible = true;
    end
end
obj_surface = struct('v', v_mesh, 'f', struct('v', f_mesh));
[v_target, finalScaleFactors, info_rescale] = rescale_target_edges(v_initial, f_out, obj_surface, v_target, rescale_facs);
[~, finalMaxAnglePerTri, ~] = triangle_internal_angles(v_target, f_out, alphaMax);
lam_final = finalScaleFactors(:);
lam_final = lam_final(isfinite(lam_final));

% Keep v_initial on the same flattened plane as vt_mesh (z = 0).
v_initial(:,3) = 0;
% Match v_initial centroid to vt_mesh centroid on XY plane.
c_vt_xy = mean(vt_mesh, 1);
c_init_xy = mean(v_initial(:,1:2), 1);
v_initial(:,1:2) = v_initial(:,1:2) + (c_vt_xy - c_init_xy);

info = struct();
info.history.minLam = histMinLam;
info.history.maxLam = histMaxLam;
info.history.ratio = histRatio;
info.history.maxAngle = histMaxAngle;
info.history.rmsUpdate = histRmsUpdate;
info.history.stepAccepted = histStepAccepted;
info.history.obj = histObj;
info.history.admDistMean = histAdmDistMean;
info.history.admViolFrac = histAdmViolFrac;
info.history.boundViolSum = histBoundViolSum;
info.history.boundViolCount = histBoundViolCount;
info.finalScaleFactors = finalScaleFactors;
info.finalMinLam = min(lam_final);
info.finalMaxLam = max(lam_final);
info.finalRatio = info.finalMaxLam / max(info.finalMinLam, eps(class(info.finalMinLam)));
if ~isempty(opts.lamHardMin) && ~isempty(opts.lamHardMax)
    info.finalBoundViolCount = sum((lam_final < opts.lamHardMin) | (lam_final > opts.lamHardMax));
    info.finalBoundViolSum = sum(max(0, opts.lamHardMin - lam_final) + max(0, lam_final - opts.lamHardMax));
else
    info.finalBoundViolCount = 0;
    info.finalBoundViolSum = 0;
end
info.finalMaxAngle = max(finalMaxAnglePerTri);
info.v_out_transformed = v_out;
info.vt_mesh_transformed = vt_mesh3;
info.v_mesh_transformed = v_mesh;
info.q_registered = q_reg;
info.v_initial_xy_shift = c_vt_xy - c_init_xy;
info.rescale_facs = rescale_facs;
info.rescale = info_rescale;
info.rescale_feasible = rescale_feasible;
info.admissible = admissible;
if ~isempty(opts.lamHardMin) && ~isempty(opts.lamHardMax)
    info.hardBounds = [opts.lamHardMin, opts.lamHardMax];
else
    info.hardBounds = [];
end
info.status = status;
end

function [ratio_max, alphaMax, opts, admissible] = parse_constraint_inputs(anisotropy_filter, opts)
admissible = struct();
admissible.mode = "anisotropy_filter";
admissible.lamMin = NaN;
admissible.lamMax = NaN;
admissible.ratioMax = NaN;
admissible.alphaMax = NaN;

if ~(istable(anisotropy_filter) || isstruct(anisotropy_filter))
    error('anisotropy_filter must be a table or struct.');
end

[lamMinFilter, lamMaxFilter, alphaMax] = bounds_from_anisotropy_filter(anisotropy_filter);
if ~isfield(opts, 'ratioTarget') || isempty(opts.ratioTarget)
    opts.ratioTarget = 1.60 / 1.15;
end
if ~isfinite(opts.ratioTarget) || opts.ratioTarget <= 1
    error('opts.ratioTarget must be finite and > 1.');
end
ratio_max = opts.ratioTarget;

admissible.lamMin = lamMinFilter;
admissible.lamMax = lamMaxFilter;
admissible.ratioMax = ratio_max;
admissible.alphaMax = alphaMax;
admissible.hardLamMin = [];
admissible.hardLamMax = [];
end

function [lamMin, lamMax, alphaMax] = bounds_from_anisotropy_filter(anisotropy_filter)
eps_vals = pick_field_vector(anisotropy_filter, {'eps_bist', 'strain3'});
eps_vals = eps_vals(isfinite(eps_vals));
if isempty(eps_vals)
    error('anisotropy_filter must include finite eps_bist or strain3.');
end

lam_vals = 1 + eps_vals;
lam_vals = lam_vals(isfinite(lam_vals) & lam_vals > 0);
if isempty(lam_vals)
    error('anisotropy_filter strain field must produce positive scale factors.');
end

lamMin = min(lam_vals);
lamMax = max(lam_vals);
if lamMax < lamMin
    tmp = lamMin;
    lamMin = lamMax;
    lamMax = tmp;
end

a1 = pick_field_vector(anisotropy_filter, {'a1'});
a2 = pick_field_vector(anisotropy_filter, {'a2'});
a3 = pick_field_vector(anisotropy_filter, {'a3'});
ang = [a1(:); a2(:); a3(:)];
ang = ang(isfinite(ang));
if isempty(ang)
    error('anisotropy_filter must include finite a1, a2, a3.');
end
alphaMax = max(ang);
end

function v = pick_field_vector(S, names)
v = [];
for k = 1:numel(names)
    name = names{k};
    if istable(S)
        if any(strcmp(S.Properties.VariableNames, name))
            v = S.(name);
            return;
        end
    elseif isstruct(S)
        if isfield(S, name)
            v = S.(name);
            return;
        end
    else
        error('Input must be table or struct.');
    end
end
error('Missing required field(s): %s', strjoin(names, ', '));
end

function opts = fill_default_opts(opts)
if ~isfield(opts, 'nIter') || isempty(opts.nIter), opts.nIter = 60; end
if ~isfield(opts, 'step') || isempty(opts.step), opts.step = 0.1; end
if ~isfield(opts, 'w_adm') || isempty(opts.w_adm), opts.w_adm = 4.0; end
if ~isfield(opts, 'w_ratio') || isempty(opts.w_ratio), opts.w_ratio = 0.4; end
if ~isfield(opts, 'w_bound') || isempty(opts.w_bound), opts.w_bound = 8.0; end
if ~isfield(opts, 'w_ang') || isempty(opts.w_ang), opts.w_ang = 0.5; end
if ~isfield(opts, 'w_ang_local') || isempty(opts.w_ang_local), opts.w_ang_local = 0.35; end
if ~isfield(opts, 'w_smooth') || isempty(opts.w_smooth), opts.w_smooth = 0.02; end
if ~isfield(opts, 'verbose') || isempty(opts.verbose), opts.verbose = false; end
if ~isfield(opts, 'qFrac') || isempty(opts.qFrac), opts.qFrac = 0.10; end
if ~isfield(opts, 'maxBacktrack') || isempty(opts.maxBacktrack), opts.maxBacktrack = 10; end
if ~isfield(opts, 'tol_ratio') || isempty(opts.tol_ratio), opts.tol_ratio = 1e-4; end
if ~isfield(opts, 'tol_ang') || isempty(opts.tol_ang), opts.tol_ang = 1e-6; end
if ~isfield(opts, 'tol_adm') || isempty(opts.tol_adm), opts.tol_adm = 1e-3; end
if ~isfield(opts, 'tol_bound') || isempty(opts.tol_bound), opts.tol_bound = 1e-9; end
if ~isfield(opts, 'tol_rms') || isempty(opts.tol_rms), opts.tol_rms = 1e-8; end
if ~isfield(opts, 'maxDisp') || isempty(opts.maxDisp), opts.maxDisp = inf; end
if ~isfield(opts, 'admTol') || isempty(opts.admTol), opts.admTol = 0.08; end
if ~isfield(opts, 'acceptTol') || isempty(opts.acceptTol), opts.acceptTol = 1e-10; end
if ~isfield(opts, 'nnChunk') || isempty(opts.nnChunk), opts.nnChunk = 2048; end
if ~isfield(opts, 'maxFilter') || isempty(opts.maxFilter), opts.maxFilter = 30000; end
if ~isfield(opts, 'ratioTarget') || isempty(opts.ratioTarget), opts.ratioTarget = 1.60 / 1.15; end
if ~isfield(opts, 'lamHardMin'), opts.lamHardMin = []; end
if ~isfield(opts, 'lamHardMax'), opts.lamHardMax = []; end
if ~isfield(opts, 'lamMinTarget') || isempty(opts.lamMinTarget), opts.lamMinTarget = 1.15; end
end

function state = evaluate_state(v_out, f_out, v_target, ratio_max, alphaMax, vertNbrs, admModel, opts)
scale_facs = calculate_scale_facs(v_out, v_target, f_out);
lam_all = scale_facs(:);
lam_all = lam_all(isfinite(lam_all));
if isempty(lam_all)
    error('No finite scale factors available in reparameterization.');
end
minLam = min(lam_all);
maxLam = max(lam_all);
ratio = maxLam / max(minLam, eps(class(minLam)));
ratioViol = max(0, log(ratio / ratio_max));
if ~isempty(opts.lamHardMin) && ~isempty(opts.lamHardMax)
    boundLow = max(0, opts.lamHardMin - lam_all);
    boundHigh = max(0, lam_all - opts.lamHardMax);
    boundViol = boundLow + boundHigh;
    boundViolSum = sum(boundViol);
    boundViolCount = sum(boundViol > 0);
else
    boundViol = zeros(size(lam_all), 'like', lam_all);
    boundViolSum = 0;
    boundViolCount = 0;
end

[triAngles, maxAnglePerTri, ~] = triangle_internal_angles(v_target, f_out, alphaMax);
maxAngle = max(maxAnglePerTri);
angleViol = max(0, maxAngle - alphaMax);

X = descriptor_from_scale_facs(scale_facs, triAngles);
[idxNN, distRaw] = nearest_admissible_rows(X, admModel, opts.nnChunk);
distViol = max(0, distRaw - opts.admTol);
phi_adm = mean(distViol.^2);

phi_smooth = smoothing_energy(v_target, vertNbrs);
J = opts.w_adm * phi_adm + opts.w_ratio * ratioViol^2 + opts.w_bound * (boundViolSum^2 / max(numel(lam_all),1)) ...
    + opts.w_ang * angleViol^2 + opts.w_smooth * phi_smooth;

state = struct();
state.J = J;
state.scale_facs = scale_facs;
state.lam_all = lam_all;
state.minLam = minLam;
state.maxLam = maxLam;
state.ratio = ratio;
state.ratioViol = ratioViol;
state.boundViolSum = boundViolSum;
state.boundViolCount = boundViolCount;
state.boundViol = boundViol;
state.triAngles = triAngles;
state.maxAnglePerTri = maxAnglePerTri;
state.maxAngle = maxAngle;
state.angleViol = angleViol;
state.admDistRaw = distRaw;
state.distViol = distViol;
state.phi_adm = phi_adm;
state.admDistMean = mean(distViol);
state.admViolFrac = mean(distViol > 0);
state.idxNN = idxNN;
state.targetLam = admModel.lamRef(idxNN);
state.targetAlpha = admModel.alphaRef(idxNN);
end

function X = descriptor_from_scale_facs(scale_facs, triAnglesFallback)
strainFallback = mean(scale_facs, 2) - 1;
try
    anisotropy_level = scale_facs_to_angles(scale_facs);
    a1 = pick_field_vector(anisotropy_level, {'a1'});
    a2 = pick_field_vector(anisotropy_level, {'a2'});
    a3 = pick_field_vector(anisotropy_level, {'a3'});
    strain3 = pick_field_vector(anisotropy_level, {'strain3', 'eps_bist'});
    angSort = sort([a1(:), a2(:), a3(:)], 2, 'ascend');
    strain = strain3(:);
    bad = ~all(isfinite(angSort),2) | ~isfinite(strain);
    if any(bad)
        angFb = sort(triAnglesFallback, 2, 'ascend');
        angSort(bad,:) = angFb(bad,:);
        strain(bad) = strainFallback(bad);
    end
    X = [angSort, strain];
catch
    angSort = sort(triAnglesFallback, 2, 'ascend');
    X = [angSort, strainFallback];
end
end

function g = admissible_edge_correction(v_target, edgeI, edgeJ, L0_safe, targetLam, triWeight, w)
triW = max(0, triWeight(:));
if ~any(triW > 0)
    g = zeros(size(v_target), 'like', v_target);
    return;
end

iIdx = edgeI(:);
jIdx = edgeJ(:);
xi = v_target(iIdx,:);
xj = v_target(jIdx,:);
dij = xi - xj;
Lcur = sqrt(sum(dij.^2, 2));
Lsafe = max(Lcur, eps(class(Lcur)));
dir = dij ./ Lsafe;

targetL = repelem(targetLam(:), 3, 1) .* L0_safe(:);
delta = targetL - Lcur;
wEdge = w .* repelem(triW, 3, 1);
corr = (wEdge .* delta) .* dir;
g = scatter_add_edges(size(v_target,1), iIdx, jIdx, corr);
end

function g = local_angle_correction(v_target, f_out, triAngles, targetAlpha, triWeight, w)
g = zeros(size(v_target), 'like', v_target);
maxTri = max(triAngles, [], 2);
excess = max(0, maxTri - targetAlpha(:));
wTri = w .* max(0, triWeight(:)) .* excess;
violTri = find(wTri > 0);
for tt = 1:numel(violTri)
    t = violTri(tt);
    tri = f_out(t,:);
    [~, kLoc] = max(triAngles(t,:));
    if kLoc == 1
        iV = tri(1); jV = tri(2); kV = tri(3);
    elseif kLoc == 2
        iV = tri(2); jV = tri(3); kV = tri(1);
    else
        iV = tri(3); jV = tri(1); kV = tri(2);
    end
    midJK = 0.5 * (v_target(jV,:) + v_target(kV,:));
    dirI = midJK - v_target(iV,:);
    nrm = norm(dirI);
    if nrm <= 0
        continue;
    end
    dirI = dirI / nrm;
    wt = wTri(t);
    g(iV,:) = g(iV,:) + wt * dirI;
    g(jV,:) = g(jV,:) - 0.5 * wt * dirI;
    g(kV,:) = g(kV,:) - 0.5 * wt * dirI;
end
end

function g = global_angle_correction(v_target, f_out, triAngles, alphaMax, w_ang)
g = zeros(size(v_target), 'like', v_target);
maxTri = max(triAngles, [], 2);
violTri = find(maxTri > alphaMax);
for tt = 1:numel(violTri)
    t = violTri(tt);
    tri = f_out(t,:);
    [~, kLoc] = max(triAngles(t,:));
    if kLoc == 1
        iV = tri(1); jV = tri(2); kV = tri(3);
    elseif kLoc == 2
        iV = tri(2); jV = tri(3); kV = tri(1);
    else
        iV = tri(3); jV = tri(1); kV = tri(2);
    end
    excess = maxTri(t) - alphaMax;
    midJK = 0.5 * (v_target(jV,:) + v_target(kV,:));
    dirI = midJK - v_target(iV,:);
    nrm = norm(dirI);
    if nrm <= 0
        continue;
    end
    dirI = dirI / nrm;
    wt = w_ang * 2 * excess;
    g(iV,:) = g(iV,:) + wt * dirI;
    g(jV,:) = g(jV,:) - 0.5 * wt * dirI;
    g(kV,:) = g(kV,:) - 0.5 * wt * dirI;
end
end

function g = global_ratio_tail_correction(v_target, edgeI, edgeJ, L0_safe, scale_facs, lam_ref, w_ratio, qFrac)
g = zeros(size(v_target), 'like', v_target);
validMask = isfinite(scale_facs);
lam_vec = scale_facs(validMask);
if isempty(lam_vec)
    return;
end
validLin = find(validMask);
qCount = max(1, ceil(qFrac * numel(lam_vec)));
[~, ord] = sort(lam_vec, 'ascend');
lowLin = validLin(ord(1:qCount));
highLin = validLin(ord(max(1, numel(ord)-qCount+1):end));
g = g + ratio_edge_correction(v_target, edgeI, edgeJ, L0_safe, scale_facs, highLin, lam_ref, w_ratio, -1);
g = g + ratio_edge_correction(v_target, edgeI, edgeJ, L0_safe, scale_facs, lowLin,  lam_ref, w_ratio, +1);
end

function model = build_admissible_model(anisotropy_filter, opts)
a1 = pick_field_vector(anisotropy_filter, {'a1'});
a2 = pick_field_vector(anisotropy_filter, {'a2'});
a3 = pick_field_vector(anisotropy_filter, {'a3'});
epsVals = pick_field_vector(anisotropy_filter, {'eps_bist', 'strain3'});

ang = sort([a1(:), a2(:), a3(:)], 2, 'ascend');
strain = epsVals(:);
lamRef = 1 + strain;
alphaRef = max(ang, [], 2);

X = [ang, strain];
good = all(isfinite(X), 2) & isfinite(lamRef) & (lamRef > 0);
X = X(good, :);
lamRef = lamRef(good);
alphaRef = alphaRef(good);
if isempty(X)
    error('anisotropy_filter has no finite admissible rows.');
end

if size(X,1) > opts.maxFilter
    idx = round(linspace(1, size(X,1), opts.maxFilter));
    X = X(idx,:);
    lamRef = lamRef(idx);
    alphaRef = alphaRef(idx);
end

mu = mean(X, 1);
sigma = std(X, 0, 1);
sigma(sigma < 1e-12) = 1;
XN = (X - mu) ./ sigma;

model = struct();
model.mu = mu;
model.sigma = sigma;
model.XN = XN;
model.lamRef = lamRef;
model.alphaRef = alphaRef;
end

function [idxNN, dist] = nearest_admissible_rows(X, model, chunkSize)
XN = (X - model.mu) ./ model.sigma;
Y = model.XN;
nX = size(XN,1);
nY = size(Y,1);
idxNN = ones(nX,1);
best = inf(nX,1);

Y2 = sum(Y.^2, 2)';
for s = 1:chunkSize:nX
    e = min(s + chunkSize - 1, nX);
    XC = XN(s:e, :);
    X2 = sum(XC.^2, 2);
    D2 = X2 + Y2 - 2 * (XC * Y');
    D2 = max(D2, 0);
    [b, id] = min(D2, [], 2);
    best(s:e) = b;
    idxNN(s:e) = id;
end
dist = sqrt(best);
end

function g = ratio_edge_correction(v_target, edgeI, edgeJ, L0_safe, scale_facs, linIdx, lam_ref, w_ratio, signMode)
g = zeros(size(v_target), 'like', v_target);
if isempty(linIdx)
    return;
end

[triIdx, edgeIdx] = ind2sub(size(scale_facs), linIdx);
iIdx = edgeI(linIdx);
jIdx = edgeJ(linIdx);

xi = v_target(iIdx,:);
xj = v_target(jIdx,:);
dij = xi - xj;
Lcur = sqrt(sum(dij.^2, 2));
Lsafe = max(Lcur, eps(class(Lcur)));
dir = dij ./ Lsafe;

L0sel = L0_safe(sub2ind(size(L0_safe), triIdx, edgeIdx));
targetL = lam_ref .* L0sel;

if signMode < 0
    delta = max(0, Lcur - targetL);
    corr = -w_ratio * delta .* dir;
else
    delta = max(0, targetL - Lcur);
    corr = w_ratio * delta .* dir;
end

g = scatter_add_edges(size(v_target,1), iIdx, jIdx, corr);
end

function g = scatter_add_edges(nV, iIdx, jIdx, corr)
g = zeros(nV, 3, 'like', corr);
for d = 1:3
    g(:,d) = g(:,d) + accumarray(iIdx, corr(:,d), [nV,1], @sum, 0);
    g(:,d) = g(:,d) - accumarray(jIdx, corr(:,d), [nV,1], @sum, 0);
end
end

function [J, ratioViol, angleViol, phi_smooth] = objective_value(v_target, ratio, maxAngle, ratio_max, alphaMax, vertNbrs, opts)
ratioViol = max(0, log(ratio / ratio_max));
angleViol = max(0, maxAngle - alphaMax);
phi_smooth = smoothing_energy(v_target, vertNbrs);
J = opts.w_ratio * ratioViol^2 + opts.w_ang * angleViol^2 + opts.w_smooth * phi_smooth;
end

function m = constraint_merit(ratioViol, angleViol, admViol, boundViol)
if nargin < 3, admViol = 0; end
if nargin < 4, boundViol = 0; end
m = max([ratioViol, angleViol, admViol, boundViol]);
end

function phi_smooth = smoothing_energy(v_target, vertNbrs)
phi_smooth = 0;
nCount = 0;
for vi = 1:size(v_target,1)
    nb = vertNbrs{vi};
    if ~isempty(nb)
        d = v_target(vi,:) - mean(v_target(nb,:), 1);
        phi_smooth = phi_smooth + sum(d.^2);
        nCount = nCount + 1;
    end
end
if nCount > 0
    phi_smooth = phi_smooth / nCount;
end
end

function lap = laplacian_smoothing(v_target, vertNbrs)
lap = zeros(size(v_target), 'like', v_target);
for vi = 1:size(v_target,1)
    nb = vertNbrs{vi};
    if ~isempty(nb)
        lap(vi,:) = mean(v_target(nb,:), 1) - v_target(vi,:);
    end
end
end

function dv = cap_row_norms(dv, maxDisp)
if ~isfinite(maxDisp) || maxDisp <= 0
    return;
end
nrm = sqrt(sum(dv.^2, 2));
scale = ones(size(nrm), 'like', nrm);
mask = nrm > maxDisp;
scale(mask) = maxDisp ./ nrm(mask);
dv = dv .* scale;
end

function [edgeI, edgeJ] = tri_edge_pairs(f_out)
edgeI = [f_out(:,1), f_out(:,2), f_out(:,3)];
edgeJ = [f_out(:,2), f_out(:,3), f_out(:,1)];
end

function Ltri = tri_edge_lengths(V, f_out)
edge12 = sqrt(sum((V(f_out(:,2),:) - V(f_out(:,1),:)).^2, 2));
edge23 = sqrt(sum((V(f_out(:,3),:) - V(f_out(:,2),:)).^2, 2));
edge31 = sqrt(sum((V(f_out(:,1),:) - V(f_out(:,3),:)).^2, 2));
Ltri = [edge12, edge23, edge31];
end

function [E, vertNbrs] = build_edge_list_and_vertex_neighbors(f_out, nV)
Eall = [f_out(:,[1 2]); f_out(:,[2 3]); f_out(:,[3 1])];
Eall = sort(Eall, 2);
E = unique(Eall, 'rows');

vertNbrs = cell(nV,1);
for k = 1:size(E,1)
    i = E(k,1);
    j = E(k,2);
    vertNbrs{i}(end+1) = j;
    vertNbrs{j}(end+1) = i;
end
for i = 1:nV
    vertNbrs{i} = unique(vertNbrs{i});
end
end

function [ang, maxAngPerTri, maxViol] = triangle_internal_angles(V, f_out, alphaMax)
A = V(f_out(:,1),:);
B = V(f_out(:,2),:);
C = V(f_out(:,3),:);

AB = B - A; AC = C - A;
BA = A - B; BC = C - B;
CA = A - C; CB = B - C;

angA = safe_acos_row(dot_rows(AB, AC) ./ max(row_norm(AB).*row_norm(AC), eps(class(A))));
angB = safe_acos_row(dot_rows(BA, BC) ./ max(row_norm(BA).*row_norm(BC), eps(class(A))));
angC = safe_acos_row(dot_rows(CA, CB) ./ max(row_norm(CA).*row_norm(CB), eps(class(A))));

ang = [angA, angB, angC];
maxAngPerTri = max(ang, [], 2);
maxViol = max(max(maxAngPerTri - alphaMax, 0));
end

function d = dot_rows(X, Y)
d = sum(X .* Y, 2);
end

function n = row_norm(X)
n = sqrt(sum(X.^2, 2));
end

function a = safe_acos_row(x)
x = min(1, max(-1, x));
a = acos(x);
end

function triData = precompute_surface_data(surfV, surfF)
F = size(surfF,1);

p1 = surfV(surfF(:,1),:);
p2 = surfV(surfF(:,2),:);
p3 = surfV(surfF(:,3),:);

fn = cross(p2-p1, p3-p1, 2);
fnNorm = sqrt(sum(fn.^2,2));
valid = fnNorm > 0;
fn(valid,:) = fn(valid,:) ./ fnNorm(valid);
fn(~valid,:) = repmat([0 0 1], sum(~valid), 1);
fc = (p1 + p2 + p3) / 3;

nV = size(surfV,1);
v2f = cell(nV,1);
for f = 1:F
    vv = surfF(f,:);
    v2f{vv(1)}(end+1) = f;
    v2f{vv(2)}(end+1) = f;
    v2f{vv(3)}(end+1) = f;
end

E = [surfF(:,[1 2]); surfF(:,[2 3]); surfF(:,[3 1])];
fid = [(1:F)'; (1:F)'; (1:F)'];
E = sort(E,2);
[~, ~, ic] = unique(E, 'rows');
faceNbrs = cell(F,1);
for k = 1:numel(ic)
    grp = find(ic == ic(k));
    if numel(grp) > 1
        faces = unique(fid(grp));
        for i = 1:numel(faces)
            f0 = faces(i);
            faceNbrs{f0} = unique([faceNbrs{f0}; faces(:)]);
        end
    end
end
for f = 1:F
    faceNbrs{f}(faceNbrs{f} == f) = [];
end

triData = struct();
triData.faceNormals = fn;
triData.faceCenters = fc;
triData.v2f = v2f;
triData.faceNbrs = faceNbrs;
end

function VN = compute_vertex_normals(surfV, surfF)
NV = size(surfV,1);
VN = zeros(NV,3, 'like', surfV);

p1 = surfV(surfF(:,1),:);
p2 = surfV(surfF(:,2),:);
p3 = surfV(surfF(:,3),:);
FN = cross(p2-p1, p3-p1, 2);

for k = 1:size(surfF,1)
    tri = surfF(k,:);
    VN(tri(1),:) = VN(tri(1),:) + FN(k,:);
    VN(tri(2),:) = VN(tri(2),:) + FN(k,:);
    VN(tri(3),:) = VN(tri(3),:) + FN(k,:);
end

nrm = sqrt(sum(VN.^2,2));
idx = nrm > 0;
VN(idx,:) = VN(idx,:) ./ nrm(idx);
VN(~idx,:) = repmat([0 0 1], sum(~idx), 1);
end

function [v_proj, anchorFace, n_at_point] = project_points_to_surface(vq, surfV, surfF, triData, anchorFace)
nq = size(vq,1);
v_proj = zeros(nq,3, 'like', vq);
n_at_point = zeros(nq,3, 'like', vq);

if nargin < 5 || isempty(anchorFace)
    anchorFace = zeros(nq,1);
end

fc = triData.faceCenters;
faceNormals = triData.faceNormals;
faceNbrs = triData.faceNbrs;

for i = 1:nq
    q = vq(i,:);

    if anchorFace(i) >= 1 && anchorFace(i) <= size(surfF,1)
        cand = [anchorFace(i); faceNbrs{anchorFace(i)}(:)];
        cand = unique(cand);

        if numel(cand) < 6
            ring2 = cand;
            for ii = 1:numel(cand)
                ring2 = [ring2; faceNbrs{cand(ii)}(:)]; %#ok<AGROW>
            end
            cand = unique(ring2);
        end
    else
        d2 = sum((fc - q).^2, 2);
        [~, idx] = sort(d2, 'ascend');
        K = min(30, numel(idx));
        cand = idx(1:K);
    end

    [cp, fBest, ~] = closest_point_on_mesh_faces(q, surfV, surfF, cand);

    if fBest == 0
        [cp, fBest, ~] = closest_point_on_mesh_faces(q, surfV, surfF, (1:size(surfF,1)).');
    end

    v_proj(i,:) = cp;
    anchorFace(i) = fBest;
    n_at_point(i,:) = faceNormals(fBest,:);
end

nrm = sqrt(sum(n_at_point.^2,2));
bad = nrm < eps(class(vq));
if any(bad)
    VN = compute_vertex_normals(surfV, surfF);
    d2 = pdist2_local(v_proj(bad,:), surfV);
    [~, idv] = min(d2, [], 2);
    n_at_point(bad,:) = VN(idv,:);
end

nrm = sqrt(sum(n_at_point.^2,2));
good = nrm > 0;
n_at_point(good,:) = n_at_point(good,:) ./ nrm(good);
n_at_point(~good,:) = repmat([0 0 1], sum(~good), 1);
end

function [cpBest, fBest, d2Best] = closest_point_on_mesh_faces(q, surfV, surfF, faceIdx)
cpBest = [0 0 0];
fBest = 0;
d2Best = inf;

for kk = 1:numel(faceIdx)
    f = faceIdx(kk);
    tri = surfF(f,:);
    a = surfV(tri(1),:);
    b = surfV(tri(2),:);
    c = surfV(tri(3),:);

    cp = closest_point_on_triangle(q, a, b, c);
    d2 = sum((q - cp).^2);

    if d2 < d2Best
        d2Best = d2;
        cpBest = cp;
        fBest = f;
    end
end
end

function cp = closest_point_on_triangle(p, a, b, c)
ab = b - a;
ac = c - a;
ap = p - a;
d1 = dot(ab, ap);
d2 = dot(ac, ap);
if d1 <= 0 && d2 <= 0
    cp = a;
    return;
end

bp = p - b;
d3 = dot(ab, bp);
d4 = dot(ac, bp);
if d3 >= 0 && d4 <= d3
    cp = b;
    return;
end

vc = d1*d4 - d3*d2;
if vc <= 0 && d1 >= 0 && d3 <= 0
    v = d1 / (d1 - d3);
    cp = a + v * ab;
    return;
end

cpv = p - c;
d5 = dot(ab, cpv);
d6 = dot(ac, cpv);
if d6 >= 0 && d5 <= d6
    cp = c;
    return;
end

vb = d5*d2 - d1*d6;
if vb <= 0 && d2 >= 0 && d6 <= 0
    w = d2 / (d2 - d6);
    cp = a + w * ac;
    return;
end

va = d3*d6 - d5*d4;
if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
    bc = c - b;
    w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    cp = b + w * bc;
    return;
end

n = cross(ab, ac);
n2 = dot(n, n);
if n2 <= eps
    cp = a;
    return;
end
cp = p - (dot(p - a, n) / n2) * n;
end

function D2 = pdist2_local(A, B)
AA = sum(A.^2, 2);
BB = sum(B.^2, 2)';
D2 = max(AA + BB - 2*(A*B'), 0);
end

function [v_initial, v_target] = initialize_pair_from_surface(v_out, q_reg, uv, surfV, surfF)
q = q_reg;

TR = triangulation(surfF, uv);
[ti, bc] = pointLocation(TR, q);

v_target = zeros(size(q,1), 3, 'like', surfV);
v_initial = [v_out(:,1:2), zeros(size(v_out,1),1, 'like', v_out)];
inside = ~isnan(ti);
if any(inside)
    faces = surfF(ti(inside), :);
    v_target(inside,:) = ...
        bc(inside,1) .* surfV(faces(:,1),:) + ...
        bc(inside,2) .* surfV(faces(:,2),:) + ...
        bc(inside,3) .* surfV(faces(:,3),:);
end

outside = find(~inside);
if isempty(outside)
    return;
end

uv_centers = (uv(surfF(:,1),:) + uv(surfF(:,2),:) + uv(surfF(:,3),:)) / 3;
for kk = 1:numel(outside)
    idx = outside(kk);
    d2 = sum((uv_centers - q(idx,:)).^2, 2);
    [~, fIdx] = min(d2);
    tri_uv = uv(surfF(fIdx,:), :);
    tri_uv3 = [tri_uv, zeros(3,1, 'like', tri_uv)];
    p3 = [q(idx,:), 0];
    bc_loc = cart2barycentric(tri_uv3, p3);
    v_target(idx,:) = bc_loc(1) * surfV(surfF(fIdx,1),:) + ...
                      bc_loc(2) * surfV(surfF(fIdx,2),:) + ...
                      bc_loc(3) * surfV(surfF(fIdx,3),:);
end
end
