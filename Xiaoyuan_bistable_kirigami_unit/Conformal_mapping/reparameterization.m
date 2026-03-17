function [v_initial, v_target, max_ang, v_target0] = reparameterization(v_out, f_out, vt_mesh, v_mesh, f_mesh)
%REPARAMETERIZATION Surface-constrained reparameterization with robust objective minimization.
%   [v_initial, v_target, max_ang, v_target0] = reparameterization( ...
%       v_out, f_out, vt_mesh, v_mesh, f_mesh)
%
% Inputs
%   v_out   : Nx2 or Nx3 planar grid vertices
%   f_out   : Mx3 grid connectivity
%   vt_mesh : Ns x2 UV coordinates of target mesh
%   v_mesh  : Ns x3 target surface vertices
%   f_mesh  : Kx3 target surface connectivity
%
% Outputs
%   v_initial : Nx3 planar grid
%   v_target  : Nx3 optimized curved grid constrained on target surface
%   max_ang   : maximum internal angle (radians)
%   v_target0 : barycentric initial curved grid (before optimization)

%% ------------------------------- checks -------------------------------
if nargin < 5
    error('reparameterization requires five inputs.');
end
if size(f_out,2) ~= 3 || size(f_mesh,2) ~= 3
    error('f_out and f_mesh must be Mx3 triangle connectivity.');
end
if size(v_mesh,2) ~= 3
    error('v_mesh must be Ns x3.');
end
if size(vt_mesh,2) > 2
    vt_mesh = vt_mesh(:,1:2);
end
if size(vt_mesh,2) ~= 2 || size(vt_mesh,1) ~= size(v_mesh,1)
    error('vt_mesh must be Ns x2 and match v_mesh rows.');
end

if size(v_out,2) == 2
    v_out = [v_out, zeros(size(v_out,1),1, 'like', v_out)];
elseif size(v_out,2) ~= 3
    error('v_out must be Nx2 or Nx3.');
end

%% ------------------------- constants / weights ------------------------
lam_min = 1.15;
lam_max = 1.65;
ratio_lim = lam_max / lam_min;

% Objective weights: anisotropy is primary.
weights = struct();
weights.w_iso_mean = 1.00;
weights.w_iso_max = 1.50;     % strong worst-unit control
weights.w_ratio = 80.0;       % strong global ratio enforcement
weights.w_bounds = 20.0;      % soft absolute bound helper
weights.w_smooth = 0.01;      % mild regularization
weights.alpha_cmp = 1.00;     % anisotropy comparison score: mean + alpha*max

nIter = 60;
lineMax = 12;
step0 = 0.30;
objTol = 1e-12;
gradTol = 1e-9;
fdRel = 2e-3;

%% ----------------------- planar registration in UV --------------------
q_src = v_out(:,1:2);
q_reg = register_planar_to_uv(q_src, vt_mesh);

%% -------------------- barycentric initialization ----------------------
% Keep v_target0 from barycentric interpolation (old method).
[v_initial, v_target0] = initialize_pair_from_surface_old(v_out, q_reg, vt_mesh, v_mesh, f_mesh);

triData = precompute_surface_data(v_mesh, f_mesh);
[v_target, anchorFace, normals] = project_points_to_surface(v_target0, v_mesh, f_mesh, triData, zeros(size(v_target0,1),1)); %#ok<ASGLU>

L0 = tri_edge_lengths(v_initial, f_out);
L0_safe = max(L0, eps(class(L0)));
[~, vertNbrs] = build_edge_list_and_vertex_neighbors(f_out, size(v_target,1));

state0 = evaluate_state(v_target0, f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, weights);
state = state0;

L0v = L0_safe(:);
L0v = L0v(isfinite(L0v) & L0v > 0);
if isempty(L0v)
    h_fd = 1e-4;
else
    h_fd = fdRel * median(L0v);
end
h_fd = max(h_fd, 1e-6);

%% --------------------- projected gradient descent ---------------------
for it = 1:nIter
    [grad_tan, normals, anchorFace] = tangent_fd_gradient(v_target, anchorFace, normals, ...
        f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, weights, h_fd, ...
        v_mesh, f_mesh, triData);

    gnorm = sqrt(sum(grad_tan(:).^2));
    if ~isfinite(gnorm) || gnorm < gradTol
        break;
    end

    dir = -grad_tan;
    step = step0;
    accepted = false;

    for bt = 1:lineMax
        v_try = v_target + step * dir;
        [v_try, anchor_try, normals_try] = project_points_to_surface(v_try, v_mesh, f_mesh, triData, anchorFace);
        state_try = evaluate_state(v_try, f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, weights);

        ratio_ok_progress = (state_try.ratioViol <= state.ratioViol + 1e-14) || (state_try.ratioViol <= 1e-14);
        if isfinite(state_try.obj) && (state_try.obj < state.obj - objTol) && ratio_ok_progress
            v_target = v_try;
            state = state_try;
            anchorFace = anchor_try;
            normals = normals_try;
            accepted = true;
            break;
        end

        step = 0.5 * step;
    end

    if ~accepted
        break;
    end
end

%% ------------------------- final evaluation ---------------------------
[v_target, ~, ~] = project_points_to_surface(v_target, v_mesh, f_mesh, triData, anchorFace);
state_final = evaluate_state(v_target, f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, weights);

% Fallback rule: if anisotropy score is worse than barycentric, return v_target0.
if state_final.anisoScore > state0.anisoScore + objTol
    v_target = v_target0;
    state_final = state0;
end

max_ang = state_final.maxAngle;

end

%% =====================================================================
%% Objective and anisotropy metrics

function state = evaluate_state(v_target, f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, w)
Lcur = tri_edge_lengths(v_target, f_out);
scale_facs = Lcur ./ L0_safe;
triAngles = triangle_internal_angles(v_target, f_out);

% Per-triangle geometric isotropy (edge equality in deformed metric).
E_iso_len_tri = tri_isotropy_metric(Lcur);

% Per-triangle scale isotropy (edge-scale equality).
E_iso_scale_tri = tri_isotropy_metric(scale_facs);

% Combined anisotropy per triangle.
E_iso_tri = 0.5 * (E_iso_len_tri + E_iso_scale_tri);
E_iso_mean = mean(E_iso_tri, 'omitnan');
E_iso_max = max(E_iso_tri, [], 'omitnan');

% Global scale-ratio constraint.
lam = scale_facs(:);
lam = lam(isfinite(lam));
if isempty(lam)
    minLam = inf;
    maxLam = inf;
    ratio = inf;
else
    minLam = min(lam);
    maxLam = max(lam);
    ratio = maxLam / max(minLam, eps(class(maxLam)));
end
ratioViol = max(0, ratio - ratio_lim);
E_ratio = ratioViol^2;

% Mild helper penalty to keep absolute extremes close to [lam_min, lam_max].
violMin = max(0, lam_min - minLam);
violMax = max(0, maxLam - lam_max);
E_bounds = violMin^2 + violMax^2;

% Mild smoothness.
lap = laplacian_residual(v_target, vertNbrs);
E_smooth = mean(sum(lap.^2, 2), 'omitnan');

% Primary objective.
obj = w.w_iso_mean * E_iso_mean + ...
      w.w_iso_max * E_iso_max + ...
      w.w_ratio * E_ratio + ...
      w.w_bounds * E_bounds + ...
      w.w_smooth * E_smooth;

if ~isfinite(obj)
    obj = inf;
end

% Scalar anisotropy score used for final fallback comparison.
anisoScore = E_iso_mean + w.alpha_cmp * E_iso_max;

state = struct();
state.obj = obj;
state.anisoScore = anisoScore;
state.E_iso_mean = E_iso_mean;
state.E_iso_max = E_iso_max;
state.E_ratio = E_ratio;
state.E_bounds = E_bounds;
state.E_smooth = E_smooth;
state.ratio = ratio;
state.ratioViol = ratioViol;
state.minScale = minLam;
state.maxScale = maxLam;
state.scale_facs = scale_facs;
state.maxAngle = max(triAngles(:));
end

function E_tri = tri_isotropy_metric(T)
% T is Mx3 with triangle edge descriptors.
a = T(:,1);
b = T(:,2);
c = T(:,3);
num = (a-b).^2 + (b-c).^2 + (c-a).^2;
den = a.^2 + b.^2 + c.^2 + eps(class(T));
E_tri = num ./ den;
end

%% =====================================================================
%% Finite-difference tangent gradient

function [grad_tan, normals, anchorFace] = tangent_fd_gradient(v_target, anchorFace, normals, ...
    f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, w, h_fd, ...
    v_mesh, f_mesh, triData)

nV = size(v_target,1);
grad_tan = zeros(size(v_target), 'like', v_target);

[v_target, anchorFace, normals] = project_points_to_surface(v_target, v_mesh, f_mesh, triData, anchorFace);
[t1, t2] = tangent_basis_from_normals(normals);

for i = 1:nV
    [fp1, fm1, anchorFace] = directional_objective(v_target, i, t1(i,:), h_fd, anchorFace, ...
        f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, w, v_mesh, f_mesh, triData);

    [fp2, fm2, anchorFace] = directional_objective(v_target, i, t2(i,:), h_fd, anchorFace, ...
        f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, w, v_mesh, f_mesh, triData);

    d1 = (fp1 - fm1) / (2*h_fd);
    d2 = (fp2 - fm2) / (2*h_fd);

    if isfinite(d1)
        grad_tan(i,:) = grad_tan(i,:) + d1 * t1(i,:);
    end
    if isfinite(d2)
        grad_tan(i,:) = grad_tan(i,:) + d2 * t2(i,:);
    end
end
end

function [f_plus, f_minus, anchorFace] = directional_objective(v_base, vi, dir, h, anchorFace, ...
    f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, w, v_mesh, f_mesh, triData)

v_p = v_base;
v_m = v_base;

[p_plus, face_plus] = project_single_point(v_base(vi,:) + h * dir, anchorFace(vi), v_mesh, f_mesh, triData);
[p_minus, face_minus] = project_single_point(v_base(vi,:) - h * dir, anchorFace(vi), v_mesh, f_mesh, triData);

v_p(vi,:) = p_plus;
v_m(vi,:) = p_minus;

s_p = evaluate_state(v_p, f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, w);
s_m = evaluate_state(v_m, f_out, L0_safe, lam_min, lam_max, ratio_lim, vertNbrs, w);

f_plus = s_p.obj;
f_minus = s_m.obj;

if face_plus > 0
    anchorFace(vi) = face_plus;
elseif face_minus > 0
    anchorFace(vi) = face_minus;
end
end

function [t1, t2] = tangent_basis_from_normals(normals)
nV = size(normals,1);
t1 = zeros(nV,3, 'like', normals);
t2 = zeros(nV,3, 'like', normals);

for i = 1:nV
    n = normals(i,:);
    nn = norm(n);
    if nn <= 0
        n = [0 0 1];
    else
        n = n / nn;
    end

    if abs(n(3)) < 0.9
        r = [0 0 1];
    else
        r = [1 0 0];
    end

    a = cross(n, r);
    na = norm(a);
    if na <= 0
        r = [0 1 0];
        a = cross(n, r);
        na = max(norm(a), eps);
    end
    a = a / na;

    b = cross(n, a);
    nb = max(norm(b), eps);
    b = b / nb;

    t1(i,:) = a;
    t2(i,:) = b;
end
end

%% =====================================================================
%% Basic geometry and mesh helpers

function q_reg = register_planar_to_uv(q_src, uv)
% Similarity registration (translation + uniform scale) to target UV frame.
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
end

function lap = laplacian_residual(v_target, vertNbrs)
lap = zeros(size(v_target), 'like', v_target);
for vi = 1:size(v_target,1)
    nb = vertNbrs{vi};
    if ~isempty(nb)
        lap(vi,:) = v_target(vi,:) - mean(v_target(nb,:), 1);
    end
end
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
    vertNbrs{i}(end+1) = j; %#ok<AGROW>
    vertNbrs{j}(end+1) = i; %#ok<AGROW>
end
for i = 1:nV
    vertNbrs{i} = unique(vertNbrs{i});
end
end

function ang = triangle_internal_angles(V, f_out)
A = V(f_out(:,1),:);
B = V(f_out(:,2),:);
C = V(f_out(:,3),:);

AB = B - A; AC = C - A;
BA = A - B; BC = C - B;
CA = A - C; CB = B - C;

denA = max(row_norm(AB).*row_norm(AC), eps(class(A)));
denB = max(row_norm(BA).*row_norm(BC), eps(class(A)));
denC = max(row_norm(CA).*row_norm(CB), eps(class(A)));

angA = safe_acos_row(dot_rows(AB, AC) ./ denA);
angB = safe_acos_row(dot_rows(BA, BC) ./ denB);
angC = safe_acos_row(dot_rows(CA, CB) ./ denC);
ang = [angA, angB, angC];
end

function d = dot_rows(X, Y)
d = sum(X .* Y, 2);
end

function n = row_norm(X)
n = sqrt(sum(X.^2, 2));
end

function y = safe_acos_row(x)
x = min(max(x, -1), 1);
y = acos(x);
end

%% =====================================================================
%% Surface projection helpers

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
triData.faceNbrs = faceNbrs;
end

function [v_proj, anchorFace, n_at_point] = project_points_to_surface(vq, surfV, surfF, triData, anchorFace)
nq = size(vq,1);
v_proj = zeros(nq,3, 'like', vq);
n_at_point = zeros(nq,3, 'like', vq);

if nargin < 5 || isempty(anchorFace)
    anchorFace = zeros(nq,1);
end

for i = 1:nq
    [cp, fBest] = project_single_point(vq(i,:), anchorFace(i), surfV, surfF, triData);
    v_proj(i,:) = cp;
    anchorFace(i) = fBest;
    n_at_point(i,:) = triData.faceNormals(fBest,:);
end
end

function [cp, fBest] = project_single_point(q, faceHint, surfV, surfF, triData)
fc = triData.faceCenters;
faceNbrs = triData.faceNbrs;

if faceHint >= 1 && faceHint <= size(surfF,1)
    cand = [faceHint; faceNbrs{faceHint}(:)];
    cand = unique(cand);
else
    d2 = sum((fc - q).^2, 2);
    [~, idx] = sort(d2, 'ascend');
    cand = idx(1:min(40, numel(idx)));
end

[cp, fBest, ~] = closest_point_on_mesh_faces(q, surfV, surfF, cand);
if fBest == 0
    [cp, fBest, ~] = closest_point_on_mesh_faces(q, surfV, surfF, (1:size(surfF,1)).');
end
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

%% =====================================================================
%% Initialization helpers

function [v_initial, v_target] = initialize_pair_from_surface_old(v_out, q_reg, uv, surfV, surfF)
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
