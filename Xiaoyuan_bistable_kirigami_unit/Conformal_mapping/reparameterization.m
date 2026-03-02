function [v_target, info] = reparameterization(v_out, f_out, obj_2D, v_target, ratio_max, alphaMax, opts)
%REPARAMETERIZATION Surface-constrained reparameterization with scale-ratio/angle bounds.
%   [v_target, info] = reparameterization(v_out, f_out, obj_2D, v_target, ...
%       ratio_max, alphaMax, opts)
%
% Inputs (kept as requested):
%   v_out    : Nx2 or Nx3 flattened configuration
%   f_out    : Mx3 triangles on v_out/v_target
%   obj_2D.v : Ns x 3 surface vertices
%   obj_2D.f.v : Ms x 3 surface faces
%   v_target : Nx3 initial deployed positions
%   ratio_max: max allowed anisotropy ratio (lambda_max/lambda_min)
%   alphaMax : maximum internal angle (radians)
%   opts     : struct with fields nIter, step, w_len, w_ang, w_smooth, verbose
%              optional: qFrac
%
% Output:
%   v_target : updated deployed positions on the target surface
%   info     : diagnostics

if nargin < 7
    opts = struct();
end
opts = fill_default_opts(opts);
if nargin < 5 || isempty(ratio_max)
    ratio_max = 1.7 / 1.1;
end

if size(v_out, 2) == 2
    v_out = [v_out, zeros(size(v_out,1),1, 'like', v_out)];
elseif size(v_out,2) ~= 3
    error('v_out must be Nx2 or Nx3.');
end
if size(v_target,2) ~= 3
    error('v_target must be Nx3.');
end
if size(f_out,2) ~= 3
    error('f_out must be Mx3.');
end

% Surface data
surfV = obj_2D.v;
surfF = obj_2D.f.v;
if size(surfV,2) ~= 3 || size(surfF,2) ~= 3
    error('obj_2D.v must be Ns x 3 and obj_2D.f.v must be Ms x 3.');
end

% Precompute mesh data
[E, vertNbrs] = build_edge_list_and_vertex_neighbors(f_out, size(v_target,1));
L0 = edge_lengths(v_out, E);
L0_safe = max(L0, eps(class(L0)));

triData = precompute_surface_data(surfV, surfF);

% Initial projection to surface + anchor faces
[v_target, anchorFace, faceNormalsAtPts] = project_points_to_surface(v_target, surfV, surfF, triData, zeros(size(v_target,1),1));

nIter = opts.nIter;
histMaxScale = zeros(nIter,1);
histMaxAngle = zeros(nIter,1);
histRmsUpdate = zeros(nIter,1);

for it = 1:nIter
    % --------- scale constraints on edges ---------
    [lam_e, L1] = compute_scale_factors(v_target, E, L0_safe);
    ratioMax = ratio_max;
    qFrac = opts.qFrac;
    lam_min_cur = min(lam_e);
    lam_max_cur = max(lam_e);
    ratio = lam_max_cur / max(lam_min_cur, eps(class(lam_e)));
    maxScaleViol = max(0, ratio - ratioMax);

    g = zeros(size(v_target), 'like', v_target);

    % Apply anisotropy-ratio correction only when violated.
    if maxScaleViol > 0
        nE = numel(lam_e);
        nSel = max(1, ceil(qFrac * nE));
        [~, idSort] = sort(lam_e, 'ascend');
        idLo = idSort(1:nSel);
        idHi = idSort(max(1, nE-nSel+1):nE);

        lam_med = median(lam_e);
        lam_med = max(lam_med, eps(class(lam_e)));

        targetLam = lam_med;

        % High-stretch edges: push shorter.
        for kk = 1:numel(idHi)
            eIdx = idHi(kk);
            i = E(eIdx,1);
            j = E(eIdx,2);

            xi = v_target(i,:);
            xj = v_target(j,:);
            dij = xi - xj;
            L = max(L1(eIdx), eps(class(L1)));
            dir = dij / L;

            w_hi = max(0, lam_e(eIdx) / lam_med - 1);
            targetL = targetLam * L0_safe(eIdx);
            delta = L - targetL;
            corr = -opts.w_len * w_hi * delta * dir;

            g(i,:) = g(i,:) + corr;
            g(j,:) = g(j,:) - corr;
        end

        % Low-stretch edges: push longer.
        for kk = 1:numel(idLo)
            eIdx = idLo(kk);
            i = E(eIdx,1);
            j = E(eIdx,2);

            xi = v_target(i,:);
            xj = v_target(j,:);
            dij = xi - xj;
            L = max(L1(eIdx), eps(class(L1)));
            dir = dij / L;

            w_lo = max(0, 1 - lam_e(eIdx) / lam_med);
            targetL = targetLam * L0_safe(eIdx);
            delta = L - targetL;
            corr = -opts.w_len * w_lo * delta * dir;

            g(i,:) = g(i,:) + corr;
            g(j,:) = g(j,:) - corr;
        end
    end

    % --------- angle constraints on triangles ---------
    [triAngles, maxTriAngle, maxAngleViol] = triangle_internal_angles(v_target, f_out, alphaMax);

    violTri = find(maxTriAngle > alphaMax);
    for tt = 1:numel(violTri)
        t = violTri(tt);
        tri = f_out(t,:);
        a = tri(1); b = tri(2); c = tri(3);

        [~, kLoc] = max(triAngles(t,:));
        if kLoc == 1
            iV = a; jV = b; kV = c;
        elseif kLoc == 2
            iV = b; jV = c; kV = a;
        else
            iV = c; jV = a; kV = b;
        end

        excess = maxTriAngle(t) - alphaMax;
        midJK = 0.5 * (v_target(jV,:) + v_target(kV,:));
        dirI = midJK - v_target(iV,:);

        nrm = norm(dirI);
        if nrm > 0
            dirI = dirI / nrm;
            w = opts.w_ang * excess;
            g(iV,:) = g(iV,:) + w * dirI;
            g(jV,:) = g(jV,:) - 0.5 * w * dirI;
            g(kV,:) = g(kV,:) - 0.5 * w * dirI;
        end
    end

    % --------- Laplacian smoothing (optional) ---------
    if opts.w_smooth > 0
        lap = zeros(size(v_target), 'like', v_target);
        for vi = 1:size(v_target,1)
            nb = vertNbrs{vi};
            if ~isempty(nb)
                lap(vi,:) = mean(v_target(nb,:), 1) - v_target(vi,:);
            end
        end
        g = g + opts.w_smooth * lap;
    end

    % --------- project update to surface tangent plane ---------
    n = faceNormalsAtPts;
    dotgn = sum(g .* n, 2);
    g_t = g - dotgn .* n;

    dv = opts.step * g_t;
    v_try = v_target + dv;

    % --------- reproject back to surface ---------
    [v_new, anchorFace, faceNormalsAtPts] = project_points_to_surface(v_try, surfV, surfF, triData, anchorFace);

    rmsUpdate = sqrt(mean(sum((v_new - v_target).^2, 2)));
    v_target = v_new;

    histMaxScale(it) = maxScaleViol;
    histMaxAngle(it) = maxAngleViol;
    histRmsUpdate(it) = rmsUpdate;

    if opts.verbose
        fprintf('[reparameterization] iter %d/%d | maxScaleViol=%.4e | maxAngleViol=%.4e | rmsUpdate=%.4e\n', ...
            it, nIter, maxScaleViol, maxAngleViol, rmsUpdate);
    end

    if maxScaleViol < 1e-6 && maxAngleViol < 1e-6 && rmsUpdate < 1e-8
        histMaxScale = histMaxScale(1:it);
        histMaxAngle = histMaxAngle(1:it);
        histRmsUpdate = histRmsUpdate(1:it);
        break;
    end
end

[finalScaleFactors, ~] = compute_scale_factors(v_target, E, L0_safe);
[~, finalMaxAnglePerTri, ~] = triangle_internal_angles(v_target, f_out, alphaMax);

info = struct();
info.history.maxScaleViolation = histMaxScale;
info.history.maxAngleViolation = histMaxAngle;
info.history.rmsUpdate = histRmsUpdate;
info.finalScaleFactors = finalScaleFactors;
info.finalMaxAngle = max(finalMaxAnglePerTri);
info.edgeList = E;

end

function opts = fill_default_opts(opts)
if ~isfield(opts, 'nIter') || isempty(opts.nIter), opts.nIter = 60; end
if ~isfield(opts, 'step') || isempty(opts.step), opts.step = 0.2; end
if ~isfield(opts, 'w_len') || isempty(opts.w_len), opts.w_len = 1.0; end
if ~isfield(opts, 'w_ang') || isempty(opts.w_ang), opts.w_ang = 0.5; end
if ~isfield(opts, 'w_smooth') || isempty(opts.w_smooth), opts.w_smooth = 0.05; end
if ~isfield(opts, 'verbose') || isempty(opts.verbose), opts.verbose = false; end
if ~isfield(opts, 'qFrac') || isempty(opts.qFrac), opts.qFrac = 0.10; end
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

function L = edge_lengths(V, E)
d = V(E(:,1),:) - V(E(:,2),:);
L = sqrt(sum(d.^2, 2));
end

function [lam_e, L1] = compute_scale_factors(v_target, E, L0)
L1 = edge_lengths(v_target, E);
lam_e = L1 ./ max(L0, eps(class(L0)));
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

% vertex -> faces adjacency
nV = size(surfV,1);
v2f = cell(nV,1);
for f = 1:F
    vv = surfF(f,:);
    v2f{vv(1)}(end+1) = f;
    v2f{vv(2)}(end+1) = f;
    v2f{vv(3)}(end+1) = f;
end

% face -> face adjacency via shared edges
E = [surfF(:,[1 2]); surfF(:,[2 3]); surfF(:,[3 1])];
fid = [(1:F)'; (1:F)'; (1:F)'];
E = sort(E,2);
[Eu, ~, ic] = unique(E, 'rows'); %#ok<ASGLU>
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

        % Expand one more ring if candidate list is tiny.
        if numel(cand) < 6
            ring2 = cand;
            for ii = 1:numel(cand)
                ring2 = [ring2; faceNbrs{cand(ii)}(:)]; %#ok<AGROW>
            end
            cand = unique(ring2);
        end
    else
        % global coarse search by face-center distance (one-time expensive fallback)
        d2 = sum((fc - q).^2, 2);
        [~, idx] = sort(d2, 'ascend');
        K = min(30, numel(idx));
        cand = idx(1:K);
    end

    [cp, fBest, ~] = closest_point_on_mesh_faces(q, surfV, surfF, cand);

    if fBest == 0
        % Absolute fallback over all faces
        [cp, fBest, ~] = closest_point_on_mesh_faces(q, surfV, surfF, (1:size(surfF,1)).');
    end

    v_proj(i,:) = cp;
    anchorFace(i) = fBest;
    n_at_point(i,:) = faceNormals(fBest,:);
end

% Safety fallback if some normals are zero
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
% Closest point on triangle (Ericson, Real-Time Collision Detection).
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

% Inside face region
n = cross(ab, ac);
n2 = dot(n, n);
if n2 <= eps
    cp = a;
    return;
end
cp = p - (dot(p - a, n) / n2) * n;
end

function D2 = pdist2_local(A, B)
% Squared pairwise distances (toolbox-free)
AA = sum(A.^2, 2);
BB = sum(B.^2, 2)';
D2 = max(AA + BB - 2*(A*B'), 0);
end
