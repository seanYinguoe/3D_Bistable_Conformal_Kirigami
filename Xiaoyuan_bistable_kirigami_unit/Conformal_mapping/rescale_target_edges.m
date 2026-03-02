function [v_target, scale_facs, info] = rescale_target_edges(v_out, f_out, obj_2D, v_target, s, opts)
%RESCALE_TARGET_EDGES Uniformly rescale deployed geometry and remove drift.
%   [v_target, scale_facs, info] = rescale_target_edges(v_out, f_out, obj_2D, v_target, s, opts)

if nargin < 6 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'do_project') || isempty(opts.do_project)
    opts.do_project = false;
end

if size(v_target,2) ~= 3
    error('v_target must be Nx3.');
end
if size(f_out,2) ~= 3
    error('f_out must be Mx3.');
end
if ~isscalar(s) || ~isfinite(s) || s <= 0
    error('s must be a finite positive scalar.');
end

v_target0 = v_target;
c0 = mean(v_target0, 1);
v_scaled = (v_target0 - c0) * s + c0;

[R, t, rms_fit] = kabsch_rigid(v_scaled, v_target0);
v_target = v_scaled * R + t;

if opts.do_project
    surfV = obj_2D.v;
    surfF = obj_2D.f.v;
    triData = precompute_surface_data(surfV, surfF);
    anchorFace0 = zeros(size(v_target,1), 1);
    [v_target, anchorFace0] = project_points_to_surface(v_target, surfV, surfF, triData, anchorFace0); %#ok<NASGU>
end

scale_facs = calculate_scale_facs(v_out, v_target, f_out);

E = build_edge_list(f_out);
Lold = edge_lengths(v_target0, E);
Lnew = edge_lengths(v_target, E);
ratio = Lnew ./ max(Lold, eps(class(Lold)));

info = struct();
info.scale_target = s;
info.R = R;
info.t = t;
info.rms_fit = rms_fit;
info.edge_scale_check = [min(ratio) max(ratio)];
end

function [R, t, rms_fit] = kabsch_rigid(A, B)
if size(A,2) ~= 3 || size(B,2) ~= 3 || size(A,1) ~= size(B,1)
    error('A and B must be Nx3 with matching size.');
end

cA = mean(A, 1);
cB = mean(B, 1);
A0 = A - cA;
B0 = B - cB;

H = A0' * B0;
[U, ~, V] = svd(H);
R = U * V';
if det(R) < 0
    U(:,end) = -U(:,end);
    R = U * V';
end

t = cB - cA * R;
A_fit = A * R + t;
diffAB = A_fit - B;
rms_fit = sqrt(mean(sum(diffAB.^2, 2)));
end

function E = build_edge_list(f_out)
E = [f_out(:,[1 2]); f_out(:,[2 3]); f_out(:,[3 1])];
E = sort(E, 2);
E = unique(E, 'rows');
end

function L = edge_lengths(V, E)
d = V(E(:,1),:) - V(E(:,2),:);
L = sqrt(sum(d.^2, 2));
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

function [v_proj, anchorFace] = project_points_to_surface(vq, surfV, surfF, triData, anchorFace)
nq = size(vq,1);
v_proj = zeros(nq,3, 'like', vq);

if nargin < 5 || isempty(anchorFace)
    anchorFace = zeros(nq,1);
end

fc = triData.faceCenters;
faceNbrs = triData.faceNbrs;

for i = 1:nq
    q = vq(i,:);
    if anchorFace(i) >= 1 && anchorFace(i) <= size(surfF,1)
        cand = [anchorFace(i); faceNbrs{anchorFace(i)}(:)];
        cand = unique(cand);
    else
        d2 = sum((fc - q).^2, 2);
        [~, idx] = sort(d2, 'ascend');
        cand = idx(1:min(30, numel(idx)));
    end

    [cp, fBest] = closest_point_on_mesh_faces(q, surfV, surfF, cand);
    if fBest == 0
        [cp, fBest] = closest_point_on_mesh_faces(q, surfV, surfF, (1:size(surfF,1)).');
    end

    v_proj(i,:) = cp;
    anchorFace(i) = fBest;
end
end

function [cpBest, fBest] = closest_point_on_mesh_faces(q, surfV, surfF, faceIdx)
cpBest = [0 0 0];
fBest = 0;
d2Best = inf;
for kk = 1:numel(faceIdx)
    f = faceIdx(kk);
    tri = surfF(f,:);
    cp = closest_point_on_triangle(q, surfV(tri(1),:), surfV(tri(2),:), surfV(tri(3),:));
    d2 = sum((q - cp).^2);
    if d2 < d2Best
        d2Best = d2;
        cpBest = cp;
        fBest = f;
    end
end
end

function cp = closest_point_on_triangle(p, a, b, c)
ab = b - a; ac = c - a; ap = p - a;
d1 = dot(ab, ap); d2 = dot(ac, ap);
if d1 <= 0 && d2 <= 0, cp = a; return; end

bp = p - b;
d3 = dot(ab, bp); d4 = dot(ac, bp);
if d3 >= 0 && d4 <= d3, cp = b; return; end

vc = d1*d4 - d3*d2;
if vc <= 0 && d1 >= 0 && d3 <= 0
    v = d1 / (d1 - d3);
    cp = a + v * ab;
    return;
end

cpv = p - c;
d5 = dot(ab, cpv); d6 = dot(ac, cpv);
if d6 >= 0 && d5 <= d6, cp = c; return; end

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
if n2 <= eps, cp = a; return; end
cp = p - (dot(p - a, n) / n2) * n;
end
