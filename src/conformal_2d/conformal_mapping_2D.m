function [uv, f_uv, info] = conformal_mapping_2D(v, f)
% CONFORMAL_MAPPING_2D  Harmonic parameterization of an open triangle mesh
%                       to a rectangle with arc-length proportional boundary.
%
% [uv, f_uv, info] = conformal_mapping_2D(v, f)
%
% INPUTS
%   v  : (#V x 2) or (#V x 3) vertex positions
%   f  : (#F x 3) triangle indices, 1-based
%
% OUTPUTS
%   uv     : (#V x 2) 2-D parameterization in [0,W] x [0,H]
%   f_uv   : (#F x 3) same connectivity as f
%   info   : struct with fields
%              .bnd_loop   - ordered boundary vertex indices  (#B x 1)
%              .rect       - [W H] dimensions of output rectangle
%              .corners    - 4 boundary vertex indices nearest rectangle corners
%              .angle_def  - mean absolute angle defect of interior vertices
%              .flipped    - number of inverted triangles in UV
%
% ALGORITHM
%   1. Detect & order the single boundary loop.
%   2. Compute cumulative arc-length along boundary.
%   3. Map boundary vertices onto a rectangle by arc-length proportion.
%      Aspect ratio is chosen from the ratio of opposite side arc-lengths.
%   4. Build cotangent-weight Laplacian L (sparse, #V x #V).
%   5. Solve  L_ii * uv_i = -L_ib * uv_b  for x and y independently.
%   6. Normalize result to [0,W] x [0,H] and compute diagnostics.
%
% REQUIREMENTS  Pure MATLAB, no toolboxes beyond the sparse solver.

% =========================================================================
% 0. Input validation & dimension handling
% =========================================================================
assert(ismatrix(v) && size(v,2) >= 2, ...
    'v must be (#V x 2) or (#V x 3).');
assert(ismatrix(f) && size(f,2) == 3, ...
    'f must be (#F x 3) with 1-based indices.');
assert(all(f(:) >= 1) && all(f(:) <= size(v,1)), ...
    'f contains out-of-range vertex indices.');

nV = size(v, 1);
nF = size(f,  1);

% Work in 3-D internally (pad with zeros if 2-D input)
if size(v,2) == 2
    v3 = [v, zeros(nV,1)];
else
    v3 = v(:,1:3);
end

% =========================================================================
% 1. Detect boundary loop
% =========================================================================
bnd_loop = detect_boundary_loop(f, nV);   % may error internally

% =========================================================================
% 2. Arc-length parameterisation of boundary
% =========================================================================
nB = numel(bnd_loop);
bnd_pts = v3(bnd_loop, :);                % (#B x 3) positions in order

% Edge lengths around the loop (closed: last vertex back to first)
diffs   = [bnd_pts(2:end,:); bnd_pts(1,:)] - bnd_pts;
seg_len = sqrt(sum(diffs.^2, 2));         % (#B x 1)
cum_len = [0; cumsum(seg_len)];           % (#B+1 x 1), cum_len(end)=total
total_L = cum_len(end);

if total_L < eps
    error('conformal_mapping_2D: boundary has zero arc-length.');
end

t = cum_len(1:end-1) / total_L;          % normalised parameter in [0,1)

% -------------------------------------------------------------------------
% Choose rectangle aspect ratio
%   The boundary is traversed counter-clockwise.  We split it into 4 sides
%   at the parametric midpoints 0, 0.25, 0.5, 0.75 and use the arc-length
%   ratio of horizontal vs vertical sides to set W and H.
%   Simple heuristic: side0+side2 (width-like), side1+side3 (height-like).
% -------------------------------------------------------------------------
breaks   = [0, 0.25, 0.5, 0.75, 1.0];   % parametric corners
side_len = zeros(1,4);
for s = 1:4
    mask = (t >= breaks(s)) & (t < breaks(s+1));
    side_len(s) = sum(seg_len(mask));
end
% Avoid degenerate aspect ratios
horiz = max(side_len(1) + side_len(3), eps);
vert  = max(side_len(2) + side_len(4), eps);
ratio = horiz / vert;                     % W/H

% Normalise so that area ≈ 1 (cosmetic; will be renormalised later)
H = 1.0;
W = ratio;

% -------------------------------------------------------------------------
% Map each boundary vertex to rectangle edge by arc-length
%   Corners (t): 0 → (0,0), 0.25 → (W,0), 0.5 → (W,H), 0.75 → (0,H)
% -------------------------------------------------------------------------
uv_bnd = zeros(nB, 2);
corner_t = [0, 0.25, 0.5, 0.75];
corner_xy = [0,0; W,0; W,H; 0,H];

for k = 1:nB
    tk = t(k);
    if tk < 0.25          % bottom edge: (0,0) → (W,0)
        s = tk / 0.25;
        uv_bnd(k,:) = (1-s)*corner_xy(1,:) + s*corner_xy(2,:);
    elseif tk < 0.5       % right edge:  (W,0) → (W,H)
        s = (tk - 0.25) / 0.25;
        uv_bnd(k,:) = (1-s)*corner_xy(2,:) + s*corner_xy(3,:);
    elseif tk < 0.75      % top edge:    (W,H) → (0,H)
        s = (tk - 0.5) / 0.25;
        uv_bnd(k,:) = (1-s)*corner_xy(3,:) + s*corner_xy(4,:);
    else                  % left edge:   (0,H) → (0,0)
        s = (tk - 0.75) / 0.25;
        uv_bnd(k,:) = (1-s)*corner_xy(4,:) + s*corner_xy(1,:);
    end
end

% =========================================================================
% 3. Cotangent Laplacian
% =========================================================================
L = build_cotan_laplacian(v3, f, nV);

% =========================================================================
% 4. Partition: boundary (b) and interior (i)
% =========================================================================
is_bnd = false(nV, 1);
is_bnd(bnd_loop) = true;
idx_all  = (1:nV)';
interior = idx_all(~is_bnd);
nI = numel(interior);

if nI == 0
    % Trivial: every vertex is on boundary, no interior to solve
    uv = zeros(nV, 2);
    uv(bnd_loop, :) = uv_bnd;
    f_uv = f;
    uv = normalise_uv(uv, W, H);
    info = build_info(bnd_loop, uv, f, W, H, t);
    return
end

% Index map: global → interior row
g2i = zeros(nV, 1);
g2i(interior) = (1:nI)';

% =========================================================================
% 5. Solve harmonic system:  L_ii * uv_i = -L_ib * uv_b
% =========================================================================
L_ii = L(interior, interior);          % (nI x nI) sparse
L_ib = L(interior, bnd_loop);          % (nI x nB) sparse

rhs_x = -L_ib * uv_bnd(:,1);
rhs_y = -L_ib * uv_bnd(:,2);

% Use backslash (MATLAB's sparse LU / Cholesky fallback)
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
sol_x = L_ii \ rhs_x;
sol_y = L_ii \ rhs_y;
warning('on','MATLAB:singularMatrix');
warning('on','MATLAB:nearlySingularMatrix');

if ~all(isfinite(sol_x)) || ~all(isfinite(sol_y))
    error('conformal_mapping_2D: linear solve produced non-finite values. Check mesh connectivity.');
end

% =========================================================================
% 6. Assemble full UV
% =========================================================================
uv = zeros(nV, 2);
uv(bnd_loop,  :) = uv_bnd;
uv(interior, 1)  = sol_x;
uv(interior, 2)  = sol_y;

% =========================================================================
% 7. Post-process: normalise to [0,W] x [0,H]
% =========================================================================
uv = normalise_uv(uv, W, H);

% =========================================================================
% 8. Outputs
% =========================================================================
f_uv = f;
info = build_info(bnd_loop, uv, f, W, H, t);

end % ── main function ──────────────────────────────────────────────────────


% =========================================================================
% LOCAL HELPER: detect_boundary_loop
% =========================================================================
function bnd_loop = detect_boundary_loop(f, nV)
% Returns an ordered (#B x 1) list of boundary vertex indices.
% Errors if mesh has 0 or >1 boundary loops, or is non-manifold.

    nF = size(f, 1);

    % Build all half-edges and count occurrences of undirected edges
    % A boundary edge appears exactly once.
    all_edges = [f(:,1), f(:,2);
                 f(:,2), f(:,3);
                 f(:,3), f(:,1)];          % (3*nF x 2) directed

    % Sort each edge so a < b (undirected key)
    edge_key = sort(all_edges, 2);         % (3*nF x 2)

    % Find edges that appear exactly once → boundary
    [uniq_edges, ~, ic] = unique(edge_key, 'rows');
    counts = accumarray(ic, 1);
    bnd_mask = (counts == 1);

    if ~any(bnd_mask)
        error('conformal_mapping_2D: no boundary found. Mesh may be closed (genus 0) or degenerate.');
    end

    % Recover directed boundary half-edges
    % Among the directed half-edges, pick those whose undirected key is boundary
    is_bnd_he = bnd_mask(ic);             % (3*nF x 1) logical
    bnd_he    = all_edges(is_bnd_he, :);  % (#bnd_edges x 2) directed

    % Build adjacency list: next[a] = b
    nBE = size(bnd_he, 1);
    adj = containers.Map('KeyType','int32','ValueType','int32');
    for k = 1:nBE
        a = int32(bnd_he(k,1));
        b = int32(bnd_he(k,2));
        if isKey(adj, a)
            error('conformal_mapping_2D: non-manifold boundary detected at vertex %d.', a);
        end
        adj(a) = b;
    end

    % Walk all loops
    visited = false(nV, 1);
    all_loops = {};
    starts = int32(bnd_he(:,1));

    for s = 1:numel(starts)
        v0 = starts(s);
        if visited(v0), continue; end
        loop = [];
        cur = v0;
        for steps = 1:nBE+1
            if visited(cur) && cur ~= v0
                break;   % hit a previously completed loop
            end
            if ~isKey(adj, cur)
                error('conformal_mapping_2D: boundary is not a closed loop (vertex %d has no successor).', cur);
            end
            loop(end+1) = cur; %#ok<AGROW>
            visited(cur) = true;
            cur = adj(cur);
            if cur == v0, break; end
        end
        all_loops{end+1} = loop(:); %#ok<AGROW>
    end

    if numel(all_loops) == 0
        error('conformal_mapping_2D: failed to detect any boundary loop.');
    elseif numel(all_loops) > 1
        sizes = cellfun(@numel, all_loops);
        error('conformal_mapping_2D: mesh has %d boundary loops (sizes: %s). Requires exactly 1.', ...
              numel(all_loops), num2str(sizes));
    end

    bnd_loop = double(all_loops{1});
end


% =========================================================================
% LOCAL HELPER: build_cotan_laplacian
% =========================================================================
function L = build_cotan_laplacian(v, f, nV)
% Assembles the cotangent-weight Laplacian (sparse, symmetric, nV x nV).
% L(i,j) = -0.5*(cot(alpha_ij) + cot(beta_ij))  for adjacent i,j
% L(i,i) = -sum_{j~i} L(i,j)
%
% Robustness: clamp cotangents to avoid blow-up on degenerate triangles.

    nF = size(f, 1);
    II = zeros(6*nF, 1);
    JJ = zeros(6*nF, 1);
    VV = zeros(6*nF, 1);
    ptr = 0;

    for fi = 1:nF
        idx = f(fi, :);                 % [i j k]
        P   = v(idx, :);               % 3x3, rows are vertices

        % Edge vectors
        e = [P(3,:)-P(2,:);            % edge opposite vertex 1
             P(1,:)-P(3,:);            % edge opposite vertex 2
             P(2,:)-P(1,:)];           % edge opposite vertex 3  (= -e(1)-e(2))

        % Cotangent of each angle via dot/cross
        % angle at vertex k is between edges from k to the other two
        cots = zeros(1,3);
        for k = 1:3
            a = k;               % vertex index in local {1,2,3}
            b = mod(k,3)+1;
            c = mod(k+1,3)+1;
            ea = P(b,:) - P(a,:);
            ec = P(c,:) - P(a,:);
            cos_a = dot(ea, ec);
            sin_a = norm(cross(ea, ec));
            cots(k) = cos_a / (sin_a + eps);
        end
        % Clamp to avoid degenerate triangles
        cots = max(min(cots, 1e6), -1e6);

        % Cotangent weight for each edge:
        %   edge (i,j) opposite vertex k  → weight = 0.5 * cot(angle_k)
        % Edges: (2,3) opposite k=1, (1,3) opposite k=2, (1,2) opposite k=3
        pairs = [2,3; 1,3; 1,2];
        for e_loc = 1:3
            a = idx(pairs(e_loc,1));
            b = idx(pairs(e_loc,2));
            w = 0.5 * cots(e_loc);
            % Off-diagonal entries: L(a,b) += -w, L(b,a) += -w
            % Diagonal entries:     L(a,a) += +w, L(b,b) += +w
            ptr = ptr + 1;
            II(ptr) = a; JJ(ptr) = b; VV(ptr) = -w;
            ptr = ptr + 1;
            II(ptr) = b; JJ(ptr) = a; VV(ptr) = -w;
            ptr = ptr + 1;
            II(ptr) = a; JJ(ptr) = a; VV(ptr) = +w;
            ptr = ptr + 1;
            II(ptr) = b; JJ(ptr) = b; VV(ptr) = +w;
        end
    end

    II = II(1:ptr);
    JJ = JJ(1:ptr);
    VV = VV(1:ptr);
    L  = sparse(II, JJ, VV, nV, nV);
end


% =========================================================================
% LOCAL HELPER: normalise_uv
% =========================================================================
function uv = normalise_uv(uv, W, H)
% Shift and scale UV so it fits exactly in [0,W] x [0,H].

    mn = min(uv, [], 1);
    mx = max(uv, [], 1);
    rng = mx - mn;
    rng(rng < eps) = 1;      % avoid division by zero

    uv(:,1) = (uv(:,1) - mn(1)) / rng(1) * W;
    uv(:,2) = (uv(:,2) - mn(2)) / rng(2) * H;
end


% =========================================================================
% LOCAL HELPER: build_info
% =========================================================================
function info = build_info(bnd_loop, uv, f, W, H, t)
% Assemble diagnostic struct.

    info.bnd_loop = bnd_loop;
    info.rect     = [W, H];

    % Corner indices: boundary vertices closest to t = 0, 0.25, 0.5, 0.75
    corner_t = [0, 0.25, 0.5, 0.75];
    corners  = zeros(1,4);
    for c = 1:4
        [~, ci]    = min(abs(t - corner_t(c)));
        corners(c) = bnd_loop(ci);
    end
    info.corners = corners;

    % Angle defect at interior vertices (Gaussian curvature proxy)
    %   For a flat map, angle defect should be ~0 at interior vertices.
    nV = size(uv,1);
    is_bnd = false(nV,1);
    is_bnd(bnd_loop) = true;
    interior = find(~is_bnd);

    % Sum of triangle angles at each interior vertex in UV
    angle_sum = zeros(nV,1);
    for fi = 1:size(f,1)
        idx = f(fi,:);
        P   = uv(idx,:);
        for k = 1:3
            a = k; b = mod(k,3)+1; c = mod(k+1,3)+1;
            ea = P(b,:) - P(a,:);
            ec = P(c,:) - P(a,:);
            na = norm(ea); nc = norm(ec);
            if na > eps && nc > eps
                cos_a = dot(ea,ec)/(na*nc);
                cos_a = max(min(cos_a,1),-1);
                angle_sum(idx(a)) = angle_sum(idx(a)) + acos(cos_a);
            end
        end
    end
    defects = abs(2*pi - angle_sum(interior));
    info.angle_def = mean(defects);

    % Count flipped triangles (signed area < 0 in UV)
    n_flip = 0;
    for fi = 1:size(f,1)
        idx = f(fi,:);
        A = uv(idx(1),:); B = uv(idx(2),:); C = uv(idx(3),:);
        area2 = (B(1)-A(1))*(C(2)-A(2)) - (B(2)-A(2))*(C(1)-A(1));
        if area2 < 0
            n_flip = n_flip + 1;
        end
    end
    info.flipped = n_flip;
end


% =========================================================================
%
%  USAGE DEMO (run this block as a script after saving this file):
%
% -------------------------------------------------------------------------
%
%   % --- Create a simple triangulated disk mesh ---
%   n_rings   = 10;
%   n_sectors = 32;
%   [vd, fd]  = make_disk_mesh(n_rings, n_sectors);   % see sub-function below
%
%   % --- Run conformal mapping ---
%   [uv, f_uv, info] = conformal_mapping_2D(vd, fd);
%
%   fprintf('Rectangle: %.4f x %.4f\n', info.rect(1), info.rect(2));
%   fprintf('Mean angle defect: %.4e rad\n', info.angle_def);
%   fprintf('Flipped triangles: %d\n', info.flipped);
%
%   % --- Plot ---
%   figure('Name','Conformal Mapping Demo','Color','w');
%
%   subplot(1,2,1); hold on; axis equal; title('Original Mesh (XY)');
%   triplot(fd, vd(:,1), vd(:,2), 'b-', 'LineWidth', 0.3);
%   bnd = info.bnd_loop;
%   plot(vd([bnd;bnd(1)],1), vd([bnd;bnd(1)],2), 'r-', 'LineWidth', 2);
%   xlabel('X'); ylabel('Y');
%
%   subplot(1,2,2); hold on; axis equal;
%   title(sprintf('UV Parameterization [%.2f x %.2f]', info.rect(1), info.rect(2)));
%   triplot(f_uv, uv(:,1), uv(:,2), 'g-', 'LineWidth', 0.3);
%   bnd_uv = info.bnd_loop;
%   plot(uv([bnd_uv;bnd_uv(1)],1), uv([bnd_uv;bnd_uv(1)],2), 'r-', 'LineWidth', 2);
%   rectangle('Position',[0 0 info.rect(1) info.rect(2)],'EdgeColor','k','LineWidth',2);
%   xlabel('U'); ylabel('V');
%
% -------------------------------------------------------------------------
% Helper: make_disk_mesh
%   function [v, f] = make_disk_mesh(n_rings, n_sectors)
%       v = [0 0 0];
%       f = [];
%       for r = 1:n_rings
%           radius = r / n_rings;
%           for s = 1:n_sectors
%               angle = 2*pi*(s-1)/n_sectors;
%               v(end+1,:) = [radius*cos(angle), radius*sin(angle), 0];
%           end
%       end
%       % Fan from center to first ring
%       for s = 1:n_sectors
%           a = 1 + s;
%           b = 1 + mod(s, n_sectors) + 1;
%           f(end+1,:) = [1, a, b];
%       end
%       % Ring-to-ring quads split into triangles
%       for r = 1:n_rings-1
%           base_i = 1 + (r-1)*n_sectors;
%           base_o = 1 +  r   *n_sectors;
%           for s = 1:n_sectors
%               i0 = base_i + s - 1;
%               i1 = base_i + mod(s, n_sectors);
%               o0 = base_o + s - 1;
%               o1 = base_o + mod(s, n_sectors);
%               f(end+1,:) = [i0+1, o0+1, i1+1];
%               f(end+1,:) = [i1+1, o0+1, o1+1];
%           end
%       end
%   end
%
% =========================================================================