function info = plot_edge_scale_factors(v_target, f_out, scale_facs)
%PLOT_EDGE_SCALE_FACTORS Edge-wise scale-factor visualization on deployed 3D mesh.
%   info = plot_edge_scale_factors(v_target, f_out, scale_facs)
%
% Inputs
%   v_target   : #V x 3 deployed vertices
%   f_out      : #F x 3 triangle connectivity
%   scale_facs : edge-wise scale values, accepted formats:
%                (A) #F x 3, with per-face edge order:
%                    col1 -> edge (n1,n2)
%                    col2 -> edge (n2,n3)
%                    col3 -> edge (n3,n1)
%                (B) (3*#F) x 1 using stacked order [e12; e23; e31]
%                (C) #E_unique x 1 already mapped to unique edges
%
% Notes on edge correspondence:
%   Unique edges are built from f_out by sorting node indices in each edge
%   and calling unique(...,'rows','stable'). If scale_facs is face-wise, shared
%   edges are averaged (ignoring NaNs) to get one value per unique edge.

validateattributes(v_target, {'numeric'}, {'2d','ncols',3,'finite','real'});
validateattributes(f_out, {'numeric'}, {'2d','ncols',3,'integer','positive'});

nF = size(f_out,1);

% Build triangle-edge list in the declared per-face order.
e12 = f_out(:, [1 2]);
e23 = f_out(:, [2 3]);
e31 = f_out(:, [3 1]);
edges_face = [e12; e23; e31];              % (3*nF) x 2
edges_key  = sort(edges_face, 2);          % undirected key
[edges_unique, ~, ic] = unique(edges_key, 'rows', 'stable');
nE = size(edges_unique,1);

% Map scale_facs onto unique edges.
edge_vals = map_scale_to_unique_edges(scale_facs, nF, nE, ic);

% -------- Plot style (clean, publication-friendly) --------
figure('Color','w');
hold on; box on;

% Very light surface context layer.
patch('Vertices', v_target, 'Faces', f_out, ...
      'FaceColor', [0.92 0.95 0.98], 'FaceAlpha', 0.10, ...
      'EdgeColor', 'none');

% Colormap and value range.
cm = summer(256);
colormap(cm);

finite_mask = isfinite(edge_vals);
if any(finite_mask)
    vmin = min(edge_vals(finite_mask));
    vmax = max(edge_vals(finite_mask));
else
    vmin = 0; vmax = 1;
end
if abs(vmax - vmin) < eps
    vmax = vmin + 1;
end
caxis([vmin vmax]);

% Draw each unique edge with edge-wise color.
P1 = v_target(edges_unique(:,1), :);
P2 = v_target(edges_unique(:,2), :);
for e = 1:nE
    val = edge_vals(e);
    if ~isfinite(val)
        c = [0.75 0.75 0.75]; % fallback color for NaN edges
    else
        t = (val - vmin) / (vmax - vmin);
        t = min(max(t, 0), 1);
        idx = 1 + round(t * (size(cm,1)-1));
        c = cm(idx, :);
    end
    plot3([P1(e,1) P2(e,1)], [P1(e,2) P2(e,2)], [P1(e,3) P2(e,3)], ...
        '-', 'Color', c, 'LineWidth', 1.1);
end

% Visual cleanup.
axis equal;
axis off;
view(3);
camlight('headlight');
lighting gouraud;

c = colorbar;
c.FontSize = 16;
c.Label.String = 'Edge Scale Factor';
c.Label.FontSize = 16;

info = struct();
info.edges_unique = edges_unique;
info.edge_values = edge_vals;
info.assumed_face_edge_order = '[n1,n2], [n2,n3], [n3,n1]';
info.n_unique_edges = nE;
info.value_range = [vmin, vmax];
end

function edge_vals = map_scale_to_unique_edges(scale_facs, nF, nE, ic)
% Convert various scale_facs layouts to one value per unique edge.

if ismatrix(scale_facs) && size(scale_facs,1) == nF && size(scale_facs,2) == 3
    vals_face = [scale_facs(:,1); scale_facs(:,2); scale_facs(:,3)];
elseif isvector(scale_facs) && numel(scale_facs) == 3*nF
    vals_face = scale_facs(:);
elseif isvector(scale_facs) && numel(scale_facs) == nE
    edge_vals = scale_facs(:);
    return;
else
    error(['scale_facs size is ambiguous. Expected #F x 3, (3*#F)x1, ' ...
           'or #E_unique x 1.']);
end

vals_face = vals_face(:);
ok = isfinite(vals_face);
sum_vals = accumarray(ic(ok), vals_face(ok), [nE 1], @sum, NaN);
cnt_vals = accumarray(ic(ok), 1, [nE 1], @sum, 0);

edge_vals = sum_vals ./ max(cnt_vals, 1);
edge_vals(cnt_vals == 0) = NaN;
end

