function plot_edge_stretch(scale_data, v_out, conn)
% Plot 2D mesh edges and label stretch factor on each unique edge.
%
% Usage 1 (new):
%   plot_edge_stretch(scale_facs, v_out, f_out)
%   scale_facs : nt×3 triangle-wise edge stretch factors
%   f_out      : nt×3 triangle connectivity
%
% Usage 2 (legacy):
%   plot_edge_stretch(lambda, v_out, E)
%   lambda     : ne×1 edge stretch factors
%   E          : ne×2 edge connectivity

if size(v_out,2) > 2
    v2 = v_out(:,1:2);
else
    v2 = v_out;
end

if size(conn,2) == 2
    % Legacy mode: already edge-wise
    E = sort(conn, 2);
    lambda_e = scale_data(:);
else
    % New mode: triangle-wise -> edge-wise (average on shared edges)
    f_out = conn;
    if size(f_out,2) ~= 3 || size(scale_data,2) ~= 3 || size(scale_data,1) ~= size(f_out,1)
        error('For triangle input, require scale_facs(nt,3) and f_out(nt,3).');
    end

    e12 = sort(f_out(:,[1 2]), 2);
    e23 = sort(f_out(:,[2 3]), 2);
    e31 = sort(f_out(:,[3 1]), 2);

    E_all = [e12; e23; e31];
    lam_all = [scale_data(:,1); scale_data(:,2); scale_data(:,3)];

    [E, ~, ic] = unique(E_all, 'rows');
    sum_lam = accumarray(ic, lam_all, [], @sum);
    cnt_lam = accumarray(ic, 1, [], @sum);
    lambda_e = sum_lam ./ max(cnt_lam, 1);
end

figure;
hold on;
axis equal;
axis off;

for k = 1:size(E,1)
    p1 = v2(E(k,1),:);
    p2 = v2(E(k,2),:);
    plot([p1(1) p2(1)], [p1(2) p2(2)], 'k-');
end

title('Stretch Factor \lambda = L_{3D} / L_{2D}');
xlabel('x');
ylabel('y');

mid = 0.5 * (v2(E(:,1),:) + v2(E(:,2),:));
for k = 1:size(E,1)
    text(mid(k,1), mid(k,2), sprintf('%.2f', lambda_e(k)), ...
        'FontSize', 8, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Color','r');
end
end
