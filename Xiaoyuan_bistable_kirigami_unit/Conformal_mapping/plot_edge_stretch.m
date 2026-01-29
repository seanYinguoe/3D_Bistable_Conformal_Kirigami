function plot_edge_stretch(lambda, v_out, E)
% Plot 2D mesh edges and label stretch factor on each edge
%
% lambda : ne×1   stretch factor for each edge
% v_out  : n×2    2D vertex coordinates
% E      : ne×2   edge connectivity (i,j)

figure;
hold on; 
axis equal;
axis off;

% draw edges
for k = 1:size(E,1)
    p1 = v_out(E(k,1),:);
    p2 = v_out(E(k,2),:);
    plot([p1(1) p2(1)], [p1(2) p2(2)], 'k-');
end

title('Stretch Factor lambda = L_{3D} / L_{2D}');
xlabel('x'); ylabel('y');

% label strech factor on edge midpoints
mid = 0.5 * (v_out(E(:,1),:) + v_out(E(:,2),:));

for k = 1:size(E,1)
    text(mid(k,1), mid(k,2), sprintf('%.2f', lambda(k)), ...
        'FontSize', 8, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'Color','r');
end
end