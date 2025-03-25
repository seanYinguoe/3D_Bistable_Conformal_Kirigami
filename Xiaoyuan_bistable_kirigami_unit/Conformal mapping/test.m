% rows = 3;
% cols = 3;
% l = 1;
% ox = 0;
% oy = 0;
%[v_grid, f_grid, c_grid, i_grid, x_grid] = generate_triangular_grid(rows, cols, l, ox, oy);
vt = obj_2D.vt;
[v_grid, f_grid, c_grid, i_grid, x_grid] = generate_overlay_grid(vt, 10);

% Plot the vertices and faces of the grid in a figure window
figure;
hold on;

% Plot vertices as red dots
scatter(v_grid(:, 1), v_grid(:, 2), 'r.');

% Plot each triangular face as blue lines connecting vertices
for i = 1:size(f_grid, 1)
    v_idx = f_grid(i, :); % Vertex indices of the current triangle face
    patch('Vertices', v_grid(:, 1:2), 'Faces', v_idx, 'FaceColor', 'none', 'EdgeColor', 'b');
end
patch('Vertices', obj_2D.vt, 'Faces', obj_2D.f.vt, 'FaceColor', 'none', 'EdgeColor', 'b'); % Plot 2D figure
% Label axes and set title for visualization clarity
xlabel('X');
ylabel('Y');
title('Triangular Grid');
axis equal;
hold off;