% Create interpolants for x, y, z coordinates
F_x = scatteredInterpolant(obj_2D.vt(:,1), obj_2D.vt(:,2), obj_2D.v(:,1),'natural');
F_y = scatteredInterpolant(obj_2D.vt(:,1), obj_2D.vt(:,2), obj_2D.v(:,2),'natural');
F_z = scatteredInterpolant(obj_2D.vt(:,1), obj_2D.vt(:,2), obj_2D.v(:,3),'natural');

% Interpolate grid points to the deployed surface
v_grid_deploy = zeros(size(v_grid,1), 3);
v_grid_deploy(:,1) = F_x(v_grid(:,1), v_grid(:,2));
v_grid_deploy(:,2) = F_y(v_grid(:,1), v_grid(:,2));
v_grid_deploy(:,3) = F_z(v_grid(:,1), v_grid(:,2));

figure()
patch('Vertices', v_grid_deploy, 'Faces', f_out,'FaceColor', 'none', 'EdgeColor', 'b','LineWidth',1);
c.FontSize = 18;
axis equal; 
axis off



