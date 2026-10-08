function plot_surface_result(result)
%PLOT_SURFACE_RESULT Compact overview; computational routines create no UI.
fig=figure('Color','w','Position',[100 100 1200 400]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile;patch('Faces',result.faces,'Vertices',result.target_vertices, ...
    'FaceVertexCData',min(result.edge_stretch,[],2),'FaceColor','flat', ...
    'EdgeColor',[.3 .35 .4]);view(35,25);axis equal off;colorbar;
title('Discrete target: minimum edge stretch');
nexttile;hold on;
colours={'white',[.25 .35 .43],[.08 .18 .24],[.25 .35 .43]};
for k=1:numel(result.compact),plot_triangle(result.compact{k},colours);end
axis equal off;title('Flat cut pattern');
nexttile;
patch('Faces',result.faces,'Vertices',result.flat_vertices, ...
    'FaceVertexCData',result.beta*180/pi,'FaceColor','flat','EdgeColor','none');
axis equal off;colorbar;title('Assigned tilt angle (degrees)');
exportgraphics(fig,fullfile(result.output_dir,'surface_overview.png'),'Resolution',160);
end
