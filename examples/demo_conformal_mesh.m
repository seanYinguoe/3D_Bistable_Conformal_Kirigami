function result = demo_conformal_mesh(model, makePlot)
%DEMO_CONFORMAL_MESH Inspect the supplied BFF map and local area expansion.
% The UV map is precomputed; this does not run BFF or prove bistability.
if nargin < 1, model = 'quarter_dome'; end
if nargin < 2, makePlot = true; end
model = validatestring(model,{'quarter_dome','double_dome','hemisphere'});
root = setup_project;
obj = readObj(fullfile(root,'data','meshes',model),[model '_flat.obj']);
uv = vertice_sort(obj.vt,obj.f.v,obj.f.vt);
Aflat = triangle_area_2D(uv,obj.f.v);
Acurved = triangle_area_3D(obj.v,obj.f.v);
assert(all(Aflat>0) && all(Acurved>0),'The mesh contains degenerate faces.');
lambda_area = sqrt(Acurved./Aflat);
assert(all(isfinite(lambda_area)),'Nonfinite area scale.');
out = tempname(fullfile(root,'results')); mkdir(out);
result = struct('model',model,'uv',uv,'vertices',obj.v,'faces',obj.f.v, ...
    'lambda_area',lambda_area,'output_dir',out,'matlab_version',version);
save(fullfile(out,'conformal_mesh.mat'),'result');
if makePlot
    fig = figure('Color','w','Position',[100 100 1000 430]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile; patch('Faces',obj.f.v,'Vertices',uv,'FaceVertexCData',lambda_area, ...
        'FaceColor','flat','EdgeColor','none');
    axis equal off; title('Flattened coordinates: local area scale'); colorbar;
    nexttile; patch('Faces',obj.f.v,'Vertices',obj.v,'FaceVertexCData',lambda_area, ...
        'FaceColor','flat','EdgeColor','none');
    view(35,25); axis equal off; title('Target surface'); colorbar;
    colormap(parula);
    exportgraphics(fig,fullfile(out,'conformal_mesh.png'),'Resolution',180);
end
fprintf('%s: %d vertices, %d faces. Area scale range %.3f to %.3f.\n', ...
    model,size(obj.v,1),size(obj.f.v,1),min(lambda_area),max(lambda_area));
fprintf('Saved to %s\n',out);
end
