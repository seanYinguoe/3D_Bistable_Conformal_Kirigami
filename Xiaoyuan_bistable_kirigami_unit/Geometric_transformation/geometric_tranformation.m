%% Input 3D Geometric information
% Read the file path
modelname = 'hemisphere';
path = strcat('/Users/sean/Desktop/Project 2/3D_Bistable_Conformal_Kirigami/Xiaoyuan_bistable_kirigami_unit/Input_model/Reference_model/',modelname,'/');
filename_3D = strcat(modelname,'_flat.obj'); % 2D model
obj_3D = readObj(path,filename_3D); % read 2D object, vertices, connectivity

% Change the order the nodes in vt the connectivity doesn't match
sorted_uv = vertice_sort(obj_3D.vt,obj_3D.f.v,obj_3D.f.vt);
obj_3D.vt = sorted_uv;

%% Rotate the geometry on xy plane
flattened_surface = obj_3D.vt;
flattened_surface = [flattened_surface,zeros(size(flattened_surface,1),1)];
deployed_surface = obj_3D.v;    % Deployed 3D vertex positions
[obj_3D.vt, obj_3D.v] = model_rotate(modelname, flattened_surface, deployed_surface);

% Read the nodes and connectivity information
v_mesh = obj_3D.vt;
f_mesh = obj_3D.f.v;

%% Remesh the 3D geometry(triangulation)
% Isotropic remeshing  


%% Plot the geometry
figure()
patch('Vertices', obj_3D.v, 'Faces', obj_3D.f.v, ...
    'FaceColor', 'none', 'EdgeColor', 'black', 'FaceAlpha', 0.6);
view(3)
grid off;
axis off
axis equal;