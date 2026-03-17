function export_obj(V, F, filename)
%EXPORT_OBJ Export vertices/faces to a Wavefront OBJ file.
%   export_obj(V, F)
%   export_obj(V, F, filename)
%
% Inputs
%   V        : Nx3 (or Nx2) vertices
%   F        : Mx3 faces (1-based indexing)
%   filename : output OBJ file path (optional, default: 'mesh.obj')

if nargin < 3 || isempty(filename)
    filename = 'mesh.obj';
end

if size(V,2) == 2
    V = [V, zeros(size(V,1),1, 'like', V)];
elseif size(V,2) ~= 3
    error('V must be Nx2 or Nx3.');
end

if size(F,2) ~= 3
    error('F must be Mx3.');
end

fid = fopen(filename, 'w');
if fid < 0
    error('Cannot open output file: %s', filename);
end

cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '# Exported by export_obj\n');
fprintf(fid, '# Vertices: %d\n', size(V,1));
fprintf(fid, '# Faces: %d\n', size(F,1));

for i = 1:size(V,1)
    fprintf(fid, 'v %.12g %.12g %.12g\n', V(i,1), V(i,2), V(i,3));
end

for i = 1:size(F,1)
    fprintf(fid, 'f %d %d %d\n', F(i,1), F(i,2), F(i,3));
end

end
