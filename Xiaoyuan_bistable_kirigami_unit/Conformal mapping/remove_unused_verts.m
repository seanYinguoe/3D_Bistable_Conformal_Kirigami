function [v_new, f_new] = remove_unused_verts(v, f)
%{
    Removes unused vertices from a mesh and reindexes the faces.
    
    Inputs:
    - v: Nx3 matrix of vertex coordinates.
    - f: Mx3 matrix of face indices.

    Outputs:
    - v_new: Updated vertex list with only referenced vertices.
    - f_new: Updated face list with new indices.
%}

v_new = [];   % New list of vertices
f_new = zeros(size(f));  % New list of faces (preallocate for speed)
index_map = containers.Map('KeyType', 'double', 'ValueType', 'double'); % Mapping old -> new indices
new_index = 1;

for i = 1:size(f, 1)
    vTemp = round(v(f(i, :), :), 8); % Round for approximate comparison
    face_indices = zeros(1, 3);  % Placeholder for new face indices

    for j = 1:3
        vi = vTemp(j, :);  % Extract the vertex

        % Handle case where v_new is empty
        if isempty(v_new)
            v_new = vi;  % Initialize v_new
            face_indices(j) = new_index;
            index_map(size(v_new, 1)) = new_index;
            new_index = new_index + 1;
        else
            % Check if the vertex already exists in v_new
            [found, idx] = ismember(vi, v_new, 'rows');
            if found
                face_indices(j) = idx;
            else
                v_new = [v_new; vi]; % Append new vertex
                index_map(size(v_new, 1)) = new_index; % Store mapping
                face_indices(j) = new_index;
                new_index = new_index + 1;
            end
        end
    end

    f_new(i, :) = face_indices; % Store new face indices
end
end
