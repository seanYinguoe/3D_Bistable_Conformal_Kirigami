function [E,lambda] = edge_stretch(V3,L,F)
% Input:
% V3: nodes coordinates in 3D(deployed)
% L:  edge length in 2D grid
% F: connectivity of faces
% Output:
% lambda = 3D edge length / 2D edge length(stretch factor)

% Sort out the indices for edges
E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
E = unique(sort(E,2),'rows');        % undirected + unique

% Calculate the length of edge in 3D
L3 = vecnorm(V3(E(:,1),:) - V3(E(:,2),:), 2, 2);

% Calculate the stretch factor
lambda = L3 ./ L;
end