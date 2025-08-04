function res = ifoverlapping(nodes)
% nodes include node1, node2,node3 are 1x2 matrices representing the (x,y) coordinates of each node

% Calculate the vectors formed by the nodes
vector1 = nodes(1,:) - nodes(2,:);
vector2 = nodes(3,:) - nodes(2,:);


% Calculate the cross product of the vectors
cross_nodes = cross([vector1,0],[vector2,0]);

% Determine whether overlap occurs
res = dot(cross_nodes,[0 0 1]);
end
