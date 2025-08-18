function angle = angle_calculate(nodes)
% node1, node2, node3 are 2x1 matrices representing the (x,y) coordinates of each node

% Calculate the vectors formed by the nodes
vector1 = nodes(2,:) - nodes(1,:);
vector2 = nodes(3,:) - nodes(2,:);
vector3 = nodes(1,:) - nodes(3,:);

% Calculate the magnitudes of the vectors
mag1 = norm(vector1);
mag2 = norm(vector2);
mag3 = norm(vector3);

% Calculate the angles
angle = acos((mag1^2 + mag2^2 - mag3^2) / (2 * mag1 * mag2));
end
