function neighbors_idx = pointCloudClosestPoints(uv_c, query_point, search_radius)
% Finds the indices of mesh centroids (uv_c) that are within a given
% search radius from a grid centroid (query_point).
% Inputs:
% - uv_c: (N x 2) or (N x 3) matrix of mesh centroids (N points in 2D/3D)
% - query_point: (1 x 2) or (1 x 3) row vector representing the grid centroid
% - search_radius: Scalar specifying the maximum allowed distance
% 
% Output:
% - neighbors_idx: Indices of uv_c points that are within the search radius
%

    % Compute Euclidean distances from query_point to all uv_c points
    distances = vecnorm(uv_c - query_point, 2, 2); 
    
    % Find indices where the distance is within the search radius
    neighbors_idx = find(distances <= search_radius);
    
end
