function [R, T] = align_triangles(p1, p2, p3, q1, q2, q3)
% Inputs: Corresponding 1x3 row vectors (p1↔q1, p2↔q2, p3↔q3)
% Outputs: Rotation matrix R (3x3), translation vector t (3x1)

% Source and target points (3x3 matrices)
source = [p1; p2; p3];
target = [q1; q2; q3];

% Centroids
centroid_src = mean(source, 1);
centroid_tgt = mean(target, 1);

% Center points
source_centered = source - centroid_src;
target_centered = target - centroid_tgt;

% Covariance matrix
H = source_centered' * target_centered;

% SVD decomposition
[U, ~, V] = svd(H);

% Rotation matrix
R = V * U';
if det(R) < 0
    V(:,3) = -V(:,3);
    R = V * U';
end

% Translation vector
T = centroid_tgt' - R * centroid_src';
end
