% Deploy a flattened surface onto the deployed surface (conformal mapping)
% using interpolation and analysis.
% Input: obj_2D(flattened_surface and deployed_surface)
function [flattened_surface,deployed_surface] = mesh_deployment(obj_2D)
% Input flattened and deployed surface
flattened_surface = obj_2D.vt;  % Flattened 2D vertex positions
flattened_surface = [flattened_surface,zeros(size(flattened_surface,1),1)];
deployed_surface = obj_2D.v;    % Deployed 3D vertex positions

% Keep flattened platform fixed, and rigidly align deployed mesh to it.
% This avoids model-dependent global rotation drift during interpolation.
[R, t] = kabsch_rigid(deployed_surface, flattened_surface);
deployed_surface = deployed_surface * R + t;

% Keep deployed surface above the flattened plane for clearer deployment.
if mean(deployed_surface(:,3)) < 0
    deployed_surface(:,3) = -deployed_surface(:,3);
end

% Recenter platform: move flattened centroid to (0,0) in XY,
% and apply the same XY shift to deployed surface.
c_flat_xy = mean(flattened_surface(:,1:2), 1);
flattened_surface(:,1:2) = flattened_surface(:,1:2) - c_flat_xy;
deployed_surface(:,1:2) = deployed_surface(:,1:2) - c_flat_xy;

end

function [R, t] = kabsch_rigid(A, B)
% Solve rigid map A*R + t ~= B (no scaling, no reflection).
if size(A,2) ~= 3 || size(B,2) ~= 3 || size(A,1) ~= size(B,1)
    error('kabsch_rigid expects A,B as Nx3 with matching N.');
end

cA = mean(A, 1);
cB = mean(B, 1);
A0 = A - cA;
B0 = B - cB;

H = A0' * B0;
[U, ~, V] = svd(H);
R = U * V';
if det(R) < 0
    U(:,end) = -U(:,end);
    R = U * V';
end
t = cB - cA * R;
end
