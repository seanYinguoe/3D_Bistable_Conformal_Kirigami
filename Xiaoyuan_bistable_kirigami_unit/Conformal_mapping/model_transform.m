function [flattened_out, deployed_out, tf] = model_transform(flattened_surface, deployed_surface)
%MODEL_TRANSFORM Automatic rigid alignment between flattened and deployed surfaces.
%   [flattened_out, deployed_out, tf] = model_transform(flattened_surface, deployed_surface)
%   - Removes global translation and rotation automatically.
%   - Rotates both surfaces so their best-fit planes align with XY plane.
%   - Enforces the flattened surface on z=0 and deployed surface above XY plane.
%   - Aligns flattened in-plane orientation to deployed via 2D rigid fit.
%
% Inputs:
%   flattened_surface : Nx2 or Nx3 (obj_2D.vt)
%   deployed_surface  : Nx3        (obj_2D.v)
%
% Outputs:
%   flattened_out : Nx3 aligned flattened surface (z=0)
%   deployed_out  : Nx3 aligned deployed surface
%   tf            : struct with transforms for applying same mapping to other
%                   flattened-domain points (see tf_apply_flat_points)

if size(flattened_surface,2) == 2
    flattened_surface = [flattened_surface, zeros(size(flattened_surface,1),1)];
elseif size(flattened_surface,2) ~= 3
    error('flattened_surface must be Nx2 or Nx3');
end
if size(deployed_surface,2) ~= 3
    error('deployed_surface must be Nx3');
end
if size(flattened_surface,1) ~= size(deployed_surface,1)
    error('flattened_surface and deployed_surface must have same number of vertices.');
end

% 1) Remove centroids
c_flat = mean(flattened_surface, 1);
c_dep = mean(deployed_surface, 1);
flat0 = flattened_surface - c_flat;
dep0  = deployed_surface  - c_dep;

% 2) Rotate each best-fit plane to XY plane
R_flat_plane = plane_to_xy_rotation(flat0);
R_dep_plane  = plane_to_xy_rotation(dep0);

flat1 = rotate_points(flat0, R_flat_plane);
dep1  = rotate_points(dep0,  R_dep_plane);

% 2.5) Enforce a consistent positive-Z convention for the deployed shape.
% Flip if needed, then shift so the surface sits on/above the XY plane.
if median(dep1(:,3)) < 0
    dep1(:,3) = -dep1(:,3);
    dep_z_flip = -1;
else
    dep_z_flip = 1;
end

dep_z_shift = -min(dep1(:,3));
if dep_z_shift ~= 0
    dep1(:,3) = dep1(:,3) + dep_z_shift;
end

% 3) In-plane rigid alignment (flattened XY -> deployed XY)
Fxy = flat1(:,1:2);
Dxy = dep1(:,1:2);

% Remove tiny mean drift before in-plane solve
Fxy = Fxy - mean(Fxy,1);
Dxy = Dxy - mean(Dxy,1);

Q = inplane_kabsch(Fxy, Dxy); % 2x2
flat_xy_aligned = Fxy * Q;

t_xy = mean(dep1(:,1:2),1) - mean(flat_xy_aligned,1);

% 4) Build outputs
flattened_out = [flat_xy_aligned + t_xy, zeros(size(flat_xy_aligned,1),1)];
deployed_out  = dep1;

% 5) Save transform for other flattened-domain points
% For row point p:
%   p0   = p - c_flat
%   p1   = p0 * R_flat_plane'
%   pxy2 = p1(1:2) * Q + t_xy
%   pout = [pxy2, 0]
tf = struct();
tf.c_flat = c_flat;
tf.c_dep = c_dep;
tf.R_flat_plane = R_flat_plane;
tf.R_dep_plane = R_dep_plane;
tf.Q2 = Q;
tf.t_xy = t_xy;
tf.dep_z_flip = dep_z_flip;
tf.dep_z_shift = dep_z_shift;
end

function R = plane_to_xy_rotation(P)
% Fit normal via SVD and rotate it to +Z.
[~,~,V] = svd(P, 'econ');
n = V(:,3);
if dot(n, [0;0;1]) < 0
    n = -n;
end
R = align_vec_to_vec(n, [0;0;1]);
end

function R = align_vec_to_vec(a, b)
% Rodrigues rotation that maps vector a -> b.
a = a(:); b = b(:);
a = a / max(norm(a), eps);
b = b / max(norm(b), eps);

v = cross(a,b);
s = norm(v);
c = dot(a,b);

if s < 1e-12
    if c > 0
        R = eye(3);
    else
        % 180 deg: pick any axis orthogonal to a
        if abs(a(1)) < 0.9
            k = [1;0;0];
        else
            k = [0;1;0];
        end
        v = cross(a,k);
        v = v / max(norm(v), eps);
        K = [   0   -v(3)  v(2);
              v(3)   0    -v(1);
             -v(2)  v(1)   0   ];
        R = eye(3) + 2*(K*K); % pi rotation
    end
    return;
end

v = v / s;
K = [   0   -v(3)  v(2);
      v(3)   0    -v(1);
     -v(2)  v(1)   0   ];
ang = atan2(s,c);
R = eye(3) + sin(ang)*K + (1-cos(ang))*(K*K);
end

function Q = inplane_kabsch(X, Y)
% Orthogonal 2D matrix Q minimizing ||X*Q - Y||_F.
H = X' * Y;
[U,~,V] = svd(H);
Q = U * V';
if det(Q) < 0
    U(:,end) = -U(:,end);
    Q = U * V';
end
end

function Pout = rotate_points(P, R)
% Row-wise 3D rotation: Pout = (R * P')'
Pout = (R * P')';
end
