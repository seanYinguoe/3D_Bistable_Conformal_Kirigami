function vq_target = map_points_barycentric_2D(vq, v_source, f_source, v_target)
%MAP_POINTS_BARYCENTRIC_2D Transfer points through a planar mesh mapping.
%   vq_target = map_points_barycentric_2D(vq, v_source, f_source, v_target)
%
% Inputs
%   vq       : Nq x2 query points in source domain
%   v_source : Ns x2 source mesh vertices
%   f_source : Mx3 source mesh faces
%   v_target : Ns x2 target mesh vertices (same indexing/connectivity)
%
% Output
%   vq_target: Nq x2 mapped query points in target domain

if size(vq,2) ~= 2 || size(v_source,2) ~= 2 || size(v_target,2) ~= 2
    error('vq, v_source, and v_target must be Nx2 arrays.');
end
if size(v_source,1) ~= size(v_target,1)
    error('v_source and v_target must have the same number of vertices.');
end
if size(f_source,2) ~= 3
    error('f_source must be Mx3 triangle connectivity.');
end

TR = triangulation(f_source, v_source);
[tid, bc] = pointLocation(TR, vq);

vq_target = zeros(size(vq), 'like', vq);
for i = 1:size(vq,1)
    if isfinite(tid(i)) && tid(i) >= 1
        tri = f_source(tid(i), :);
        vq_target(i,:) = bc(i,:) * v_target(tri,:);
    else
        % Robust fallback for points just outside triangulation due to numerical
        % tolerances: nearest source vertex correspondence.
        idx = dsearchn(v_source, vq(i,:));
        vq_target(i,:) = v_target(idx,:);
    end
end
end
