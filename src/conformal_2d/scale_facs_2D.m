function scale_facs = scale_facs_2D(v_initial, f_initial, v_target, f_target)
%SCALE_FACS_2D Edge-wise scale factors for planar-to-planar mapping.
%   scale_facs = scale_facs_2D(v_initial, f_initial, v_target, f_target)
%
% Inputs
%   v_initial : Nv x2 or Nv x3 initial vertices
%   f_initial : Nf x3 initial connectivity
%   v_target  : Nv x2 or Nv x3 target vertices
%   f_target  : Nf x3 target connectivity
%
% Output
%   scale_facs : Nf x3 edge stretch ratios per triangle

if size(f_initial,2) ~= 3 || size(f_target,2) ~= 3
    error('f_initial and f_target must be Mx3 triangle connectivity.');
end
if size(f_initial,1) ~= size(f_target,1) || any(f_initial(:) ~= f_target(:))
    error('f_initial and f_target must be identical for edge-wise comparison.');
end
if size(v_initial,1) ~= size(v_target,1)
    error('v_initial and v_target must have the same number of vertices.');
end

if size(v_initial,2) == 2
    v_initial = [v_initial, zeros(size(v_initial,1),1)];
end
if size(v_target,2) == 2
    v_target = [v_target, zeros(size(v_target,1),1)];
end

scale_facs = calculate_scale_facs(v_initial, v_target, f_initial);
end
