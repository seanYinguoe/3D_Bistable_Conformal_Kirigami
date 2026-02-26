function scale_facs = calculate_scale_facs(v_out, v_target, f_out)
%CALCULATE_SCALE_FACS Triangle edge-wise scale factors from 2D->3D meshes.
%   scale_facs = calculate_scale_facs(v_out, v_target, f_out)
%
% Inputs:
%   v_out    : Nx2 or Nx3 flattened vertices (if Nx2, z=0 assumed)
%   v_target : Nx3 deployed vertices, same vertex ordering as v_out
%   f_out    : Mx3 triangle connectivity (1-based)
%
% Output:
%   scale_facs : Mx3 edge stretch ratios per triangle:
%                col1 -> edge (n1,n2)
%                col2 -> edge (n2,n3)
%                col3 -> edge (n3,n1)

if size(v_out,2) == 2
    v_out = [v_out, zeros(size(v_out,1),1, 'like', v_out)];
elseif size(v_out,2) ~= 3
    error('v_out must be Nx2 or Nx3.');
end

if size(v_target,2) ~= 3
    error('v_target must be Nx3.');
end

if size(v_out,1) ~= size(v_target,1)
    error('v_out and v_target must have the same number of vertices.');
end

if size(f_out,2) ~= 3
    error('f_out must be Mx3 triangle connectivity.');
end

% Triangle vertex indices
n1 = f_out(:,1);
n2 = f_out(:,2);
n3 = f_out(:,3);

% Flattened-edge lengths L0
L0_12 = vecnorm(v_out(n1,:) - v_out(n2,:), 2, 2);
L0_23 = vecnorm(v_out(n2,:) - v_out(n3,:), 2, 2);
L0_31 = vecnorm(v_out(n3,:) - v_out(n1,:), 2, 2);

% Deployed-edge lengths L1
L1_12 = vecnorm(v_target(n1,:) - v_target(n2,:), 2, 2);
L1_23 = vecnorm(v_target(n2,:) - v_target(n3,:), 2, 2);
L1_31 = vecnorm(v_target(n3,:) - v_target(n1,:), 2, 2);

% Robust degenerate-edge handling
tol = 1e-12;
scale_facs = nan(size(f_out,1), 3, 'like', v_target);

ok12 = L0_12 > tol;
ok23 = L0_23 > tol;
ok31 = L0_31 > tol;

scale_facs(ok12,1) = L1_12(ok12) ./ L0_12(ok12);
scale_facs(ok23,2) = L1_23(ok23) ./ L0_23(ok23);
scale_facs(ok31,3) = L1_31(ok31) ./ L0_31(ok31);

if any(~ok12 | ~ok23 | ~ok31)
    warning('calculate_scale_facs:DegenerateEdge', ...
        'Degenerate reference edge(s) detected (L0 ~ 0). Corresponding scale factors set to NaN.');
end

%{
% Usage test (manual verification for one triangle)
% v_out    = [0 0; 1 0; 0 1];
% v_target = [0 0 0; 2 0 0; 0 3 0];
% f_out    = [1 2 3];
% sf = calculate_scale_facs(v_out, v_target, f_out)
% Expected: sf(1,1)=2, sf(1,2)=sqrt(13)/sqrt(2), sf(1,3)=3
%}

end
