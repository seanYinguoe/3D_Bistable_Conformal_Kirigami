function [triangle,E_total] = deform_triangle_semi(delta,edgeLen,l1,l4,beta,t)
% DEFORM_TRIANGLE Deforms a isotropic triangle based on input node positions and edge lengths.
%
% Inputs:
%   q1, q2, q3   - Unit node coordinates (1x3 vectors)
%   edgeLen      - Original triangle edge length
%   l1, l4, beta   - l1: length of flanks l4:thickness of flanks
%   beta:tilting angle
%   t            - Thickness of filaments
%   i_out        - Orientation flag (0 for upwards, 1 for downwards)
%
% Output:
%   triangle_new - Deformed triangle coordinates

%% Define helpers
rotation = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)]; % rotation matrix  
rotrow = @(v, ang) (rotation(ang) * v(:))';   % row -> column -> rotate -> row
nrm1   = @(v) v./norm(v);          % normalize
angle = @(v1,v2) acos(dot(v1,v2)/(norm(v1)*norm(v2))); 

%% Calculate geometric parameters
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(pi/3-beta); % length of filament
N=30; E=4.3e11; b=1.0;
%% Generate uniformaly deployed triangle as initial guess
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;
[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);

%% Calculate energy based on given geometry
% Nodes in flanks
d1 = triangle(20,:);
d1_ = triangle(19,:);
% Extract triangle
B = triangle(44,:);
A = triangle(43,:);
% Ratate beam end direction as energy free state
vecA0 = (d1_ - d1);   % normal at the fixed flank end
vecB0 = (A  - B );    % edge-normal at inner triangle vertex
% norm_beam = nrm1(B-d1)*rotation(-pi/2);
% tilt_A = angle(vecA0,norm_beam);
% tilt_B = angle(vecB0,norm_beam);
tilt_A = pi/6;
tilt_B = pi/6;
vecA  = nrm1( rotrow(vecA0,  tilt_A) );
vecB  = nrm1( rotrow(vecB0,  tilt_B) );
% Calculate single ligament energy from hbm_energy (L0 = l2)
L0 = norm(B-d1);
[E_lig,XYdef,springs] = hbm_energy(L0, d1, B, vecA, vecB, N, E, b, t, ...
    'VectorsAreNormals', true);
E_total = E_lig;

%% Plot the results
% figure(); axis equal;
% title('Hencky bar-chain: deformed');
% 
% % deformed beam
% plot(XYdef(:,1),XYdef(:,2),'-o','Color',[1 0 0],'DisplayName','Deformed');
% plot_triangle(triangle);
% 
% legend show
% xlabel('x'); ylabel('y');
end
