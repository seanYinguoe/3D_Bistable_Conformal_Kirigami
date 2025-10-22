%% bistable unit
% % Helpers
% nrm1  = @(v) v./norm(v); % normlize vectors
% rotation = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)]; % rotation matrix
% angle = @(v1,v2) acos(dot(v1,v2)/(norm(v1)*norm(v2)));
% rotrow = @(v, ang) (rotation(ang) * v(:))';   % rotate a vector
%
% % In bisatble unit
% triangle_test = triangle;
% A1 = triangle_test(20,:); B1 = triangle_test(44,:);
% vecA0 = nrm1(triangle_test(19,:) - triangle_test(20,:));
% vecB0 = nrm1(triangle_test(43,:) - triangle_test(44,:));
% norm_beam = nrm1(B1-A1)*rotation(-pi/2);
% tilt_A1 = angle(vecA0,norm_beam);
% tilt_B1 = angle(vecB0,norm_beam);
% vecA1  = nrm1( rotrow(vecA0,  tilt_A1) );
% vecB1  = nrm1( rotrow(vecB0,  tilt_B1) );
% L0 = norm(B1-A1);
% N=20; E=4.3e11; b=1.0; t=0.015*edgeLen;I = b*t^3/12;EI = E*I;
%% simply support beam
% Example data
L0 = 1;
A0=[0,0]; B0=[1,0];
A1=[0,0]; B1=[1,0];
vecA0 = [0,1];
vecB0 = [0,1];
vecA1=[cos(pi/2-pi/10),sin(pi/2-pi/10)]; vecB1=[cos(pi/2+pi/10),sin(pi/2+pi/10)];      % clamped end directions (normals)
%vecA1=[cos(pi/2),sin(pi/2)]; vecB1=[cos(pi/2),sin(pi/2)];      % clamped end directions (normals)
N=30;
E=4.3e11; b=1.0; t=0.01;I = b*t^3/12;

% Run solver (returns deformed config only)
[Etotal, XYdef, springs] = hbm_energy(L0, A1,B1, vecA1,vecB1, N, E,b, t, ...
    'VectorsAreNormals', true);

fprintf('N = %d and Total elastic energy = %g J\n', N, Etotal);

%% build undeformed geometry
L0 = norm(B0 - A0);
a0 = L0/N;
phi0 = atan2((B0(2)-A0(2)), (B0(1)-A0(1)));   % angle of undeformed bar
XYund = zeros(N+1,2);
XYund(1,:) = A0;
for i=1:N
    XYund(i+1,:) = XYund(i,:) + a0*[cos(phi0), sin(phi0)];
end

%% plot
figure(); hold on; axis equal;
title('Hencky bar-chain: deformed');

% undeformed beam
plot(XYund(:,1),XYund(:,2),'-o','Color',[0 0.4 1],'DisplayName','Undeformed');

% deformed beam
plot(XYdef(:,1),XYdef(:,2),'-o','Color',[1 0 0],'DisplayName','Deformed');
%plot_triangle(triangle);

legend show
xlabel('x'); ylabel('y');