% Example data
L0 = l2;
A1 = triangle(20,:); B1 = triangle(44,:);
vecA1 = triangle(19,:) - triangle(20,:);
vecB1 = triangle(43,:) - triangle(44,:);
% L0 = 1;
% A1=[0,0]; B1=[1,0];
%vecA1=[cos(pi/2-pi/10),sin(pi/2-pi/10)]; vecB1=[cos(pi/2+pi/10),sin(pi/2+pi/10)];      % clamped end directions (normals)
%vecA1=[cos(pi/2),sin(pi/2)]; vecB1=[cos(pi/2),sin(pi/2)];      % clamped end directions (normals)
N=20; E=4.3e11; b=1.0; t=0.015*edgeLen;I = b*t^3/12;EI = E*I;

% Run solver (returns deformed config only)
[Etotal, XYdef, springs] = hbm_energy(L0, A1,B1, vecA1,vecB1, N, E,b, t, ...
                                      'VectorsAreNormals', true);

fprintf('N = %d and Total elastic energy = %.4f J\n', N, Etotal);

% ---- build undeformed geometry ----
% L0 = norm(B0 - A0);
% a0 = L0/N;
% phi0 = atan2((B0(2)-A0(2)), (B0(1)-A0(1)));   % angle of undeformed bar
% XYund = zeros(N+1,2);
% XYund(1,:) = A0;
% for i=1:N
%     XYund(i+1,:) = XYund(i,:) + a0*[cos(phi0), sin(phi0)];
% end

% ---- plot ----
figure; hold on; axis equal;
title('Hencky bar-chain: undeformed vs deformed');

% undeformed beam
% plot(XYund(:,1),XYund(:,2),'-o','Color',[0 0.4 1],'DisplayName','Undeformed');

% deformed beam
plot(XYdef(:,1),XYdef(:,2),'-o','Color',[1 0 0],'DisplayName','Deformed');

legend show
xlabel('x'); ylabel('y');