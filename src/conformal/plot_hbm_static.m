L = 1.0; A0=[0,0]; B0=[L,0];
A1=A0; B1=B0;                     % only impose rotations
thetaA = 10*pi/180;               % 10 deg at A
thetaB = -10*pi/180;              % -10 deg at B
N = 50;

E  = 7.3e9;    % Young's modulus (杨氏模量)
b  = 1.0;    % width (宽度)
t  = 0.01;   % thickness (厚度)
A  = b*t;
I  = b*t^3/12;

a0 = L/N;                        % segment rest length
Ks = E*A / a0;                   % per-bar axial stiffness
Kb = E*I / a0;                   % per-hinge bending stiffness

[Etotal, XY, out] = hbm_static(A0,B0,A1,B1,thetaA,thetaB,N,Kb,Ks, ...
                               'w_pos',1e6, 'w_ang',1e8, 'use_fminunc', true);

figure; plot(XY(:,1),XY(:,2),'-o'); axis equal; grid on;
title(sprintf('Hencky bar-chain shape, E=%.2g J',Etotal));
fprintf('N = %d and Total elastic energy = %.4f J\n', N, Etotal);