function [triangle_new,E_total] = deform_triangle_isotropic(delta,edgeLen,l1,l4,beta,t)
% DEFORM_TRIANGLE Deforms a isotropic triangle based on input node positions and edge lengths.
%
% Inputs:
%   delta     : scalar or vector of signed edge-length changes (+ extend, - shorten)
%   edgeLen   : original equilateral outer edge length
%   l1, l4    : flank geometry parameters (your notation)
%   beta      : flank tilt angle (rad)
%   t         : filament thickness (used for visual offset of outer clamp normal only
% Output:
%   triangle_new - Deformed triangle coordinates

%% Define helpers
rotation = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)]; % rotation matrix  
rotrow = @(v, ang) (rotation(ang) * v(:))';   % row -> column -> rotate -> row
nrm1   = @(v) v / norm(v);          % normalize

%% Calculate geometric parameters
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(pi/3-beta); % length of filament
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) ...
     - cos(pi/3 - beta) .* l4;
l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 ...
     - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) ); % length of inner triangle
R = sqrt(3)/3 * l3;
N=10; E=4.3e11; b=1.0;
%% Generate uniformaly deployed triangle as initial guess
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;
[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);
triangle_new = triangle;

%% Define optimise parameters
% Nodes in flanks
d1 = triangle(20,:);
d1_ = triangle(21,:);

% Extract original inner triangle vertices as initial guess
B_orig = triangle(44,:);

edge = edgeLen + delta; % deformed length of a unit
centroid = [-sqrt(3)/6*edge,-1/2*edge];
theta0 = atan2(B_orig(2)-centroid(2), B_orig(1)-centroid(1));

%% Difine the objective function HBM and add constraints as penalty
    function [E_oneLig, pack,XYdef, springs] = energy(th)
        % place B on the circle around centroid
        B = centroid + R*[cos(th), sin(th)];
        % get C by rotating B around centroid by -120°
        C = centroid + (B - centroid) * rotation(-2*pi/3);
        % build boundary vectors as
        vecA = (d1 - d1_);   % tangent at the fixed flank end
        vecB = (C  - B );    % edge-normal at inner triangle vertex
        % single ligament energy from hbm_energy (L0 = l2)
        [E_oneLig,XYdef, springs] = hbm_energy(l2, d1, B, vecA, vecB, N, E, b, t*sqrt(3)/2, ...
                'VectorsAreNormals',false);  % get projection t
        if nargout>1
            pack.B = B; pack.C = C;
        end
    end

%% minimize over theta
% search a window over theta
options = optimset('TolX',1e-6,'TolFun',1e-10,'MaxFunEvals',10,'Display','off');
[theta_opt, E_total] = fminbnd(@(th) energy(th), theta0 - pi/20, theta0 + pi/20, options);
E_total = E_total * 3;

%% Create optimised unit
% rebuild final geometry at theta_opt
[~, pack,XYdef, ~] = energy(theta_opt);
B_new = pack.B;
C_new = pack.C;
A_new = centroid + (B_new - centroid) * rotation(2*pi/3); 

triangle_new(43,:) = A_new;
triangle_new(44,:) = B_new;
triangle_new(45,:) = C_new;

% Create the new filaments
triangle_new(32,:) = triangle_new(20,:);
triangle_new(33,:) = triangle_new(44,:);
triangle_new(34,:) = triangle_new(33,:) + t/(norm(A_new-B_new))*(A_new-B_new);

triangle_new(37,:) = triangle_new(45,:);
triangle_new(36,:) = triangle_new(24,:);
triangle_new(38,:) = triangle_new(37,:) + t/(norm(B_new-C_new))*(B_new-C_new);

triangle_new(40,:) = triangle_new(28,:);
triangle_new(41,:) = triangle_new(43,:);
triangle_new(42,:) = triangle_new(41,:) + t/(norm(C_new-A_new))*(C_new-A_new);

% Plot results
figure()
plot_triangle(triangle_new);
plot(XYdef(:,1),XYdef(:,2),'-o','Color',[1 0 0],'DisplayName','Deformed');

end

