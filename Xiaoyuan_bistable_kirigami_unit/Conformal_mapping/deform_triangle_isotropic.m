function [triangle_new,E_total] = deform_triangle_isotropic(delta,edgeLen,l1,l4,beta,t)
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

%% Define local function
rotation = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)]; % rotation matrix  
rotrow = @(v, ang) (rotation(ang) * v(:))';   % row -> column -> rotate -> row
nrm1   = @(v) v / norm(v);          % normalize

%% Calculate geometric parameters
A = pi/3-beta;
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(A); % length of filament
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) ...
     - cos(pi/3 - beta) .* l4;
l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 ...
     - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) );
R = sqrt(3)/3 * l3;
N=20; E=4.3e11; b=1.0;
%% Generate uniformaly deployed triangle as initial guess
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;
[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);
triangle_new = triangle;

%% Define optimise parameters
% Nodes in flanks
d1 = triangle(20,:);
d1_ = triangle(19,:);

% Extract original inner triangle vertices as initial guess
B_orig = triangle(44,:);

edge = edgeLen + delta; % deformed length of a unit
centroid = [-sqrt(3)/6*edge,-1/2*edge];
theta0 = atan2(B_orig(2)-centroid(2), B_orig(1)-centroid(1));

%% Difine the objective function HBM and add constraints as penalty
    function [E_oneLig, pack] = energy(th)
        % place B on the circle around centroid
        B = centroid + R*[cos(th), sin(th)];
        % get A by rotating B around centroid by -120°
        A = centroid + (B - centroid) * rotation(2*pi/3);

        % build boundary vectors as
        vecA0 = (d1_ - d1);   % normal at the fixed flank end
        vecB0 = (A  - B );    % edge-normal at inner triangle vertex
        deltaTilt = -pi/3;    % rotate to initial state
        vecA  = nrm1( rotrow(vecA0,  deltaTilt) );
        vecB  = nrm1( rotrow(vecB0,  deltaTilt) );
        % single ligament energy from hbm_energy (L0 = l2)
        [E_lig, ~, ~, ~] = hbm_energy(l2, d1, B, vecA, vecB, N, E, b, t, ...
                'VectorsAreNormals', true);

        % non-overlap penalty
        res = ifoverlapping([d1_; d1; B]);     % res < 0 → OK, >0 → overlap
        P_ol = 1e6 * max(0,res)^2;

        E_oneLig = E_lig + P_ol;     % one-ligament energy with penalty

        if nargout>1
            pack.B = B; pack.A = A;
        end
    end

%% minimize over theta

% search a window around theta0; widen if needed
[theta_opt, E_total, ~] = fminbnd(@(th) energy(th), theta0 - pi/4, theta0 + pi/4);

%% Create optimised unit
% rebuild final geometry at theta_opt
[~, pack] = energy(theta_opt);
B_new = pack.B;
A_new = pack.A;
C_new = centroid + (B_new - centroid) * rotation( -2*pi/3); 

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
end
