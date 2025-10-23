function [triangle_new,E_total,theta_opt] = deform_triangle_batch(delta,edgeLen,l1,l4,beta,t, varargin)
% ---- parse optional inputs ----
p = inputParser;
addParameter(p,'theta_guess',[],@(x) isempty(x) || isscalar(x));
addParameter(p,'do_plot',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
theta_guess = p.Results.theta_guess;

% helpers
rotation = @(theta) [cos(theta),-sin(theta);sin(theta),cos(theta)];
rotrow   = @(v, ang) (rotation(ang) * v(:))';
nrm1     = @(v) v / norm(v);

% geometry
l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) .* sin(pi/3-beta);
l6 = (2/sqrt(3)) .* sin(pi/3 - beta) .* (l1 - 0.5*l4) - cos(pi/3 - beta) .* l4;
l5 = ( (sqrt(3)/2) .* l4 + sin(beta) .* l6 ) ./ sin(pi/3 - beta);
l3 = l6 - l5 - 1.5*l2 - l2 .* ( (sqrt(3)/2) .* (cos(pi/3 - beta) ./ sin(pi/3 - beta)) );
R  = sqrt(3)/3 * l3;
N=10; E=4.3e11; b=1.0;

prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;
[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);
triangle_new = triangle;

% flank nodes
d1  = triangle(20,:);
d1_ = triangle(21,:);

% inner vertex seed
B_orig  = triangle(44,:);
edge    = edgeLen + delta;
centroid= [-sqrt(3)/6*edge,-1/2*edge];
theta0  = atan2(B_orig(2)-centroid(2), B_orig(1)-centroid(1));

% if provided, overwrite initial guess
if ~isempty(theta_guess), theta0 = theta_guess; end

    function [E_oneLig, pack,XYdef, springs] = energy(th)
        B = centroid + R*[cos(th), sin(th)];
        C = centroid + (B - centroid) * rotation(-2*pi/3);
        vecA = (d1 - d1_);
        vecB = (C  - B );
        [E_oneLig,XYdef, springs] = hbm_energy(l2, d1, B, vecA, vecB, N, E, b, t*sqrt(3)/2, ...
                'VectorsAreNormals',false);
        if nargout>1, pack.B = B; pack.C = C; end
    end

% narrow bracket around theta0 for fminbnd
% Small window accelerates as delta changes slowly.
win = pi/25;                            
options = optimset('TolX',1e-6,'TolFun',1e-10,'MaxFunEvals',10,'Display','off');
[theta_opt, E_one] = fminbnd(@(th) energy(th), theta0 - win, theta0 + win, options);
E_total = 3*E_one;  % three ligaments

% rebuild final geometry
[~, pack,XYdef, ~] = energy(theta_opt);
B_new = pack.B; C_new = pack.C;
A_new = centroid + (B_new - centroid) * rotation(2*pi/3); 
triangle_new(43,:) = A_new;
triangle_new(44,:) = B_new;
triangle_new(45,:) = C_new;

triangle_new(32,:) = triangle_new(20,:);
triangle_new(33,:) = triangle_new(44,:);
triangle_new(34,:) = triangle_new(33,:) + t/(norm(A_new-B_new))*(A_new-B_new);
triangle_new(37,:) = triangle_new(45,:);
triangle_new(36,:) = triangle_new(24,:);
triangle_new(38,:) = triangle_new(37,:) + t/(norm(B_new-C_new))*(B_new-C_new);
triangle_new(40,:) = triangle_new(28,:);
triangle_new(41,:) = triangle_new(43,:);
triangle_new(42,:) = triangle_new(41,:) + t/(norm(C_new-A_new))*(C_new-A_new);

% plotting
figure();
plot_triangle(triangle_new);
hold on;
plot(XYdef(:,1),XYdef(:,2),'-o','Color',[1 0 0],'DisplayName','Deformed');
hold off;

end
