function triangle = anisotropic_triangle_expression(alpha_1, alpha_2, ...
    beta1, beta2, beta3, edgeLen, l1, l4, t)
% Use A(beta1,beta2,beta3) define a anisotropic triangle
% 2D Rotation
R = @(th)[cos(th) -sin(th); sin(th) cos(th)];

%% Compute geometry of triangle given beta, l1, l4 and t
% returns l2,l3,l5,l6 for this beta (no extra params)
len_on_flank = @(beta) deal( ...
    ... l2
    (2/sqrt(3)) * (edgeLen - l1 - l4) * sin(pi/3 - beta), ...
    ... l3
    local_l3(beta, edgeLen, l1, l4), ...
    ... l5
    ((sqrt(3)/2)*l4 + sin(beta) * ( (2/sqrt(3))*sin(pi/3 - beta)*(l1-0.5*l4) - cos(pi/3 - beta)*l4 )) ...
        / sin(pi/3 - beta), ...
    ... l6
    (2/sqrt(3))*sin(pi/3 - beta)*(l1-0.5*l4) - cos(pi/3 - beta)*l4 );

%% Build up one flank given beta
% returns: void, filament, flank polygons and the 3 inner-triangle nodes on this wedge
    function [void_i, filament_i, flank_i] = build_flank(beta)
    % lengths for this side
    [l2, l3, l5, l6] = len_on_flank(beta);

    % local anchor points for this side (same as your original)
    point_A = [0; -l4];        % fixed point
    point_O = [0;  0];

    % vectors (local frame for this side, using this beta)
    AC_ = [-l6*cos(pi/6+beta); -l6*sin(pi/6+beta)];
    C_D_ = R(-(alpha_1-pi/3)) * [ l2*cos(pi/6-beta); -l2*sin(pi/6-beta) ];
    D_D  = R(-(alpha_2-(pi-pi/3))) * R(-(alpha_1-pi/3)) * [ l3*cos(pi/6+beta);  l3*sin(pi/6+beta) ];
    DC   = R(-(alpha_1-pi/3)) * [ l2*cos(pi/6+beta);  l2*sin(pi/6+beta) ];
    CB   = [ l5*cos(pi/6+beta);  l5*sin(pi/6+beta) ];

    % node chain on this side (like your original)
    point_C_ = point_A + AC_;
    point_D_ = point_C_ + C_D_;
    point_D  = point_D_ + D_D;
    point_C  = point_D  + DC;
    point_B  = point_C  + CB;

    % filament offsets at the two hinges on this side (width t)
    D_E =  (t/l3) * D_D;
    C_F =  (t/l6) * (-AC_);

    point_E = point_D_ + D_E;
    point_F = point_C_ + C_F;

    % void polygon (white cut)
    void_i = [ ...
        point_A.'; point_F.'; point_E.'; point_D.'; point_C.'; point_B.' ...
    ];

    % filament polygon (green)
    filament_i = [ ...
        point_F.'; point_C_.'; point_D_.'; point_E.' ...
    ];

    % flank polygon (red)
    C_B_ = [-l5*sin(pi/3+beta); l5*cos(pi/3+beta)];
    point_B_ = point_C_ + C_B_;
    flank_i = [ ...
        point_A.'; point_C_.'; point_B_.'; point_O.' ...
    ];

    % % three inner-triangle nodes tied to this ligament
    % D_D__ = R(-pi/3) * D_D;      % third node around the hinge
    % point_D__ = point_D_ + D_D__;
    % inner_nodes_i = [point_D.'; point_D_.'; point_D__.' ];  % 3×2 rows
end

%% Build other two flanks by rotation the ref flank around centre
% flank 1 at theta1 = 0
theta = [0, -2*pi/3, +2*pi/3];
betas = [beta1, beta2, beta3];

void_all     = [];
flank_all    = [];
filament_all = [];
inner_all    = [];

for k = 1:3
    beta_k  = betas(k);
    theta_k = theta(k);

    [void_k, fil_k, flk_k] = build_flank(beta_k);
    C = [-sqrt(3)/6*edgeLen, -edgeLen/2];     % unit centroid

    % rotate every vertex to its global orientation
    void_k  = (void_k  - C) * R(theta_k) + C;
    fil_k   = (fil_k   - C) * R(theta_k) + C;
    flk_k   = (flk_k   - C) * R(theta_k) + C;

    void_all     = [void_all;     void_k];
    flank_all    = [flank_all;    flk_k];
    filament_all = [filament_all; fil_k];
    inner_all    = [inner_all;    fil_k(3,:)];
end

% Inner triangle vertices: one vertex from each flank (first row of inner_k)
% Using the first point (point_D) of each flank to define the inner triangle
% P1 = inner_all(1, :);
% P2 = inner_all(4, :);       % 1st row each block of 3
% P3 = inner_all(7, :);

Innertriangle = inner_all;

% Build the whole unit(flank, filament, innertriangle)
triangle = [ ...
    void_all; ...
    flank_all; ...
    filament_all; ...
    Innertriangle ...
];

end

%% Calculate the l3 based on beta,edgeLen, l1 and l4
function l3 = local_l3(beta, edgeLen, l1, l4)
    A  = pi/3 - beta;
    l2 = (2/sqrt(3)) * (edgeLen - l1 - l4) * sin(A);
    l6 = (2/sqrt(3)) * sin(A) * (l1 - 0.5*l4) - cos(A) * l4;
    l5 = ( (sqrt(3)/2) * l4 + sin(beta) * l6 ) / sin(A);
    l3 = l6 - l5 - 1.5*l2 - l2 * ( (sqrt(3)/2) * (cos(A) / sin(A)) );
end
