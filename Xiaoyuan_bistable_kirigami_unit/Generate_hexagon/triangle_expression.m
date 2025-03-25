% Using alpha_1, alpha_2 and delta to derive the geometry information
function triangle = triangle_expression(alpha_1,alpha_2,l1,l2,l3,t)
% Set the size of a unit
% Check the reference on the graph
% l1 = 2*l2+l3+l4
% The total length of a bisatble unit: L_ = l4+l2+l1+l4+t
l4 = l1 - 2*l2 - l3;
theta = pi/3;
%L = l1 + l4*cos(theta)*2 + t;

%% Defining the geometry of a unit
% Defining the geometry of a void(cut)
AF = [-l1*sin(theta),-l1*cos(theta)]';
FE = rotation(-(alpha_1-theta))*[l2*sin(theta),-l2*cos(theta)]';
ED = rotation(-(alpha_2-(pi-theta)))*rotation(-(alpha_1-theta))*[l3*sin(theta),l3*cos(theta)]';
DC = rotation(theta)*FE;
CB = [l4*sin(theta),l4*cos(theta)]';

point_A = [0,-l4]';  % fixed point, displacement control
point_O = [0,0]';
point_F = point_A + AF;
point_E = point_F + FE;
point_D = point_E + ED;
point_C = point_D + DC;
point_B = point_C + CB; % the x coordinate of B should be 0

void_1x = [point_A(1), point_F(1), point_E(1), point_D(1), point_C(1), point_B(1)];
void_1y = [point_A(2), point_F(2), point_E(2), point_D(2), point_C(2), point_B(2)];
void_1 = [void_1x;void_1y]'; % Coordinates of a void


% Defining the geometry of a filament
D_E = t/l3 * ED;
C_F = t/l1 * (-AF);

point_D_ = point_E - D_E;
point_C_ = point_F - C_F;

filament_1x = [point_F(1),point_C_(1),point_D_(1),point_E(1)];
filament_1y = [point_F(2),point_C_(2),point_D_(2),point_E(2)];
filament_1 = [filament_1x;filament_1y]';

% Defining the geometry of a flank
C_B_ = [-l4*sin(theta),l4*cos(theta)]';
point_B_ = point_C_ +C_B_;

flank_1x = [point_A(1),point_C_(1),point_B_(1),point_O(1)];
flank_1y = [point_A(2),point_C_(2),point_B_(2),point_O(2)];
flank_1 = [flank_1x;flank_1y]';

% Defining the geometry of the inner triangle
D_D = point_D - point_D_;
D_D__ = rotation(-pi/3) * D_D;
point_D__ = point_D_ + D_D__;

Innertriangle_x = [point_D(1),point_D_(1),point_D__(1)];
Innertriangle_y = [point_D(2),point_D_(2),point_D__(2)];
Innertriangle = [Innertriangle_x;Innertriangle_y]';

% Mirror to other parts to create a triangular unit
Rotation_centre = 1/3 * sum(Innertriangle);

flank_2 = (flank_1 - Rotation_centre)*rotation(-2*pi/3) + Rotation_centre;
flank_3 = (flank_1 - Rotation_centre)*rotation(2*pi/3) + Rotation_centre;

void_2 = (void_1 - Rotation_centre)*rotation(-2*pi/3) + Rotation_centre;
void_3 = (void_1 - Rotation_centre)*rotation(2*pi/3) + Rotation_centre;

filament_2 = (filament_1 - Rotation_centre)*rotation(-2*pi/3) + Rotation_centre;
filament_3 = (filament_1 - Rotation_centre)*rotation(2*pi/3) + Rotation_centre;
triangle = [void_1;void_2;void_3;flank_1;flank_2;flank_3;filament_1;filament_2;filament_3;Innertriangle];

end