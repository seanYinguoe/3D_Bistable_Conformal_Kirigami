% Main program for generating triangle-based kirigami tessellation

% Reset the environment
clc;
clear;

% Set the number of hexagonal unit
num_x = 2; 
num_y = 4;

% Set the parameters of a unit(fixed)
edgeLen = 15; % length of a unit
l4 = 0.05 * edgeLen; % thickness of flank
l1 = 0.80 * edgeLen; % length of flank(not correct)

% Set the parameters of a unit(variable)
beta = 0; % titling angle
t  = edgeLen * 0.032; % thickness of filaments

% Set the coulour of display
%colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255}; % The colour of void, flank, filament, Innertriangle
colour = {'white', [0.9216    0.8863    0.4235
], [ 0.7059    0.9608    0.4118], [0.9216    0.8863    0.4235]};
% Set the initial guess of alpha_1 and alpha_2 regarding the displacement delta
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;
delta = 0;

% Generate the triangle
[triangle,alpha_1,alpha_2] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);

% Generate the hexagon
[hexagon,~,~] = hexagon_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);

% Generate the triangle tessellation
tessellation = triangle_tessellation(num_x, num_y, beta, edgeLen, l1, l4, t);


%% Plot the result
% Plot the single triangular unit
figure(1)
plot_triangle(triangle,colour) 

% Plot the single hexagon unit
figure(2)
plot_hexagon(hexagon,colour)

% Plot the interactive hexagon unit
interactive_triangle(edgeLen,l1,l4,beta,t)

% Plot the interactive hexagon unit
interactive_hexagon(edgeLen,l1,l4,beta,t)

% Plot triangle tessellation
figure(3)
hold on
for i = 1:numel(tessellation)
    plot_triangle(tessellation{i}, colour);
end
hold off
axis equal
axis off

% Create svg cut pattern for fabrication
generate_svg(tessellation, 'triangle_tessellation.svg', false, [], 0.30);
generate_svg(hexagon, 'hexagon.svg', false, [], 0.30);