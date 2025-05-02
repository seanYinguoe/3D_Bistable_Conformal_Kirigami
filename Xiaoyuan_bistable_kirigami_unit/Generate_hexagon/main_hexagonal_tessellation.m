% Main program for generating hexagonal kirigami tessellation

% Reset the environment
clc;
clear;

% Set the number of hexagonal unit
num_x = 3; 
num_y = 3;

% Set the parameters of a unit
l1 = 1.03;
l2 = 0.2;
l3 = 0.53;
t = 0.03; 

% Set the coulour of display
colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255}; % The colour of void, flank, filament, Innertriangle

% Set the initial guess of alpha_1 and alpha_2 regarding the displacement delta
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;
delta = 0;

% Generate the triangle
[triangle,alpha_1,alpha_2] = triangle_unit(prev_alpha_1, prev_alpha_2, delta,l1,l2,l3,t);

% Generate the hexagon
[hexagon,~,~] = hexagon_unit(prev_alpha_1, prev_alpha_2, delta,l1,l2,l3,t);

% Generate the hexagon tessellation
tessellation = hexagonal_tessellation(num_x,num_y,l1,l2,l3,t);


%% Plot the result
% Plot the single triangular unit
figure(1)
plot_triangle(triangle,colour) 

% Plot the single hexagon unit
figure(2)
plot_hexagon(hexagon,colour)

% Plot the interactive hexagon unit
interactive_hexagon(l1,l2,l3,t)

% Plot hexagon tessellation
for i = 1:num_x
    for j = 1:num_y
        plot_hexagon(tessellation{i,j},colour);
    end
end

