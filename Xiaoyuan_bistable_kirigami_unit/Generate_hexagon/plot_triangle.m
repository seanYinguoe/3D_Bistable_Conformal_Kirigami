% Visulize the filament
%          the void
%          the flank
%          the Innertriangle
% Input parameters:
%          Type: triangle:filament/void/flank/Innertriangle
%          void: 6*2*3; flank:4*2*3 filament:4*2*3 Innertriangle:3*2*1
%          Colour: Colour of different part: void, flank, filament,
%          Innertriangle
function plot_triangle(triangle,colour)
if nargin < 2
    colour = {'white', 'black', 'black', 'black'};
end
axis equal;
hold on;
%% Define connecivity of filaments, voids, flanks, Innertriangle
f_void = [1  2  3  4  5  6;
 7  8  9 10 11 12;
13 14 15 16 17 18];
f_flank = [19 20 21 22;
23 24 25 26;
27 28 29 30];
f_filament = [31 32 33 34;
35 36 37 38;
39 40 41 42];
f_inner = [43 44 45];

%% Plot faces
% Plot voids (transparent fill, black edge traces the cut lines)
patch('Vertices', triangle, 'Faces', f_void, ...
    'FaceColor', 'none', 'EdgeColor', 'none');

%'FaceColor', 'none', 'EdgeColor', 'black', 'LineWidth', 0.8);
% Plot flanks
patch('Vertices', triangle, 'Faces', f_flank, ...
    'FaceColor', colour{2}, 'FaceAlpha', 0.5, ...
    'EdgeColor', 'none');

%'EdgeColor', 'none'

% Plot filaments
patch('Vertices', triangle, 'Faces', f_filament, ...
    'FaceColor', colour{3}, 'FaceAlpha', 0.5, ...
    'EdgeColor', 'none');

% Plot innertriangle
patch('Vertices', triangle, 'Faces', f_inner, ...
    'FaceColor', colour{4}, 'FaceAlpha', 0.5, ...
    'EdgeColor', 'none');

% % Plot voids
% for i = 1:3
%     fill([triangle(6*(i-1)+1:6*i,1);triangle(6*(i-1)+1,1)], [triangle(6*(i-1)+1:6*i,2);triangle(6*(i-1)+1,2)], colour{1}, 'FaceAlpha', 0.5,'EdgeColor', 'none');
% end
% % Plot flanks
% for i = 1:3
%     plot([triangle(4*(i-1)+1+18:4*i+18,1);triangle(4*(i-1)+1+18,1)], [triangle(4*(i-1)+1+18:4*i+18,2);triangle(4*(i-1)+1+18,2)], 'Color', colour{2}, 'LineWidth', 1.5, 'MarkerFaceColor', 'cyan');
%     fill([triangle(4*(i-1)+1+18:4*i+18,1);triangle(4*(i-1)+1+18,1)], [triangle(4*(i-1)+1+18:4*i+18,2);triangle(4*(i-1)+1+18,2)], colour{2}, 'FaceAlpha', 0.5);
% end
% % Plot filaments
% for i = 1:3
%     plot([triangle(4*(i-1)+1+30:4*i+30,1);triangle(4*(i-1)+1+30,1)], [triangle(4*(i-1)+1+30:4*i+30,2);triangle(4*(i-1)+1+30,2)], 'Color', colour{3}, 'LineWidth', 1.5, 'MarkerFaceColor', 'cyan');
%     fill([triangle(4*(i-1)+1+30:4*i+30,1);triangle(4*(i-1)+1+30,1)], [triangle(4*(i-1)+1+30:4*i+30,2);triangle(4*(i-1)+1+30,2)], colour{3}, 'FaceAlpha', 0.5);
% end
% % Plot Innertriangle
% plot([triangle(43:45,1);triangle(43,1)], [triangle(43:45,2);triangle(43,2)], 'Color', colour{4}, 'LineWidth', 1.5, 'MarkerFaceColor', 'cyan');
% fill([triangle(43:45,1);triangle(43,1)], [triangle(43:45,2);triangle(43,2)], colour{4}, 'FaceAlpha', 0.5);
% end
