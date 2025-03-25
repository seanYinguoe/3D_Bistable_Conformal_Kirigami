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
% Plot voids
for i = 1:3
    %plot([triangle(6*(i-1)+1:6*i,1);triangle(6*(i-1)+1:6*i,1)], [triangle(6*(i-1)+1:6*i,2);triangle(6*(i-1)+1:6*i,2)], colour{1}, 'LineWidth', 0.01, 'MarkerFaceColor', 'cyan');
    fill([triangle(6*(i-1)+1:6*i,1);triangle(6*(i-1)+1,1)], [triangle(6*(i-1)+1:6*i,2);triangle(6*(i-1)+1,2)], colour{1}, 'FaceAlpha', 0.5,'EdgeColor', 'none');
end
% Plot flanks
for i = 1:3
    plot([triangle(4*(i-1)+1+18:4*i+18,1);triangle(4*(i-1)+1+18,1)], [triangle(4*(i-1)+1+18:4*i+18,2);triangle(4*(i-1)+1+18,2)], 'Color', colour{2}, 'LineWidth', 1.5, 'MarkerFaceColor', 'cyan');
    fill([triangle(4*(i-1)+1+18:4*i+18,1);triangle(4*(i-1)+1+18,1)], [triangle(4*(i-1)+1+18:4*i+18,2);triangle(4*(i-1)+1+18,2)], colour{2}, 'FaceAlpha', 0.5);
end
% Plot filaments 
for i = 1:3 
    plot([triangle(4*(i-1)+1+30:4*i+30,1);triangle(4*(i-1)+1+30,1)], [triangle(4*(i-1)+1+30:4*i+30,2);triangle(4*(i-1)+1+30,2)], 'Color', colour{3}, 'LineWidth', 1.5, 'MarkerFaceColor', 'cyan');
    fill([triangle(4*(i-1)+1+30:4*i+30,1);triangle(4*(i-1)+1+30,1)], [triangle(4*(i-1)+1+30:4*i+30,2);triangle(4*(i-1)+1+30,2)], colour{3}, 'FaceAlpha', 0.5);
end
% Plot Innertriangle 
plot([triangle(43:45,1);triangle(43,1)], [triangle(43:45,2);triangle(43,2)], 'Color', colour{4}, 'LineWidth', 1.5, 'MarkerFaceColor', 'cyan');
fill([triangle(43:45,1);triangle(43,1)], [triangle(43:45,2);triangle(43,2)], colour{4}, 'FaceAlpha', 0.5);

end