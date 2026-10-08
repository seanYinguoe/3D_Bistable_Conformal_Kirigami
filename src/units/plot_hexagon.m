%% Visulise the geometry(static)
function plot_hexagon(hexagon,colour)
    axis equal;
    hold on;
    for i = 1:6
        triangle = hexagon{i};
        plot_triangle(triangle,colour);
    end
end
