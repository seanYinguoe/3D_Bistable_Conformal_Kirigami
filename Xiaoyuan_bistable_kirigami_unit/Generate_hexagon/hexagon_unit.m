% Create a hexogon unit from a triangular unit
% Input: deploy angle bistable unit: alpha_1, alpha_2
% Output: deformed hexogon unit

%% Create a triangular unit
function [hexagon,alpha_1_optimal,alpha_2_optimal] = hexagon_unit(prev_alpha_1, prev_alpha_2, delta, l1, l2, l3, t)

[triangle_unit_1,alpha_1_optimal,alpha_2_optimal] = triangle_unit(prev_alpha_1, prev_alpha_2, delta, l1, l2, l3,t);

% Mirror a triangular unit to a hexogon unit
triangle_unit_3 = triangle_unit_1 * rotation(-2*pi/3);
triangle_unit_4 = [-triangle_unit_3(:,1),triangle_unit_3(:,2)];
triangle_unit_6 = [-triangle_unit_1(:,1),triangle_unit_1(:,2)];
triangle_unit_2 = triangle_unit_6 * rotation(-2*pi/3);
triangle_unit_5 = [-triangle_unit_2(:,1),triangle_unit_2(:,2)];
hexagon = {triangle_unit_1,...
    triangle_unit_2,...
    triangle_unit_3,...
    triangle_unit_4,...
    triangle_unit_5,...
    triangle_unit_6};
end
