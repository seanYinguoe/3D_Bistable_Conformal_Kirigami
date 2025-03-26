function triangle_tessellation = tessellated_triangle(f_out, c_out, i_out, edgeLen)
% Find the optimised bistable unit for each triangular unit, and move it to
% cooresponding grid
% edgeLen = l2 + l1 + t + l4 + l4

%% Generate optimised bistable unit for each triangular grid
% Create standards bistable unit(upwards and downwards)
delta = 0; % deployment displacement
l1 = edgeLen * 0.65;
l2 = edgeLen * 0.15;
t  = edgeLen * 0.05;
l3 = l1 -2*l2 - (edgeLen - l2 -l1 -t)/2;
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;
[triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta,l1,l2,l3,t);
triangle_upward = triangle*rotation(-pi/6); % triangle upwards 0
centroid_upward = [0,-edgeLen/2*sqrt(3)*2/3];
triangle_downward = [triangle_upward(:,1),-triangle_upward(:,2)]; % triangle downwards 1
centroid_downward = [0,edgeLen/2*sqrt(3)*2/3];

%% Move bistable unit to responding grid
triangle_tessellation = cell(size(f_out,1),1);
for i = 1:size(f_out,1)
    if i_out(i) == 0 % upwards triangle
        triangle_tessellation{i} = triangle_upward + (c_out(i,1:2) - centroid_upward);
    else
        triangle_tessellation{i} = triangle_downward + (c_out(i,1:2) - centroid_downward);
    end
end

end




