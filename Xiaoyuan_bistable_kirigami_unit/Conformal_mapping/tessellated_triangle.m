function triangle_tessellation = tessellated_triangle(f_out, c_out, i_out, params, stretch_facs)
% Find the optimised bistable unit for each triangular unit, and move it to
% cooresponding grid
% edgeLen = l2 + l1 + t + l4 + l4 + delta

%% Generate optimised bistable unit for each triangular grid
% Create standards bistable unit(upwards and downwards)
edgeLen = params(1);
l1 = params(2);
l2 = params(3);
l3 = params(4);
t = params(5);
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;

%% Generate optimised bistable unit for each grid and move bistable unit to responding grid
triangle_tessellation = cell(size(f_out,1),1);
for i = 1:size(f_out,1)
    delta = 0; % default, closed state
    delta = edgeLen * strech_facs(i);
    if i_out(i) == 0 % upwards triangle
        [triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta,l1,l2,l3,t);
        triangle = triangle/strech_facs(i); % rescale the triangle unit to fit the grid
        triangle_upward = triangle*rotation(-pi/6); % triangle upwards 0
        centroid_upward = [0,-edgeLen/2*sqrt(3)*2/3];
        triangle_tessellation{i} = triangle_upward + (c_out(i,1:2) - centroid_upward);
    else
        [triangle,~,~] = triangle_unit(prev_alpha_1, prev_alpha_2, delta,l1,l2,l3,t);
        triangle = triangle/strech_facs(i); % rescale the triangle unit to fit the grid
        triangle_downward = [triangle_upward(:,1),-triangle_upward(:,2)]; % triangle downwards 1
        centroid_downward = [0,edgeLen/2*sqrt(3)*2/3];
        triangle_tessellation{i} = triangle_downward + (c_out(i,1:2) - centroid_downward);
    end
end

end




