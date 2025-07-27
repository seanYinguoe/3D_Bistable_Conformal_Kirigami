function triangle_tessellation = tessellated_triangle(f_out, i_out, params, v)
% Find the optimised bistable unit for each triangular unit, and move it to
% cooresponding grid
% edgeLen = l1 + 2*l4 + l2

%% Generate optimised bistable unit for each triangular grid
% Create standards bistable unit(upwards and downwards)
edgeLen = params(1);
l1 = params(2);
l2 = params(3);
l3 = params(4);
t = params(5);

%% Generate optimised bistable unit for each grid and move bistable unit to responding grid
triangle_tessellation = cell(size(f_out,1),1);
for i = 1:size(f_out,1)
    % Calculate the edge1, edge2, and edge3 of each triangular grid
    q1 = v(f_out(i,1),:);
    q2 = v(f_out(i,2),:);
    q3 = v(f_out(i,3),:);
    if i_out(i) == 0 % upwards triangle
        [triangle,~,~] = deform_triangle(q3,q1,q2,edgeLen,l1,l2,l3,t,i_out(i));
    else % downwards triangle
        [triangle,~,~] = deform_triangle(q1,q3,q2,edgeLen,l1,l2,l3,t,i_out(i));
    end
    triangle_tessellation{i} = triangle;
end
end




