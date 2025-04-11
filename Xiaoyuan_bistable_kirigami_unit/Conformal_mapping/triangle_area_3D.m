function area_3D = triangle_area_3D(v, f)
% Calculate the triangular mesh in 3D
    area_3D = zeros(size(f, 1), 1);
    for i = 1:size(f, 1)
        a = v(f(i,1),:);
        b = v(f(i,2),:);
        c = v(f(i,3),:);
        ab = b - a;
        ac = c - a;
        cross_product = cross(ab, ac);
        area2 = norm(cross_product);
        area_3D(i) = area2 / 2;
    end
end