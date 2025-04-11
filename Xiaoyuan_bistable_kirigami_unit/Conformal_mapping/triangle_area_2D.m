function area_2D = triangle_area_2D(v, f)
% Calculate the triangular mesh in 2D
    area_2D = zeros(size(f, 1), 1);
    for i = 1:size(f, 1)
        a = v(f(i,1),:);
        b = v(f(i,2),:);
        c = v(f(i,3),:);
        ab = b - a;
        ac = c - a;
        area2 = abs(ab(1)*ac(2) - ab(2)*ac(1));
        area_2D(i) = area2 / 2;
    end
end