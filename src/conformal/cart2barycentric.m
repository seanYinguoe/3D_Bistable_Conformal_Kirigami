function uv = cart2barycentric(tri,p)
% Barycentric coordinate calculation
v0 = tri(2,:) - tri(1,:);
v1 = tri(3,:) - tri(1,:);
v2 = p - tri(1,:);

d00 = dot(v0, v0);
d01 = dot(v0, v1);
d11 = dot(v1, v1);
d20 = dot(v2, v0);
d21 = dot(v2, v1);

denom = d00 * d11 - d01 * d01;
v = (d11 * d20 - d01 * d21) / denom;
w = (d00 * d21 - d01 * d20) / denom;
u = 1 - v - w;

uv = [u, v, w];
end