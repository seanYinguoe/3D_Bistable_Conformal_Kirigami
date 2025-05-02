function def_factor = def_factor(v_out,v_target,f_out)
% Calculate the stretch_factor based on deformed edge and undeformed edge
% Input:
% v_out: the coordinates of vertices of undeformed grid configuration
% v_target: the coordinates of vertices of deformed grid configuration
% f_out: connectivity of vertices to form a triangle
% Output:
% stretch_factor: ratio between deformed edge over undeformed edge

stretch_factor = zeros(size(f_out));

for i = 1:size(f_out,1)
    index = f_out(i,:);
    edge_undeformed1 = norm(v_out(index(1),:)-v_out(index(2),:));
    edge_deforrmed1 = norm(v_target(index(1),:)-v_target(index(2),:));
    stretch_facs1 = edge_deforrmed1/edge_undeformed1;

    edge_undeformed2 = norm(v_out(index(2),:)-v_out(index(3),:));
    edge_deforrmed2 = norm(v_target(index(2),:)-v_target(index(3),:));
    stretch_facs2 = edge_deforrmed2/edge_undeformed2;

    edge_undeformed3 = norm(v_out(index(1),:)-v_out(index(3),:));
    edge_deforrmed3 = norm(v_target(index(1),:)-v_target(index(3),:));
    stretch_facs3 = edge_deforrmed3/edge_undeformed3;

    stretch_factor(i,:) = [stretch_facs1 stretch_facs2 stretch_facs3];
    def_factor(i,:) = [sqrt(stretch_facs1*stretch_facs3)
        stretch_facs1*stretch_facs2
        stretch_facs2*stretch_facs3];
end
