% Caculate the in-plane energy and out-of-plane energy of tessellation
% during deployment
function [E_rotation,E_bending] = energy_calculate(initial_config, target_config, v_out, f_out, x_out)
%% Clear and organise data from mat
index = [23,36,37;
    33,37,36;
    19,32,33;
    41,33,32;
    27,40,41;
    37,41,40]; % Index to calculate the rotational energy in springs

k_r = 1; % Rotational stiffness of a single spring
k_b = 4; % Bending stiffness of a single hinge

%% Calculate the in-plane rotational energy caused by filaments
% Calculate the rotational angle of each triangle
angle_diff = zeros(size(target_config,1),size(index,1));
for i = 1:size(target_config,1)
    for j = 1:size(index,1)
        angle_diff(i,j) = angle_calculate(target_config{i}(index(j,:),:)) - ...
            angle_calculate(initial_config{i}(index(j,:),:));
    end
end

E_rotation = 1/2 * k_r * angle_diff.^2;
%% Calculate the out-of-plane bending energy caused by flanks
% Calculate the dihedral angles between triangles
dihedral_angle = dihedral_angle_calculate(v_out, f_out, x_out);
E_bending = 1/2 * k_b * (dihedral_angle(:,3)).^2;
end
