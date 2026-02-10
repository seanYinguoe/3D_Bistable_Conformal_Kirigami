function [flattened_surface, deployed_surface] = model_rotate(modelName, flattened_surface, deployed_surface)
% Apply model-specific rotations to 2D and 3D surfaces based on modelName
% Inputs:
%   - modelName: string, name of the geometric model
%   - flattened_surface: Nx3 array (2D geometry in 3D coordinates)
%   - deployed_surface: Nx3 array (3D geometry)
% Output:
%   - rotated versions of both surfaces

switch modelName
    case 'double_dome'
        flattened_surface = rotate_pts(flattened_surface, 180, 'x');
        deployed_surface = rotate_pts(deployed_surface, -90, 'x');
        deployed_surface = rotate_pts(deployed_surface, 180, 'x');

    case 'hypar'
        flattened_surface = rotate_pts(flattened_surface, 180, 'x');
        flattened_surface = rotate_pts(flattened_surface, 180, 'y');
        deployed_surface = rotate_pts(deployed_surface, 90, 'x');

    case 'hemisphere'
        deployed_surface = rotate_pts(deployed_surface, 90, 'x');

    case 'arbitrary_surface'
        deployed_surface = rotate_pts(deployed_surface, 90, 'x');

    case 'lilium_tower'
        deployed_surface = rotate_pts(deployed_surface, 90, 'x');
        flattened_surface = rotate_pts(flattened_surface, 45, 'z');

    case 'quarter_dome'
        flattened_surface = rotate_pts(flattened_surface, 180, 'x');
        flattened_surface = rotate_pts(flattened_surface, 180, 'y');
        flattened_surface = rotate_pts(flattened_surface, -35, 'z');
        deployed_surface = rotate_pts(deployed_surface, 90, 'x');

    case 'squidward'
        deployed_surface = rotate_pts(deployed_surface, 90, 'x');

    case 'tigridia'
        flattened_surface = rotate_pts(flattened_surface, -30, 'z');

    case 'heart'
        flattened_surface = rotate_pts(flattened_surface, 90, 'x');
        flattened_surface = rotate_pts(flattened_surface, 90, 'y');
        deployed_surface = rotate_pts(deployed_surface, 90, 'x');
end
    function pts_rotated = rotate_pts(pts, angle_deg, axis)
        angle_rad = deg2rad(angle_deg);
        c = cos(angle_rad);
        s = sin(angle_rad);
        switch lower(axis)
            case 'x'
                R = [1 0 0; 0 c -s; 0 s c];
            case 'y'
                R = [c 0 s; 0 1 0; -s 0 c];
            case 'z'
                R = [c -s 0; s c 0; 0 0 1];
            otherwise
                error('Invalid axis');
        end
        pts_rotated = (R * pts')'; % rotate each row
    end
end
