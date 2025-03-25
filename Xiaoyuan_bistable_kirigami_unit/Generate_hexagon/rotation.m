% Defining rotational matrix function
function R = rotation(theta)
    % rotation2D generates a rotation matrix for a given angle.
    % Input:
    %   angle - angle in degrees for the rotation
    % Output:
    %   R - 2x2 rotation matrix

    % Convert angle from degrees to radians

    % Define the 2D rotation matrix
    R = [cos(theta), -sin(theta);
         sin(theta), cos(theta)];
end