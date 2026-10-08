% Generate a hexagonal tessellation.
% Input: (m,n) the number of units in x and y direction.
%        (l1,l2,l3,t) the geometric parameters of a hexagonal unit
% Output: The hexagonal tessellation

function tessellation = hexagonal_tessellation(num_x,num_y,beta,edgeLen,l1,l4,t,delta)

% Set the number of hexagonal unit
tessellation = cell(num_x,num_y);

if nargin < 8
    delta = 0;
end

% Define the initial guess
prev_alpha_1 = pi/3;
prev_alpha_2 = 2*pi/3;

% Generate single hexagon unit
[hexagon,~,~] = hexagon_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);
%plot_hexagon(hexagon_unit,colour) % Plot the single hexagon unit

% Create a hexagon tessellation
for i = 1:num_x
    for j = 1:num_y
        for k = 1:6
            if mod(j,2) == 1
                tessellation{i,j}{k}(:,1) = hexagon{k}(:,1) + (edgeLen+delta)*sqrt(3)/2*(i-1)*2;
                tessellation{i,j}{k}(:,2) = hexagon{k}(:,2) + (edgeLen+delta)*(3/2)*(j-1);
            else
                tessellation{i,j}{k}(:,1) = hexagon{k}(:,1) + (edgeLen+delta)*sqrt(3)/2*(i-1)*2-(edgeLen+delta)*sqrt(3)/2;
                tessellation{i,j}{k}(:,2) = hexagon{k}(:,2) + (edgeLen+delta)*(3/2)*(j-1);
            end 
        end
    end
end

% % Plot the hexagon tessellation
% for i = 1:num_x
%     for j = 1:num_y
%         plot_hexagon(hexagon{i,j},colour);
%     end
% end

end