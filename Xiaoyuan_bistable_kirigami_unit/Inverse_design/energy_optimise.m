%% inverse design to maximum the energy gap
% Input:
% params: design paramter
% edgeLen = params(1): The total edge of unit, we keep it same
% l1 = params(2): The length of flanks, we keep it same with every unit
% l2 = params(3); The length of ligaments  
% l3 = params(4); The length of triangle
% t = params(5);  The thickness of ligaments(infinite small)
% q1,q2,q3: coordinates of deployed unit

%% Get the closed and deployed configuration
edgeLen = 15;
l1 = edgeLen * 0.8;
l2 = edgeLen * 0.1;
t  = edgeLen * 0.02;
l3 = l1 -2*l2 - (edgeLen - l1 -l2)/2;  % stretch length
q10 = [0,sqrt(3)/2*edgeLen,0]; 
q20 = [-1/2*edgeLen,0,0];
q30 = [1/2*edgeLen,0,0];
triangle_initial = deform_triangle(q10,q20,q30,edgeLen,l1,l2,l3,t,0);
colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255}; % The colour of void, flank, filament, Innertriangle
%plot_triangle(triangle_initial,colour)
q1 = [0,sqrt(3)/2*edgeLen*1.6,0]; 
q2 = [-1/2*edgeLen*1.6,0,0];
q3 = [1/2*edgeLen*1.6,0,0];
triangle_deploy = deform_triangle(q1,q2,q3,edgeLen,l1,l2,l3,t,0);
%plot_triangle(triangle_deploy,colour)


%% Get the energy curve during deployment(uniform deployment)
E_rotation = [];
index = [23,36,37;
    33,37,36;
    19,32,33;
    41,33,32;
    27,40,41;
    37,41,40]; % Index to calculate the rotational energy in springs
k = 1;
k_r = 1;
for delta = 0:0.1:0.55*edgeLen
    q1 = [0,sqrt(3)/2*(edgeLen+delta),0];
    q2 = [-1/2*(edgeLen+delta),0,0];
    q3 = [1/2*(edgeLen+delta),0,0];
    triangle_deploy = deform_triangle(q1,q2,q3,edgeLen,l1,l2,l3,t,0);
    for i = 1:size(index,1)
        angle_diff(i) = angle_calculate(triangle_deploy(index(i,:),:)) - ...
            angle_calculate(triangle_initial(index(i,:),:));
    end
    E_rotation(k) = sum(1/2 * k_r * angle_diff.^2);
    k = k+1;
end

% Plot the results
figure;
hold on;
delta = 0:0.1:0.55*edgeLen;
plot(delta, E_rotation,'black', 'LineWidth', 1.5); % Energy vs Delta
%Add labels, title, and legend
xlabel('Displacement (\delta)', 'FontSize', 18);
ylabel('E/K', 'FontSize', 18);
title('Energy-Displacement Curve', 'FontSize', 18);
legend({'Theory'}, 'Location', 'southeast', 'FontSize', 18);
grid on;
set(gca, 'FontSize', 18);




% Contraints: edgeLen = l1 + 2*l4 + l2  l4 = (edgeLen - l1 - l2)/2