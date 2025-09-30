function [strain, U_total] = hbm_energy_displacement(disp,l1,l4,t,beta,edgeLen)
%% Parameters
n_step = 200; % Number of interpolation steps
U_total = nan(1, n_step);          % Total energy
delta = zeros(1,n_step);
strain = zeros(1,n_step);

%% Loop through interpolation steps
for i = 1:n_step
    % Compute displacement for each step
    delta(i) = 0 + disp/(n_step-1)*(i-1);
    strain(i) = delta(i)/edgeLen;
    % Compute energy for current configuration
    try
        [~, E] = deform_triangle_semi(delta(i),edgeLen,l1,l4,beta,t);
        U_total(i) = E;
    catch ME
        warning('Step %.3f failed: %s', delta(i), ME.message);
        U_total(i) = NaN;
    end
end
