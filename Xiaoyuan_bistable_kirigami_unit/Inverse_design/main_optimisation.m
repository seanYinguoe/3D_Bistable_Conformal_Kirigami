% main program for inverse design program

% Find optimal beta for each unit based on stretch_facs
opt_beta = unit_design(stretch_facs, T2);

% Tessellate each unit with opsimised beta(undeployed)
opt_tessellation = tessellated_triangle(f_out, i_out, params,v_out,opt_beta);

% Tessellate each unit with opsimised beta(deployed)
opt_tessellation_target = tessellated_triangle(f_out, i_out, params, v_target, opt_beta);

% Create deployment from closed state to open state
tessellation_deployment(opt_tessellation,opt_tessellation_target,T,modelname); % Plot deployment

% Create deployment figure
filename = 'optimised_deployment_progress.gif';
generate_gif(opt_tessellation, opt_tessellation_target,T, modelname, filename)

% Generate energy evolvement during deployment
[E_rotation, E_bending] = deployment_energy(opt_tessellation, opt_tessellation_target, v_out, v_target, f_out, x_out);

%% Plot the tessellated configuration
colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255}; % The colour of void, flank, filament, Innertriangle
figure()
hold on
for i = 1:size(f_out,1)
    plot_triangle(opt_tessellation{i},colour);
end
hold off
axis off