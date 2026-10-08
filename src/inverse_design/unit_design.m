function opt_beta = unit_design(stretch_facs, T)
% Find optimized beta for given stretch_facs
% stretch_facs : Nx1 double (target stretch factors)
% T            : Sensitivity study with beta and strain_bist
% opt_beta     : Nx1 double (optimized beta values)

% Target strains from stretch_facs
target_strain = stretch_facs - 1;

% Remove NaNs from sensitivity results
valid_idx = ~isnan(T.strain_bist);
beta_vals = T.beta(valid_idx);
strain_vals = T.strain_bist(valid_idx);

% Interpolation
opt_beta = interp1(strain_vals, beta_vals, target_strain, 'linear', 'extrap');

% Replace negative values with 0
opt_beta(opt_beta < 0) = 0;

end