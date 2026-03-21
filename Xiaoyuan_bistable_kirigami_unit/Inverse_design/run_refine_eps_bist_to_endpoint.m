%% Run endpoint-consistent eps_bist refinement on anisotropy_filter table
clear; clc;

load('anisotropy_filter.mat', 'anisotropy_filter');

opts = struct();
opts.edgeLen = 15;
opts.l1 = opts.edgeLen * 0.85;
opts.l4 = opts.edgeLen * 0.05;
opts.t = opts.edgeLen * 0.015;
opts.etaThreshold = 0.03;
opts.alphaTarget = 0.99;
opts.maxIter = 6;
opts.epsFloor = 1e-4;
opts.epsCeil = 2.0;
opts.alphaTol = 1e-3;

anisotropy_filter_corr = refine_eps_bist_to_endpoint(anisotropy_filter, opts);
% 
save('anisotropy_filter_corr.mat', 'anisotropy_filter_corr', 'opts');
fprintf('Saved: anisotropy_filter_corr.mat\n');

% ---- Plot energy curve for the corrected row using deform_triangle_anisotropic ----
nD = 150;
Nseg = 8;

row = anisotropy_filter_corr(1,:);

a1 = row.a1;
a2 = row.a2;
a3 = row.a3;
beta_row = row.beta;

% Use refined eps if available, otherwise fall back to original
if ismember('eps_bist_new', row.Properties.VariableNames) && isfinite(row.eps_bist_new)
    eps_use = row.eps_bist_new;
else
    eps_use = row.eps_bist;
end

lam3 = 1 + eps_use;
lam1 = lam3 * sin(a1)/sin(a3);
lam2 = lam3 * sin(a2)/sin(a3);

L1 = lam1 * opts.edgeLen;   % |q2-q3|
L2 = lam2 * opts.edgeLen;   % |q1-q3|
L3 = lam3 * opts.edgeLen;   % |q1-q2|

q1 = [0, 0, 0];
q2 = [0, -L3, 0];
q3_y = (L1^2 - L2^2 - L3^2) / (2*L3);
q3_x = -sqrt(max(L2^2 - q3_y^2, 0));
q3 = [q3_x, q3_y, 0];

[E_sel, alpha_sel] = deform_triangle_anisotropic( ...
    q1, q2, q3, opts.edgeLen, opts.l1, opts.l4, beta_row, opts.t, nD, Nseg, true);
