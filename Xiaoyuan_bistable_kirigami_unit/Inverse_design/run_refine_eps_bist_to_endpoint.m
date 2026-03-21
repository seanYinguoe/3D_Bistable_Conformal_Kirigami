%% Run endpoint-consistent eps_bist refinement on existing anisotropy table
clear; clc;

load('anisotropy_study_pi12.mat', 'anisotropy_study_pi12');

opts = struct();
opts.edgeLen = 15;
opts.l1 = opts.edgeLen * 0.85;
opts.l4 = opts.edgeLen * 0.05;
opts.t = opts.edgeLen * 0.015;
opts.etaThreshold = 0.03;
opts.alphaTarget = 0.95;
opts.alphaCloseTol = 0.02;
opts.maxIter = 6;
opts.relax = 0.8;
opts.epsFloor = 1e-4;
opts.epsCeil = 2.0;
opts.alphaTol = 1e-3;

anisotropy_study_pi12_corr = refine_eps_bist_to_endpoint(anisotropy_study_pi12, opts);

save('anisotropy_study_pi12_corr.mat', 'anisotropy_study_pi12_corr', 'opts');
fprintf('Saved: anisotropy_study_pi12_corr.mat\n');
