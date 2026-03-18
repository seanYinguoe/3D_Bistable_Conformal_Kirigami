%% Extend beta sweep from existing [0, pi/15] to [0, pi/12] without rerunning old betas
% This script does NOT modify previous functions/files.
% It only computes the missing beta slices and merges them with existing data.

clear; clc;

%% User settings (must match your original ternary_study settings)
edgeLen = 15;
scale   = 1.63;
l1_geom = edgeLen * 0.85;
l4_geom = edgeLen * 0.05;
t_geom  = edgeLen * 0.015;
beta_max_target = pi/12;

inFile = 'anisotropy_study.mat';
if ~exist(inFile, 'file')
    error('Cannot find %s. Run ternary_study first.', inFile);
end

S = load(inFile, 'anisotropy_study');
if ~isfield(S, 'anisotropy_study')
    error('%s does not contain variable anisotropy_study.', inFile);
end
anisotropy_study = S.anisotropy_study;

if isempty(anisotropy_study)
    error('anisotropy_study is empty.');
end

%% Infer existing beta grid and step
beta_existing = unique(anisotropy_study.beta, 'stable');
if numel(beta_existing) < 2
    error('Need at least 2 existing beta values to infer interval.');
end
dbeta = median(diff(beta_existing));
tolb = max(1e-12, abs(dbeta)*1e-6);

% Build new beta list using same interval
beta_new = (beta_existing(end) + dbeta) : dbeta : (beta_max_target + tolb);
beta_new = beta_new(beta_new > beta_existing(end) + tolb);
beta_new = beta_new(beta_new <= beta_max_target + tolb);

if isempty(beta_new)
    fprintf('No new beta to compute. Existing beta range already reaches %.6f rad.\n', beta_existing(end));
    return;
end

fprintf('Existing beta max: %.6f rad\n', beta_existing(end));
fprintf('Target beta max:   %.6f rad\n', beta_max_target);
fprintf('Interval dbeta:    %.6f rad\n', dbeta);
fprintf('New beta count:    %d\n', numel(beta_new));

%% Recover the (a1,a2,a3) configuration list from one existing beta slice
beta0 = beta_existing(1);
cfg_mask = abs(anisotropy_study.beta - beta0) < tolb;
cfg_tbl = anisotropy_study(cfg_mask, {'a1','a2','a3'});
cfg_tbl = unique(cfg_tbl, 'rows', 'stable');

a1 = cfg_tbl.a1;
a2 = cfg_tbl.a2;
a3 = cfg_tbl.a3;
nConfig = height(cfg_tbl);
nBetaNew = numel(beta_new);

% Precompute geometry terms independent of beta
sin_a1 = sin(a1);
sin_a2 = sin(a2);
sin_a3 = sin(a3);

lam1_ratio = sin_a1 ./ sin_a3;
lam2_ratio = sin_a2 ./ sin_a3;
lam3_ratio = ones(nConfig,1);

lambda1 = lam1_ratio * scale;
lambda2 = lam2_ratio * scale;
lambda3 = lam3_ratio * scale;

L1 = lambda1 * edgeLen;
L2 = lambda2 * edgeLen;
L3 = lambda3 * edgeLen;

q3_y = (L1.^2 - L2.^2 - L3.^2) ./ (2 * L3);
inside = L2.^2 - q3_y.^2;
valid = inside > 0;
q3_x = NaN(nConfig,1);
q3_x(valid) = -sqrt(inside(valid));

%% Compute only new beta slices (with checkpoints)
outDir = 'beta_batches_extend_pi12';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

eps_new = NaN(nConfig*nBetaNew,1);
eta_new = NaN(nConfig*nBetaNew,1);
beta_all_new = kron(beta_new(:), ones(nConfig,1));

for ib = 1:nBetaNew
    beta = beta_new(ib);
    batchFile = fullfile(outDir, sprintf('beta_ext_%03d.mat', ib));

    if exist(batchFile, 'file')
        B = load(batchFile, 'eps_beta', 'eta_beta');
        eps_beta = B.eps_beta;
        eta_beta = B.eta_beta;
        fprintf('Loaded checkpoint for new beta %d/%d.\n', ib, nBetaNew);
    else
        eps_beta = NaN(nConfig,1);
        eta_beta = NaN(nConfig,1);

        parfor i = 1:nConfig
            if ~valid(i)
                continue;
            end

            q1 = [0, 0];
            q2 = [0, -L3(i)];
            q3 = [q3_x(i), q3_y(i)];

            [eps_bist, eta_val] = bistability_analysis( ...
                l1_geom, l4_geom, beta, t_geom, edgeLen, q1, q2, q3, 0);

            eps_beta(i,1) = eps_bist;
            eta_beta(i,1) = eta_val;
        end

        save(batchFile, 'ib', 'beta', 'eps_beta', 'eta_beta');
        fprintf('Saved checkpoint for new beta %d/%d.\n', ib, nBetaNew);
    end

    idx0 = (ib-1)*nConfig + 1;
    idx1 = ib*nConfig;
    eps_new(idx0:idx1,1) = eps_beta;
    eta_new(idx0:idx1,1) = eta_beta;
end

%% Build new table and merge
a1_new = repmat(a1, nBetaNew, 1);
a2_new = repmat(a2, nBetaNew, 1);
a3_new = repmat(a3, nBetaNew, 1);

anisotropy_study_new = table( ...
    a1_new, a2_new, a3_new, eps_new, eta_new, beta_all_new, ...
    'VariableNames', {'a1','a2','a3','eps_bist','eta_val','beta'});

anisotropy_study_pi12 = [anisotropy_study; anisotropy_study_new];
anisotropy_study_pi12 = sortrows(anisotropy_study_pi12, {'beta','a1','a2','a3'});

save('anisotropy_study_pi12.mat', 'anisotropy_study_pi12', 'beta_existing', 'beta_new', 'dbeta');
fprintf('Saved merged table: anisotropy_study_pi12.mat\n');

%% Optional: regenerate filled table
[anisotropy_clean, reportClean] = clean_isolated_outliers(anisotropy_study_pi12); %#ok<ASGLU>
[anisotropy_filled_pi12, reportFill] = clean_and_fill_anisotropy(anisotropy_clean); %#ok<ASGLU>
save('anisotropy_filled_pi12.mat', 'anisotropy_filled_pi12');
fprintf('Saved filled table: anisotropy_filled_pi12.mat\n');

