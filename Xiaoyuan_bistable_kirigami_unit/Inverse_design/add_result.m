%% ---- Patch: add isotropic point (a1=a2=a3=pi/3) for all beta ----
a_iso = pi/3;

% geometry constants (same as sweep)
l1_geom = edgeLen * 0.85;
l4_geom = edgeLen * 0.05;
t_geom  = edgeLen * 0.015;

% Build triangle for isotropic angles:
% In your mapping, lambda1=lambda2=lambda3=scale when sin(a1)=sin(a2)=sin(a3).
L1_iso = scale * edgeLen;
L2_iso = scale * edgeLen;
L3_iso = scale * edgeLen;

q1 = [0, 0];
q2 = [0, -L3_iso];
q3_y = (L1_iso^2 - L2_iso^2 - L3_iso^2) / (2 * L3_iso);      % = -L3/2
inside = L2_iso^2 - q3_y^2;                                  % = 3/4 L3^2
q3_x = -sqrt(max(inside, 0));
q3 = [q3_x, q3_y];

% Run bistability for each beta (no parfor needed; only 20 calls)
nBeta = numel(beta_vec);
eps_iso = NaN(nBeta,1);
eta_iso = NaN(nBeta,1);

for ib = 1:nBeta
    beta = beta_vec(ib);

    % IMPORTANT: set do_plot = 0 here, otherwise you will pop up 20 figures
    [eps_bist, eta_val] = bistability_analysis( ...
        l1_geom, l4_geom, beta, t_geom, edgeLen, ...
        q1, q2, q3, 0);

    eps_iso(ib) = eps_bist;
    eta_iso(ib) = eta_val;
end

T_iso = table( ...
    repmat(a_iso, nBeta, 1), repmat(a_iso, nBeta, 1), repmat(a_iso, nBeta, 1), ...
    eps_iso, eta_iso, beta_vec(:), ...
    'VariableNames', {'a1','a2','a3','eps_bist','eta_val','beta'});

% Append, then remove duplicates by (a1,a2,a3,beta) with a tolerance
anisotropy_study = [anisotropy_study; T_iso];

% Optional: deduplicate (robust to floating noise)
tolMerge = 1e-10;
key = [round(anisotropy_study.a1/tolMerge)*tolMerge, ...
       round(anisotropy_study.a2/tolMerge)*tolMerge, ...
       round(anisotropy_study.a3/tolMerge)*tolMerge, ...
       round(anisotropy_study.beta/tolMerge)*tolMerge];
[~, ia] = unique(key, 'rows', 'stable');
anisotropy_study = anisotropy_study(ia, :);