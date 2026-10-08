function [opt_beta, bistability, info] = assign_opt_beta(anisotropy_level, anisotropy_study, opts)
%ASSIGN_OPT_BETA Assign an optimal tilting angle beta for each unit.
%   [opt_beta, bistability, info] = assign_opt_beta(anisotropy_level, anisotropy_study, opts)
%
% Inputs
%   anisotropy_level : table/struct with fields a1, a2, a3, strain3
%                      or numeric Nx4 array [a1 a2 a3 strain3]
%   anisotropy_study : table/struct with fields a1, a2, a3, eps_bist, eta_val, beta
%                      or numeric Nx6 array [a1 a2 a3 eps_bist eta_val beta]
%   opts             : optional settings struct
%
% Outputs
%   opt_beta     : Nx1 chosen beta per unit, always clamped to [0, pi/12]
%   bistability  : Nx1 predicted eta_val at chosen beta for bistable units,
%                  NaN for monostable units
%   info         : diagnostic struct with local interpolation details
%
% The assignment uses a local interpolation strategy:
%   1) find nearest samples in (a1,a2,a3),
%   2) aggregate nearby samples by beta with distance weights,
%   3) interpolate eps_bist(beta) and eta_val(beta),
%   4) choose beta to match strain3 by eps_bist as closely as possible.

if nargin < 3
    opts = struct();
end
opts = fill_opts(opts);

[a1_level, a2_level, a3_level, strain3] = unpack_level(anisotropy_level);
[a1, a2, a3, eps_bist, eta_val, beta] = unpack_study(anisotropy_study);

nUnit = numel(strain3);
beta_max = pi / 12;
beta_grid = linspace(0, beta_max, opts.nBetaEval).';

% Rim mask: rim triangles are locked at beta_max (undeployed boundary).
is_rim = false(nUnit, 1);
if ~isempty(opts.rimMask)
    rm = logical(opts.rimMask(:));
    if numel(rm) == nUnit
        is_rim = rm;
    else
        warning('assign_opt_beta:rimMaskSize', ...
            'rimMask length (%d) does not match nUnit (%d); ignoring.', numel(rm), nUnit);
    end
end

opt_beta = zeros(nUnit, 1);
bistability = NaN(nUnit, 1);

info = struct();
info.flag = repmat("unprocessed", nUnit, 1);
info.eps_pred = NaN(nUnit, 1);
info.eta_pred = NaN(nUnit, 1);
info.eps_bist_max = NaN(nUnit, 1);
info.beta_error = NaN(nUnit, 1);
info.n_local = zeros(nUnit, 1);
info.n_feasible = zeros(nUnit, 1);
info.beta_grid = beta_grid;
info.opts = opts;
info.n_rim = sum(is_rim);

study_valid = isfinite(a1) & isfinite(a2) & isfinite(a3) & isfinite(beta) ...
    & (beta >= 0) & (beta <= beta_max);

a1 = a1(study_valid);
a2 = a2(study_valid);
a3 = a3(study_valid);
eps_bist = eps_bist(study_valid);
eta_val = eta_val(study_valid);
beta = beta(study_valid);

if isempty(beta)
    info.flag(:) = "no_valid_study_samples";
    return;
end

study_angles = [a1, a2, a3];

for i = 1:nUnit
    % Rim triangles sit on the free boundary and stay undeployed; assign
    % beta_max so they contribute no snap but don't break the pattern.
    if is_rim(i)
        opt_beta(i)    = beta_max;
        bistability(i) = NaN;
        info.flag(i)   = "rim_undeployed";
        continue;
    end

    if ~isfinite(a1_level(i)) || ~isfinite(a2_level(i)) || ~isfinite(a3_level(i)) || ~isfinite(strain3(i))
        opt_beta(i) = 0;
        bistability(i) = NaN;
        info.flag(i) = "invalid_unit_input";
        continue;
    end

    query = [a1_level(i), a2_level(i), a3_level(i)];
    d2 = sum((study_angles - query).^2, 2);
    [d2_sorted, idx_sorted] = sort(d2, 'ascend');

    k_eff = min(opts.kNN, numel(idx_sorted));
    idx_local = idx_sorted(1:k_eff);
    d_local = sqrt(d2_sorted(1:k_eff));

    info.n_local(i) = numel(idx_local);
    if numel(idx_local) < opts.minSamples
        opt_beta(i) = 0;
        bistability(i) = NaN;
        info.flag(i) = "too_few_samples";
        continue;
    end

    local = local_interp_beta( ...
        beta(idx_local), eps_bist(idx_local), eta_val(idx_local), d_local, beta_grid, opts);

    if local.n_support < opts.minSamples || all(~isfinite(local.eps_hat))
        opt_beta(i) = 0;
        bistability(i) = NaN;
        info.flag(i) = "too_few_samples";
        continue;
    end

    % Selection is driven by eps matching only (ignore eta in optimization).
    feasible = isfinite(local.eps_hat);
    info.n_feasible(i) = nnz(feasible);

    if any(feasible)
        info.eps_bist_max(i) = max(local.eps_hat(feasible));
    end

    if ~any(feasible)
        opt_beta(i) = 0;
        bistability(i) = NaN;
        info.flag(i) = "no_eps_interp";
        continue;
    end

    eps_feas = local.eps_hat(feasible);
    beta_feas = beta_grid(feasible);
    err_feas = abs(eps_feas - strain3(i));
    [~, idx_pick] = min(err_feas);

    opt_beta(i) = clamp_beta(beta_feas(idx_pick), beta_max);
    % keep bistability output for compatibility if eta interpolation exists
    eta_feas = local.eta_hat(feasible);
    if ~isempty(eta_feas) && all(isfinite(eta_feas))
        bistability(i) = eta_feas(idx_pick);
    else
        bistability(i) = NaN;
    end
    info.eps_pred(i) = eps_feas(idx_pick);
    if exist('eta_feas', 'var') && numel(eta_feas) >= idx_pick
        info.eta_pred(i) = eta_feas(idx_pick);
    else
        info.eta_pred(i) = NaN;
    end
    info.beta_error(i) = err_feas(idx_pick);
    info.flag(i) = "eps_matched";
end
end

function opts = fill_opts(opts)
if nargin < 1 || isempty(opts)
    opts = struct();
end

opts = set_default(opts, 'kNN', 80);
opts = set_default(opts, 'nBetaEval', 61);
opts = set_default(opts, 'etaMin', 0.10);
opts = set_default(opts, 'strainTol', 0.015);
opts = set_default(opts, 'w_eps', 1.0);
opts = set_default(opts, 'w_eta', 0.4);
opts = set_default(opts, 'minSamples', 6);
opts = set_default(opts, 'distPower', 2.0);
opts = set_default(opts, 'betaMergeTol', 1e-10);
opts = set_default(opts, 'epsClipMin', 0.0);
opts = set_default(opts, 'etaClipMin', 0.0);
opts = set_default(opts, 'rimMask', []);

opts.kNN = max(1, round(opts.kNN));
opts.nBetaEval = max(5, round(opts.nBetaEval));
opts.minSamples = max(3, round(opts.minSamples));
opts.etaMin = max(0, opts.etaMin);
opts.strainTol = max(0, opts.strainTol);
opts.w_eps = max(0, opts.w_eps);
opts.w_eta = max(0, opts.w_eta);
opts.distPower = max(0.5, opts.distPower);
end

function s = set_default(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end

function [a1, a2, a3, strain3] = unpack_level(anisotropy_level)
if isnumeric(anisotropy_level)
    if size(anisotropy_level, 2) < 4
        error('anisotropy_level numeric input must be Nx4+ with columns [a1 a2 a3 strain3].');
    end
    a1 = anisotropy_level(:, 1);
    a2 = anisotropy_level(:, 2);
    a3 = anisotropy_level(:, 3);
    strain3 = anisotropy_level(:, 4);
else
    req = {'a1', 'a2', 'a3', 'strain3'};
    check_fields(anisotropy_level, req, 'anisotropy_level');
    a1 = anisotropy_level.a1(:);
    a2 = anisotropy_level.a2(:);
    a3 = anisotropy_level.a3(:);
    strain3 = anisotropy_level.strain3(:);
end
if ~(numel(a1) == numel(a2) && numel(a2) == numel(a3) && numel(a3) == numel(strain3))
    error('anisotropy_level fields a1, a2, a3, strain3 must have the same length.');
end
end

function [a1, a2, a3, eps_bist, eta_val, beta] = unpack_study(anisotropy_study)
if isnumeric(anisotropy_study)
    if size(anisotropy_study, 2) < 6
        error('anisotropy_study numeric input must be Nx6+ with columns [a1 a2 a3 eps_bist eta_val beta].');
    end
    a1 = anisotropy_study(:, 1);
    a2 = anisotropy_study(:, 2);
    a3 = anisotropy_study(:, 3);
    eps_bist = anisotropy_study(:, 4);
    eta_val = anisotropy_study(:, 5);
    beta = anisotropy_study(:, 6);
else
    req = {'a1', 'a2', 'a3', 'eps_bist', 'eta_val', 'beta'};
    check_fields(anisotropy_study, req, 'anisotropy_study');
    a1 = anisotropy_study.a1(:);
    a2 = anisotropy_study.a2(:);
    a3 = anisotropy_study.a3(:);
    eps_bist = anisotropy_study.eps_bist(:);
    eta_val = anisotropy_study.eta_val(:);
    beta = anisotropy_study.beta(:);
end
if ~(numel(a1) == numel(a2) && numel(a2) == numel(a3) && numel(a3) == numel(eps_bist) ...
        && numel(eps_bist) == numel(eta_val) && numel(eta_val) == numel(beta))
    error('anisotropy_study fields a1, a2, a3, eps_bist, eta_val, beta must have the same length.');
end
end

function check_fields(S, req, label)
if ~all(isfieldish(S, req))
    error('%s must contain fields: %s', label, strjoin(req, ', '));
end
end

function tf = isfieldish(S, req)
tf = false(size(req));
if istable(S)
    vars = S.Properties.VariableNames;
    for k = 1:numel(req)
        tf(k) = ismember(req{k}, vars);
    end
elseif isstruct(S)
    for k = 1:numel(req)
        tf(k) = isfield(S, req{k});
    end
else
    error('Inputs must be table, struct, or numeric arrays with expected columns.');
end
end

function local = local_interp_beta(beta, eps_bist, eta_val, d_local, beta_grid, opts)
beta_max = pi / 20;
valid_beta = isfinite(beta) & (beta >= 0) & (beta <= beta_max);

beta = beta(valid_beta);
eps_bist = eps_bist(valid_beta);
eta_val = eta_val(valid_beta);
d_local = d_local(valid_beta);

local = struct();
local.eps_hat = NaN(size(beta_grid));
local.eta_hat = NaN(size(beta_grid));
local.n_support = numel(beta);

if numel(beta) < opts.minSamples
    return;
end

weights = 1 ./ max(d_local, 1e-12) .^ opts.distPower;

eps_mask = isfinite(eps_bist);
eta_mask = isfinite(eta_val);

eps_data = aggregate_by_beta(beta(eps_mask), eps_bist(eps_mask), weights(eps_mask), opts);
eta_data = aggregate_by_beta(beta(eta_mask), eta_val(eta_mask), weights(eta_mask), opts);

if numel(eps_data.beta_u) >= 2
    local.eps_hat = interp_on_grid(eps_data.beta_u, eps_data.val_u, beta_grid);
    local.eps_hat = max(local.eps_hat, opts.epsClipMin);
end

if numel(eta_data.beta_u) >= 2
    local.eta_hat = interp_on_grid(eta_data.beta_u, eta_data.val_u, beta_grid);
    local.eta_hat = max(local.eta_hat, opts.etaClipMin);
end

local.n_support = sum(eps_mask);
end

function out = aggregate_by_beta(beta, values, weights, opts)
out = struct('beta_u', [], 'val_u', []);
if isempty(beta)
    return;
end

[beta_sorted, ord] = sort(beta(:), 'ascend');
val_sorted = values(ord);
w_sorted = weights(ord);

beta_group = [];
val_group = [];

i = 1;
n = numel(beta_sorted);
while i <= n
    j = i;
    while j < n && abs(beta_sorted(j + 1) - beta_sorted(i)) <= opts.betaMergeTol
        j = j + 1;
    end

    w = w_sorted(i:j);
    v = val_sorted(i:j);
    wsum = sum(w);
    if wsum > 0
        beta_group(end+1, 1) = mean(beta_sorted(i:j)); %#ok<AGROW>
        val_group(end+1, 1) = sum(w .* v) / wsum; %#ok<AGROW>
    end
    i = j + 1;
end

out.beta_u = beta_group;
out.val_u = val_group;
end

function yq = interp_on_grid(x, y, xq)
yq = NaN(size(xq));
if numel(x) < 2
    return;
end

[x, ord] = sort(x(:), 'ascend');
y = y(ord);

if any(diff(x) <= 0)
    [x, ia] = unique(x, 'stable');
    y = y(ia);
end

if numel(x) < 2
    return;
end

inside = (xq >= x(1)) & (xq <= x(end));
if any(inside)
    yq(inside) = interp1(x, y, xq(inside), 'linear');
end
end

function idx = choose_best_tol(err_feas, eta_feas, within_tol)
idx_cand = find(within_tol);
[~, ord] = sortrows([-eta_feas(idx_cand), err_feas(idx_cand), idx_cand]);
idx = idx_cand(ord(1));
end

function idx = choose_best_objective(J, err_feas, eta_feas)
idx_all = (1:numel(J)).';
[~, ord] = sortrows([J, err_feas, -eta_feas, idx_all]);
idx = idx_all(ord(1));
end

function beta = clamp_beta(beta, beta_max)
beta = min(max(beta, 0), beta_max);
end
