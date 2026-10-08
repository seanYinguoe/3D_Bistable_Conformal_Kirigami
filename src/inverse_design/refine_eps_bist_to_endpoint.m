function anisotropy_filter_corr = refine_eps_bist_to_endpoint(anisotropy_tbl, opts)
%REFINE_EPS_BIST_TO_ENDPOINT Refine eps_bist so alpha_bist moves toward endpoint.
%   anisotropy_filter_corr = refine_eps_bist_to_endpoint(anisotropy_tbl, opts)
%
% Iteration rule per row:
%   eps_next = eps_curr * alpha_bist
%
% Stopping logic:
%   - converged_alpha_target: finite alpha_bist >= alphaTarget
%   - alpha_nan_endpoint_used: alpha_bist is NaN, use endpoint bistability
%   - invalid_geometry
%   - max_iter_reached
%   - stalled_alpha_change

if nargin < 2
    opts = struct();
end

if ~istable(anisotropy_tbl)
    error('Input must be a table.');
end

req = {'a1','a2','a3','eps_bist','eta_val','beta'};
if ~all(ismember(req, anisotropy_tbl.Properties.VariableNames))
    error('Input table must contain: a1,a2,a3,eps_bist,eta_val,beta');
end

% Defaults
opts = set_default(opts, 'edgeLen', 15);
opts = set_default(opts, 'l1', opts.edgeLen * 0.85);
opts = set_default(opts, 'l4', opts.edgeLen * 0.05);
opts = set_default(opts, 't', opts.edgeLen * 0.015);
opts = set_default(opts, 'etaThreshold', 0.03);
opts = set_default(opts, 'alphaTarget', 0.98);
opts = set_default(opts, 'maxIter', 6);
opts = set_default(opts, 'epsFloor', 1e-4);
opts = set_default(opts, 'epsCeil', 2.0);
opts = set_default(opts, 'alphaTol', 1e-3);

n = height(anisotropy_tbl);

% Keep table workflow and preserve original columns.
anisotropy_filter_corr = anisotropy_tbl;
anisotropy_filter_corr.eps_bist_old = anisotropy_tbl.eps_bist;
anisotropy_filter_corr.eps_bist_new = anisotropy_tbl.eps_bist;
anisotropy_filter_corr.eta_val_old = anisotropy_tbl.eta_val;
anisotropy_filter_corr.eta_val_new = anisotropy_tbl.eta_val;
anisotropy_filter_corr.alpha_bist = nan(n,1);
anisotropy_filter_corr.correction_status = strings(n,1);

% Optional iteration diagnostics (kept table-based and lightweight).
anisotropy_filter_corr.alpha_iter_history = strings(n,1);
anisotropy_filter_corr.eps_iter_history = strings(n,1);

for i = 1:n
    a1 = anisotropy_tbl.a1(i);
    a2 = anisotropy_tbl.a2(i);
    a3 = anisotropy_tbl.a3(i);
    beta = anisotropy_tbl.beta(i);

    eps_old = anisotropy_tbl.eps_bist(i);
    eta_old = anisotropy_tbl.eta_val(i);

    % Default outputs for this row
    eps_curr = eps_old;
    eta_new = eta_old;
    alpha_final = NaN;
    status = "max_iter_reached";

    alpha_hist = nan(opts.maxIter,1);
    eps_hist = nan(opts.maxIter,1);
    k_hist = 0;

    % Skip obviously invalid inputs
    if ~(isfinite(eps_old) && isfinite(eta_old) && eta_old > opts.etaThreshold)
        anisotropy_filter_corr.eps_bist_new(i) = eps_curr;
        anisotropy_filter_corr.eta_val_new(i) = eta_new;
        anisotropy_filter_corr.alpha_bist(i) = alpha_final;
        anisotropy_filter_corr.correction_status(i) = "skip_nonbistable_input";
        anisotropy_filter_corr.alpha_iter_history(i) = "[]";
        anisotropy_filter_corr.eps_iter_history(i) = "[]";
        continue;
    end

    alpha_prev = NaN;

    for it = 1:opts.maxIter
        % 1) Build geometry from current eps_curr
        [okGeom, q1, q2, q3] = row_to_q(a1, a2, a3, eps_curr, opts.edgeLen);
        if ~okGeom
            status = "invalid_geometry";
            break;
        end

        % 2) Evaluate bistability on this geometry
        [~, bistability, info] = bistability_analysis( ...
            opts.l1, opts.l4, beta, opts.t, opts.edgeLen, q1, q2, q3, false);

        % 3) Extract alpha_bist
        alpha_curr = get_alpha_bist(info);

        % Save iter history
        k_hist = k_hist + 1;
        alpha_hist(k_hist) = alpha_curr;
        eps_hist(k_hist) = eps_curr;

        if isfinite(alpha_curr)
            % Case A: detector alpha exists
            alpha_final = alpha_curr;
            eta_new = bistability;

            if alpha_curr >= opts.alphaTarget
                status = "converged_alpha_target";
                break;
            end

            if isfinite(alpha_prev) && abs(alpha_curr - alpha_prev) < opts.alphaTol
                status = "stalled_alpha_change";
                break;
            end

            % 4) Update rule
            eps_next = eps_curr * alpha_curr;
            eps_next = min(max(eps_next, opts.epsFloor), opts.epsCeil);

            if ~isfinite(eps_next)
                status = "invalid_eps_update";
                break;
            end

            eps_curr = eps_next;
            alpha_prev = alpha_curr;
        else
            % Case B: detector alpha disappeared; use endpoint-based eta
            eta_endpoint = get_endpoint_bistability(info);
            if isfinite(eta_endpoint)
                eta_new = eta_endpoint;
            else
                eta_new = NaN;
            end
            alpha_final = NaN;
            status = "alpha_nan_endpoint_used";
            break;
        end

        if it == opts.maxIter
            status = "max_iter_reached";
        end
    end

    % Write row outputs
    anisotropy_filter_corr.eps_bist_new(i) = eps_curr;
    anisotropy_filter_corr.eps_bist(i) = eps_curr;
    anisotropy_filter_corr.eta_val_new(i) = eta_new;
    anisotropy_filter_corr.eta_val(i) = eta_new;
    anisotropy_filter_corr.alpha_bist(i) = alpha_final;
    anisotropy_filter_corr.correction_status(i) = status;

    if k_hist > 0
        anisotropy_filter_corr.alpha_iter_history(i) = mat2str(alpha_hist(1:k_hist).', 4);
        anisotropy_filter_corr.eps_iter_history(i) = mat2str(eps_hist(1:k_hist).', 6);
    else
        anisotropy_filter_corr.alpha_iter_history(i) = "[]";
        anisotropy_filter_corr.eps_iter_history(i) = "[]";
    end
end

end

function opts = set_default(opts, name, value)
if ~isfield(opts, name) || isempty(opts.(name))
    opts.(name) = value;
end
end

function [ok, q1, q2, q3] = row_to_q(a1, a2, a3, eps_bist, edgeLen)
% Geometry construction exactly as requested.
q1 = [0, 0, 0];
q2 = [0, 0, 0];
q3 = [0, 0, 0];
ok = false;

if ~(isfinite(a1) && isfinite(a2) && isfinite(a3) && isfinite(eps_bist) && isfinite(edgeLen))
    return;
end
if abs(sin(a3)) < 1e-12
    return;
end

lam3 = 1 + eps_bist;
lam1 = (sin(a1)/sin(a3)) * lam3;
lam2 = (sin(a2)/sin(a3)) * lam3;

if ~(isfinite(lam1) && isfinite(lam2) && isfinite(lam3) && lam1 > 0 && lam2 > 0 && lam3 > 0)
    return;
end

L1 = lam1 * edgeLen;
L2 = lam2 * edgeLen;
L3 = lam3 * edgeLen;

if ~(isfinite(L1) && isfinite(L2) && isfinite(L3) && L1 > 0 && L2 > 0 && L3 > 0)
    return;
end

% Triangle inequality check
if ~((L1 + L2 > L3) && (L1 + L3 > L2) && (L2 + L3 > L1))
    return;
end

q1 = [0, 0, 0];
q2 = [0, -L3, 0];
q3_y = (L1^2 - L2^2 - L3^2) / (2 * L3);
arg = L2^2 - q3_y^2;
if ~(isfinite(arg) && arg >= 0)
    return;
end
q3_x = -sqrt(max(arg, 0));
q3 = [q3_x, q3_y, 0];

ok = all(isfinite([q1, q2, q3]));
end

function alpha_bist = get_alpha_bist(info)
% Prefer direct alpha_bist, fallback to min_idx/alpha_history.
alpha_bist = NaN;
if ~isstruct(info)
    return;
end

if isfield(info, 'alpha_bist') && isfinite(info.alpha_bist)
    alpha_bist = info.alpha_bist;
    return;
end

if isfield(info, 'min_idx') && isfield(info, 'alpha_history')
    idx = info.min_idx;
    ah = info.alpha_history;
    if isfinite(idx) && idx >= 1 && idx <= numel(ah)
        val = ah(idx);
        if isfinite(val)
            alpha_bist = val;
        end
    end
end
end

function eta_endpoint = get_endpoint_bistability(info)
% Endpoint fallback: eta = (E_peak - E_end) / E_peak
eta_endpoint = NaN;
if ~isstruct(info) || ~isfield(info, 'E_history')
    return;
end

E = info.E_history(:);
E = E(isfinite(E));
if isempty(E)
    return;
end

E_peak = max(E);
E_end = E(end);
if ~(isfinite(E_peak) && isfinite(E_end))
    return;
end
if abs(E_peak) < eps
    return;
end

eta_endpoint = (E_peak - E_end) / E_peak;
if ~isfinite(eta_endpoint)
    eta_endpoint = NaN;
end
end
