function anisotropy_study_corr = refine_eps_bist_to_endpoint(anisotropy_study_pi12, opts)
%REFINE_EPS_BIST_TO_ENDPOINT Post-process eps_bist to make alpha_bist -> endpoint.
%   anisotropy_study_corr = refine_eps_bist_to_endpoint(anisotropy_study_pi12, opts)
%
% This does NOT rerun full ternary study. It only refines rows that are
% already considered bistable according to existing table fields.
%
% Required table variables:
%   a1, a2, a3, eps_bist, eta_val, beta
%
% opts fields (with defaults):
%   edgeLen       : 15
%   l1            : edgeLen*0.85
%   l4            : edgeLen*0.05
%   t             : edgeLen*0.015
%   etaThreshold  : 0.03
%   alphaTarget   : 0.95
%   maxIter       : 6
%   relax         : 0.8   (0~1, damping for update)
%   epsFloor      : 1e-4
%   epsCeil       : 2.0
%   alphaTol      : 1e-3
%   alphaCloseTol : 0.02
%
% Update rule requested:
%   eps_new = eps_old * alpha_old
% Iterate until alpha is close to alphaTarget.

if nargin < 2
    opts = struct();
end

if ~istable(anisotropy_study_pi12)
    error('Input must be a table.');
end

req = {'a1','a2','a3','eps_bist','eta_val','beta'};
if ~all(ismember(req, anisotropy_study_pi12.Properties.VariableNames))
    error('Input table must contain: a1,a2,a3,eps_bist,eta_val,beta');
end

% Defaults
opts = set_default(opts, 'edgeLen', 15);
opts = set_default(opts, 'l1', opts.edgeLen * 0.85);
opts = set_default(opts, 'l4', opts.edgeLen * 0.05);
opts = set_default(opts, 't', opts.edgeLen * 0.015);
opts = set_default(opts, 'etaThreshold', 0.03);
opts = set_default(opts, 'alphaTarget', 0.95);
opts = set_default(opts, 'maxIter', 6);
opts = set_default(opts, 'relax', 0.8);
opts = set_default(opts, 'epsFloor', 1e-4);
opts = set_default(opts, 'epsCeil', 2.0);
opts = set_default(opts, 'alphaTol', 1e-3);
opts = set_default(opts, 'alphaCloseTol', 0.02);

n = height(anisotropy_study_pi12);
anisotropy_study_corr = anisotropy_study_pi12;

% Diagnostics
anisotropy_study_corr.eps_bist_old = anisotropy_study_pi12.eps_bist;
anisotropy_study_corr.eps_bist_new = anisotropy_study_pi12.eps_bist;
anisotropy_study_corr.alpha_old = nan(n,1);
anisotropy_study_corr.alpha_new = nan(n,1);
anisotropy_study_corr.eta_old = anisotropy_study_pi12.eta_val;
anisotropy_study_corr.eta_new = anisotropy_study_pi12.eta_val;
anisotropy_study_corr.n_iter_correction = zeros(n,1);
anisotropy_study_corr.correction_status = strings(n,1);
anisotropy_study_corr.corrected_flag = false(n,1);

for i = 1:n
    eps0 = anisotropy_study_pi12.eps_bist(i);
    eta0 = anisotropy_study_pi12.eta_val(i);

    % Only process rows already considered bistable by existing fields.
    if ~(isfinite(eps0) && isfinite(eta0) && eta0 > opts.etaThreshold)
        anisotropy_study_corr.correction_status(i) = "skip_nonbistable_input";
        continue;
    end

    a1 = anisotropy_study_pi12.a1(i);
    a2 = anisotropy_study_pi12.a2(i);
    a3 = anisotropy_study_pi12.a3(i);
    beta = anisotropy_study_pi12.beta(i);

    % Evaluate current row
    [okGeom, q1, q2, q3] = row_to_q(a1, a2, a3, eps0, opts.edgeLen);
    if ~okGeom
        anisotropy_study_corr.correction_status(i) = "invalid_geometry_old";
        continue;
    end

    [eps_eval, eta_eval, info] = bistability_analysis( ...
        opts.l1, opts.l4, beta, opts.t, opts.edgeLen, q1, q2, q3, false);

    alpha_curr = get_alpha_bist(info);
    anisotropy_study_corr.alpha_old(i) = alpha_curr;
    anisotropy_study_corr.eta_old(i) = eta_eval;

    if ~(isfinite(alpha_curr) && isfinite(eta_eval) && strcmp(info.status, 'bistable'))
        anisotropy_study_corr.correction_status(i) = "skip_not_bistable_rerun";
        continue;
    end

    if abs(alpha_curr - opts.alphaTarget) <= opts.alphaCloseTol
        anisotropy_study_corr.correction_status(i) = "accept_no_change";
        anisotropy_study_corr.alpha_new(i) = alpha_curr;
        anisotropy_study_corr.eta_new(i) = eta_eval;
        continue;
    end

    % Iterative refinement
    eps_curr = eps0;
    eta_curr = eta_eval;
    alpha_prev = alpha_curr;
    status = "max_iter_reached";
    n_iter = 0;

    for it = 1:opts.maxIter
        n_iter = it;
        % Requested update rule:
        % eps_new = eps_old * alpha_old
        eps_next = eps_curr * alpha_prev;
        eps_next = min(max(eps_next, opts.epsFloor), opts.epsCeil);

        [okGeom, q1, q2, q3] = row_to_q(a1, a2, a3, eps_next, opts.edgeLen);
        if ~okGeom
            status = "invalid_geometry_iter";
            break;
        end

        [eps_eval_next, eta_eval_next, info_next] = bistability_analysis( ...
            opts.l1, opts.l4, beta, opts.t, opts.edgeLen, q1, q2, q3, false);

        alpha_next = get_alpha_bist(info_next);
        if ~(isfinite(alpha_next) && isfinite(eta_eval_next) && strcmp(info_next.status, 'bistable'))
            status = "lost_bistability_iter";
            break;
        end

        eps_curr = eps_next;
        eta_curr = eta_eval_next;
        alpha_prev = alpha_next;

        if abs(alpha_prev - opts.alphaTarget) <= opts.alphaCloseTol
            status = "converged_alpha_target";
            break;
        end

        if abs(alpha_prev - alpha_curr) < opts.alphaTol
            status = "stalled_alpha_change";
            break;
        end
        alpha_curr = alpha_prev;
    end

    anisotropy_study_corr.eps_bist_new(i) = eps_curr;
    anisotropy_study_corr.eps_bist(i) = eps_curr;
    anisotropy_study_corr.alpha_new(i) = alpha_prev;
    anisotropy_study_corr.eta_new(i) = eta_curr;
    anisotropy_study_corr.eta_val(i) = eta_curr;
    anisotropy_study_corr.n_iter_correction(i) = n_iter;
    anisotropy_study_corr.correction_status(i) = status;
    anisotropy_study_corr.corrected_flag(i) = abs(eps_curr - eps0) > 1e-12;
end

end

function opts = set_default(opts, name, value)
if ~isfield(opts, name) || isempty(opts.(name))
    opts.(name) = value;
end
end

function [ok, q1, q2, q3] = row_to_q(a1, a2, a3, eps_bist, edgeLen)
q1 = [0,0,0];
q2 = [0,0,0];
q3 = [0,0,0];

if ~(isfinite(a1) && isfinite(a2) && isfinite(a3) && isfinite(eps_bist) && isfinite(edgeLen))
    ok = false;
    return;
end
if abs(sin(a3)) < 1e-12
    ok = false;
    return;
end

lam3 = 1 + eps_bist;
lam1 = (sin(a1)/sin(a3)) * lam3;
lam2 = (sin(a2)/sin(a3)) * lam3;

L1 = lam1 * edgeLen;
L2 = lam2 * edgeLen;
L3 = lam3 * edgeLen;

if ~(isfinite(L1) && isfinite(L2) && isfinite(L3) && L1 > 0 && L2 > 0 && L3 > 0)
    ok = false;
    return;
end
if ~((L1+L2>L3) && (L1+L3>L2) && (L2+L3>L1))
    ok = false;
    return;
end

q1 = [0, 0, 0];
q2 = [0, -L3, 0];
q3_y = (L1^2 - L2^2 - L3^2) / (2*L3);
inside = L2^2 - q3_y^2;
if inside <= 0
    ok = false;
    return;
end
q3_x = -sqrt(inside);
q3 = [q3_x, q3_y, 0];
ok = all(isfinite([q1,q2,q3]));
end

function alpha_bist = get_alpha_bist(info)
alpha_bist = NaN;
if isfield(info, 'alpha_bist') && isfinite(info.alpha_bist)
    alpha_bist = info.alpha_bist;
    return;
end
if isfield(info, 'min_idx') && isfield(info, 'alpha_history') ...
        && isfinite(info.min_idx) && info.min_idx >= 1 ...
        && info.min_idx <= numel(info.alpha_history)
    alpha_bist = info.alpha_history(info.min_idx);
end
end
