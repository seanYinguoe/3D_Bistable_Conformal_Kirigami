function [E_total, E_unit, info] = calculate_global_energy(params, opt_beta, opt_t, scale_facs, plot_flag)
%CALCULATE_GLOBAL_ENERGY Global in-plane energy evolution of tessellation.
%   [E_total, E_unit, info] = calculate_global_energy(params, opt_beta, ...
%       opt_t, scale_facs, plot_flag)
%
% This function assigns each unit edge-scale triplet [lambda1, lambda2,
% lambda3], calls deform_triangle_anisotropic to get the unit energy
% evolution (energy vs deployment alpha), then sums all unit curves.
%
% Inputs
%   params     : geometric parameters. Supported formats:
%                - numeric vector [edgeLen; l1; l4; t_default]
%                - struct with fields edgeLen, l1, l4
%   opt_beta   : scalar or N_units-by-1 beta values (per unit)
%   opt_t      : scalar or N_units-by-1 ligament thickness values (per unit)
%   scale_facs : N_units-by-3 edge scale factors per unit
%                [lambda1(lambda12), lambda2(lambda23), lambda3(lambda31)]
%   plot_flag  : legacy input (kept for compatibility; plotting is handled
%                externally in main script)
%
% Outputs
%   E_total    : 1-by-nD global energy curve (sum of all unit curves)
%   E_unit     : N_units-by-nD unit energy curves
%   info       : diagnostics

if nargin < 5 || isempty(plot_flag)
    plot_flag = false;
end
plot_flag = logical(plot_flag); %#ok<NASGU>

validateattributes(scale_facs, {'numeric'}, {'2d', 'ncols', 3, 'nonempty'}, mfilename, 'scale_facs', 4);
N = size(scale_facs, 1);

opt_beta = expand_to_n(opt_beta, N, 'opt_beta');
opt_t = expand_to_n(opt_t, N, 'opt_t');

[edgeLen, l1, l4] = parse_geom_params(params);
nD = 150;   % requested deployment resolution for energy curve
Nseg = 8;   % requested ligament discretisation

E_unit = nan(N, nD);
alpha_unit = nan(N, nD);
valid_unit = false(N, 1);
invalid_reason = repmat({''}, N, 1);

% Try parallel evaluation across units.
use_parallel = false;
pool_size = 0;
if N > 1 && license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
    try
        p = gcp('nocreate');
        if isempty(p)
            try
                p = parpool('threads');
            catch
                p = parpool;
            end
        end
        pool_size = p.NumWorkers;
        use_parallel = pool_size > 1;
    catch
        use_parallel = false;
        pool_size = 0;
    end
end

if use_parallel
    parfor i = 1:N
        [E_row, a_row, ok_i, reason_i] = solve_one_unit(scale_facs(i, :), opt_beta(i), opt_t(i), ...
            edgeLen, l1, l4, nD, Nseg);
        E_unit(i, :) = E_row;
        alpha_unit(i, :) = a_row;
        valid_unit(i) = ok_i;
        invalid_reason{i} = reason_i;
    end
else
    for i = 1:N
        [E_row, a_row, ok_i, reason_i] = solve_one_unit(scale_facs(i, :), opt_beta(i), opt_t(i), ...
            edgeLen, l1, l4, nD, Nseg);
        E_unit(i, :) = E_row;
        alpha_unit(i, :) = a_row;
        valid_unit(i) = ok_i;
        invalid_reason{i} = reason_i;
    end
end

E_total = sum(E_unit, 1, 'omitnan');
alpha_global = mean(alpha_unit, 1, 'omitnan');
if all(~isfinite(alpha_global))
    alpha_global = linspace(0, 1, nD);
end

% Diagnostics
info = struct();
info.max_scale = max(scale_facs(:), [], 'omitnan');
info.min_scale = min(scale_facs(:), [], 'omitnan');
info.N_units = N;
info.N_valid_units = nnz(valid_unit);
info.N_invalid_units = N - info.N_valid_units;
info.invalid_unit_index = find(~valid_unit);
info.invalid_unit_reason = invalid_reason(~valid_unit);
info.alpha = alpha_global;
info.E_unit_final = E_unit(:, end);
info.E_total_final = E_total(end);
info.mean_unit_energy_curve = mean(E_unit, 1, 'omitnan');
info.max_unit_energy_curve = max(E_unit, [], 1, 'omitnan');
info.min_unit_energy_curve = min(E_unit, [], 1, 'omitnan');
info.model = 'deform_triangle_anisotropic per unit';
info.nD = nD;
info.Nseg = Nseg;
info.parallel_used = use_parallel;
info.parallel_workers = pool_size;

end

function vec = expand_to_n(x, N, name)
validateattributes(x, {'numeric'}, {'vector', 'nonempty', 'real', 'finite'}, mfilename, name);
x = x(:);
if isscalar(x)
    vec = repmat(x, N, 1);
elseif numel(x) == N
    vec = x;
else
    error('%s must be scalar or have length N=%d.', name, N);
end
end

function [edgeLen, l1, l4] = parse_geom_params(params)
if isstruct(params)
    assert(isfield(params, 'edgeLen') && isfield(params, 'l1') && isfield(params, 'l4'), ...
        'When params is struct, fields edgeLen, l1, l4 are required.');
    edgeLen = params.edgeLen;
    l1 = params.l1;
    l4 = params.l4;
elseif isnumeric(params) && numel(params) >= 3
    edgeLen = params(1);
    l1 = params(2);
    l4 = params(3);
else
    error('params must be struct or numeric with at least [edgeLen, l1, l4].');
end

validateattributes(edgeLen, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(l1, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(l4, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
end

function [E_row, a_row, ok_i, reason_i] = solve_one_unit(lambda_i, beta_i, t_i, edgeLen, l1, l4, nD, Nseg)
E_row = nan(1, nD);
a_row = nan(1, nD);
ok_i = false;
reason_i = '';

[q1, q2, q3, is_valid, reason_i] = scale_facs_to_q(lambda_i, edgeLen);
if ~is_valid
    return;
end

[Ei, ai] = deform_triangle_anisotropic(q1, q2, q3, edgeLen, l1, l4, beta_i, t_i, nD, Nseg);

Ei = Ei(:).';
ai = ai(:).';
if numel(Ei) ~= nD || numel(ai) ~= nD
    reason_i = 'unexpected_curve_length';
    return;
end

E_row = Ei;
a_row = ai;
ok_i = true;
end
