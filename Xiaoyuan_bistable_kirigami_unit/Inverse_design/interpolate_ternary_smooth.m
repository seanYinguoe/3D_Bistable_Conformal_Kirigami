function T_plot = interpolate_ternary_smooth(T_filter, nGrid, smoothFactor, etaThreshold)
%INTERPOLATE_TERNARY_SMOOTH Dense ternary interpolation for eps_bist and eta_val.
%   T_plot = interpolate_ternary_smooth(T_filter, nGrid, smoothFactor, etaThreshold)

if nargin < 2 || isempty(nGrid)
    nGrid = 60;
end
if nargin < 3 || isempty(smoothFactor)
    smoothFactor = 1.3;
end
if nargin < 4 || isempty(etaThreshold)
    etaThreshold = 0.15;
end

requiredVars = {'a1','a2','a3','eps_bist','eta_val'};
if ~all(ismember(requiredVars, T_filter.Properties.VariableNames))
    error('T_filter must contain: a1, a2, a3, eps_bist, eta_val');
end

if isempty(T_filter)
    T_plot = table([], [], [], [], [], ...
        'VariableNames', {'a1','a2','a3','eps_bist','eta_val'});
    return;
end

% Bistable support points only
isB = isfinite(T_filter.eps_bist) & isfinite(T_filter.eta_val) & (T_filter.eta_val > etaThreshold);
TB = T_filter(isB, :);
if isempty(TB)
    T_plot = table([], [], [], [], [], ...
        'VariableNames', {'a1','a2','a3','eps_bist','eta_val'});
    return;
end

% STEP 1 — ternary (a1,a2,a3) -> 2D Cartesian
sB = TB.a1 + TB.a2 + TB.a3;
BB = TB.a2 ./ sB;
CB = TB.a3 ./ sB;
xB = BB + 0.5 * CB;
yB = (sqrt(3)/2) * CB;

% STEP 3 — dense barycentric grid
nAlloc = (nGrid + 1) * (nGrid + 2) / 2;
Aq = zeros(nAlloc,1);
Bq = zeros(nAlloc,1);
Cq = zeros(nAlloc,1);
id = 0;
for i = 0:nGrid
    for j = 0:(nGrid - i)
        k = nGrid - i - j;
        id = id + 1;
        Aq(id) = i / nGrid;
        Bq(id) = j / nGrid;
        Cq(id) = k / nGrid;
    end
end
xq = Bq + 0.5 * Cq;
yq = (sqrt(3)/2) * Cq;

% Pairwise distances query<->bistable points
Gq = xq.^2 + yq.^2;
Gb = xB.^2 + yB.^2;
D2_qb = max(Gq + Gb' - 2*(xq*xB' + yq*yB'), 0);
D_qb = sqrt(D2_qb);

% Reference kNN scale from original bistable points
nB = numel(xB);
k = min(6, nB);
if k < 1
    T_plot = table([], [], [], [], [], ...
        'VariableNames', {'a1','a2','a3','eps_bist','eta_val'});
    return;
end
if nB > 1
    D2_bb = max(Gb + Gb' - 2*(xB*xB' + yB*yB'), 0);
    D_bb = sqrt(D2_bb);
    D_bb(1:nB+1:end) = inf;
    Dsort_bb = sort(D_bb, 2, 'ascend');
    k_ref = min(k, nB - 1);
    meanNN_b = mean(Dsort_bb(:,1:k_ref), 2);
else
    meanNN_b = 0;
end

d_med = median(meanNN_b);
d_mad = median(abs(meanNN_b - d_med));
if d_mad <= eps
    supportThr = inf;
else
    supportThr = d_med + smoothFactor * 1.4826 * d_mad;
end

% Supported region mask from query-point kNN to bistable samples
Dsort_qb = sort(D_qb, 2, 'ascend');
kq = min(k, size(Dsort_qb,2));
meanNN_q = mean(Dsort_qb(:,1:kq), 2);
isSupported = meanNN_q <= supportThr;

% Also require query point inside convex hull of bistable points
if nB >= 3
    hullIdx = convhull(xB, yB);
    isSupported = isSupported & inpolygon(xq, yq, xB(hullIdx), yB(hullIdx));
end

% STEP 4 — IDW interpolation (p=2) for eps and eta
tiny = 1e-12;
p = 2;
W = 1 ./ (D_qb.^p + tiny);
W(~isSupported, :) = 0;

epsVals = TB.eps_bist(:);
etaVals = TB.eta_val(:);
wSum = sum(W, 2);

eps_q = NaN(size(xq));
eta_q = NaN(size(xq));
ok = isSupported & (wSum > 0);
eps_q(ok) = (W(ok,:) * epsVals) ./ wSum(ok);
eta_q(ok) = (W(ok,:) * etaVals) ./ wSum(ok);

% STEP 5 — mild smoothing within supported region only (no NaN filling)
% second-pass IDW over the SAME sample set (TB), using larger k_s
if smoothFactor > 0 && any(ok) && nB > 1
    [~, idxSort] = sort(D_qb, 2, 'ascend');

    k_s = max(1, min(nB, max(8, round(8 * smoothFactor))));
    idxK = idxSort(:, 1:k_s);

    eps_s = eps_q;
    eta_s = eta_q;
    ii_list = find(ok);
    for jj = 1:numel(ii_list)
        ii = ii_list(jj);

        % Ensure valid column-vector sample indices
        nb = idxK(ii, :).';
        nb = nb(isfinite(nb) & nb >= 1 & nb <= nB);
        if isempty(nb)
            continue;
        end
        nb = unique(nb, 'stable');

        d_nb = D_qb(ii, nb).';
        good = isfinite(d_nb) & ~isnan(d_nb);
        nb = nb(good);
        d_nb = d_nb(good);
        if isempty(nb)
            continue;
        end

        % Column vectors with matched sizes
        w_col = 1 ./ (d_nb.^p + tiny);
        eps_col = epsVals(nb);
        eta_col = etaVals(nb);

        ws = sum(w_col);
        if ws > 0
            eps_s(ii) = sum(w_col .* eps_col) / ws;
            eta_s(ii) = sum(w_col .* eta_col) / ws;
        end
    end

    eps_q(ok) = eps_s(ok);
    eta_q(ok) = eta_s(ok);
elseif smoothFactor > 0 && any(ok) && nB == 1
    % Graceful single-sample case: keep first-pass constant IDW values
end

% Keep only supported finite points
keep = isfinite(eps_q) & isfinite(eta_q);
T_plot = table(Aq(keep)*pi, Bq(keep)*pi, Cq(keep)*pi, eps_q(keep), eta_q(keep), ...
    'VariableNames', {'a1','a2','a3','eps_bist','eta_val'});
end
