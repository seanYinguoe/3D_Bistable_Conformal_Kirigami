function [tblClean, reportClean] = clean_isolated_outliers(tbl)
%CLEAN_ISOLATED_OUTLIERS Remove isolated/high-anisotropy support islands per beta.
%   [tblClean, reportClean] = clean_isolated_outliers(tbl)
%   Input columns: [a1 a2 a3 eps_bist eta_val beta] (table or numeric)
%
% Logic per beta slice:
%   1) Build neighbor graph in (a1,a2,a3) with
%        neighTol = 1.25 * median(min nonzero pairwise distance per node)
%   2) On bistable support nodes (finite eps & eta):
%      - Remove nodes outside largest connected component.
%      - Remove nodes that are both isolated (< minNeighbors) and high-anisotropy.
%   3) Removal means forcing eps_bist and eta_val to NaN (rows kept).

isTable = istable(tbl);

if isTable
    req = {'a1','a2','a3','eps_bist','eta_val','beta'};
    if ~all(ismember(req, tbl.Properties.VariableNames))
        error('Input table must contain: a1,a2,a3,eps_bist,eta_val,beta');
    end
    A = [tbl.a1, tbl.a2, tbl.a3, tbl.eps_bist, tbl.eta_val, tbl.beta];
else
    if ~isnumeric(tbl) || size(tbl,2) < 6
        error('Numeric input must be Nx6+ [a1 a2 a3 eps_bist eta_val beta].');
    end
    A = tbl(:,1:6);
end

n = size(A,1);
if n == 0
    tblClean = tbl;
    reportClean = struct('removed_isolated_high_anis_idx', [], ...
                         'removed_small_components_idx', [], ...
                         'neighTol_per_beta', []);
    return;
end

% Column map
cA1 = 1; cA2 = 2; cA3 = 3; cEPS = 4; cETA = 5; cBETA = 6;

% Settings
tol = 1e-10;
minNeighbors = 2;
anisTolFactor = 1.15;
neighScale = 1.25;

anis = max(A(:,[cA1 cA2 cA3]), [], 2) - min(A(:,[cA1 cA2 cA3]), [], 2);

removedIsoHigh = [];
removedSmallCC = [];

betaVals = unique(A(:,cBETA));
neighTol_per_beta = nan(numel(betaVals), 2); % [beta, neighTol]

for ib = 1:numel(betaVals)
    beta0 = betaVals(ib);
    inBeta = abs(A(:,cBETA) - beta0) <= tol;
    idxBeta = find(inBeta);
    nB = numel(idxBeta);
    if nB < 2
        continue;
    end

    Xb = A(idxBeta, [cA1 cA2 cA3]);

    % Pairwise distances for this beta slice
    G = sum(Xb.^2, 2);
    D2 = max(G + G' - 2*(Xb*Xb'), 0);
    D = sqrt(D2);
    D(1:nB+1:end) = inf;

    % Smallest nonzero distance per node (ignoring inf)
    dMin = min(D, [], 2);
    dMin = dMin(isfinite(dMin) & dMin > 0);
    if isempty(dMin)
        continue;
    end

    neighTol = neighScale * median(dMin);
    neighTol_per_beta(ib,:) = [beta0, neighTol];

    % Support nodes in this beta
    isSupportGlobal = inBeta & isfinite(A(:,cEPS)) & isfinite(A(:,cETA));
    idxSupportGlobal = find(isSupportGlobal);
    if isempty(idxSupportGlobal)
        continue;
    end

    % Map support rows into beta-local indexing
    [~, locSupport] = ismember(idxSupportGlobal, idxBeta);
    locSupport = locSupport(locSupport > 0);
    if numel(locSupport) < 1
        continue;
    end

    % Build adjacency among support nodes only
    DS = D(locSupport, locSupport);
    AS = (DS <= neighTol);
    nS = size(AS,1);
    AS(1:nS+1:end) = false;

    % Keep largest connected component among support nodes
    compId = local_conn_components(AS);
    if ~isempty(compId)
        counts = accumarray(compId, 1);
        [~, iLargest] = max(counts);
        inLargest = (compId == iLargest);

        if any(~inLargest)
            idxDrop = idxSupportGlobal(~inLargest);
            A(idxDrop, cEPS) = NaN;
            A(idxDrop, cETA) = NaN;
            removedSmallCC = [removedSmallCC; idxDrop(:)]; %#ok<AGROW>

            % Refresh support and adjacency after dropping small components
            isSupportGlobal = inBeta & isfinite(A(:,cEPS)) & isfinite(A(:,cETA));
            idxSupportGlobal = find(isSupportGlobal);
            [~, locSupport] = ismember(idxSupportGlobal, idxBeta);
            locSupport = locSupport(locSupport > 0);
            if isempty(locSupport)
                continue;
            end
            DS = D(locSupport, locSupport);
            nS = size(DS,1);
        end
    end

    % Isolated + high-anis gate on remaining support
    if nS == 1
        deg = 0;
    else
        AS = (DS <= neighTol);
        AS(1:nS+1:end) = false;
        deg = sum(AS, 2);
    end

    anisSupport = anis(idxSupportGlobal);
    if numel(anisSupport) >= 5
        anis_ref = local_percentile(anisSupport, 95);
    else
        anis_ref = max(anisSupport);
    end

    isIso = deg < minNeighbors;
    isHigh = anisSupport > (anis_ref * anisTolFactor);
    dropIsoHigh = isIso & isHigh;

    if any(dropIsoHigh)
        idxDrop = idxSupportGlobal(dropIsoHigh);
        A(idxDrop, cEPS) = NaN;
        A(idxDrop, cETA) = NaN;
        removedIsoHigh = [removedIsoHigh; idxDrop(:)]; %#ok<AGROW>
    end
end

% Build report
reportClean = struct();
reportClean.removed_isolated_high_anis_idx = unique(removedIsoHigh(:));
reportClean.removed_small_components_idx = unique(removedSmallCC(:));
reportClean.neighTol_per_beta = neighTol_per_beta;

% Preserve input format
if isTable
    tblClean = tbl;
    tblClean.eps_bist = A(:,cEPS);
    tblClean.eta_val = A(:,cETA);
else
    tblClean = tbl;
    tblClean(:,4) = A(:,cEPS);
    tblClean(:,5) = A(:,cETA);
end
end

function p = local_percentile(x, q)
% Toolbox-free percentile using linear interpolation on sorted values.
x = x(:);
x = x(isfinite(x));
if isempty(x)
    p = NaN;
    return;
end
x = sort(x);
n = numel(x);
if n == 1
    p = x;
    return;
end
pos = 1 + (q/100) * (n - 1);
i0 = floor(pos);
i1 = ceil(pos);
if i0 == i1
    p = x(i0);
else
    t = pos - i0;
    p = (1 - t) * x(i0) + t * x(i1);
end
end

function comp = local_conn_components(A)
% Connected components for undirected adjacency matrix A (logical), toolbox-free.
n = size(A,1);
comp = zeros(n,1);
cid = 0;
for i = 1:n
    if comp(i) ~= 0
        continue;
    end
    cid = cid + 1;
    queue = i;
    comp(i) = cid;
    head = 1;
    while head <= numel(queue)
        v = queue(head);
        head = head + 1;
        nbr = find(A(v,:));
        nbr = nbr(comp(nbr) == 0);
        if isempty(nbr)
            continue;
        end
        comp(nbr) = cid;
        queue = [queue, nbr]; %#ok<AGROW>
    end
end
end
