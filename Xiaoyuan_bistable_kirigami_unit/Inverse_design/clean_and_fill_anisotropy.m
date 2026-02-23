function [tblFilled, report] = clean_and_fill_anisotropy(tbl)
%CLEAN_AND_FILL_ANISOTROPY Clean outliers and locally fill missing bistable nodes.
%   [tblFilled, report] = clean_and_fill_anisotropy(tbl)
%   tbl columns: [a1 a2 a3 eps_bist eta_val beta] (table or numeric matrix)

isTable = istable(tbl);

if isTable
    req = {'a1','a2','a3','eps_bist','eta_val','beta'};
    if ~all(ismember(req, tbl.Properties.VariableNames))
        error('Input table must contain variables: a1,a2,a3,eps_bist,eta_val,beta');
    end
    A = [tbl.a1, tbl.a2, tbl.a3, tbl.eps_bist, tbl.eta_val, tbl.beta];
else
    if ~isnumeric(tbl) || size(tbl,2) < 6
        error('Numeric input must be Nx6+ with columns [a1 a2 a3 eps_bist eta_val beta].');
    end
    A = tbl(:,1:6);
end

n = size(A,1);
if n == 0
    tblFilled = tbl;
    report = struct('removed_outliers_idx', [], 'filled_idx', [], ...
        'skipped_high_anis_idx', [], 'skipped_outside_hull_idx', []);
    return;
end

% Column map
cA1 = 1; cA2 = 2; cA3 = 3; cEPS = 4; cETA = 5; cBETA = 6;

% Settings
angTol = 1e-6;
K = 8;
nGridLocal = 28;
smoothLocal = 0.6;
outlierTol = 1.15; % 15%

% Precompute anisotropy metric for all rows
anis = max(A(:,[cA1 cA2 cA3]), [], 2) - min(A(:,[cA1 cA2 cA3]), [], 2);

removed_outliers_idx = [];
filled_idx = [];
skipped_high_anis_idx = [];
skipped_outside_hull_idx = [];

betaVals = unique(A(:,cBETA));

for ib = 1:numel(betaVals)
    beta0 = betaVals(ib);
    inBeta = abs(A(:,cBETA) - beta0) <= angTol;

    % Existing bistable rows (ground truth)
    isB = inBeta & isfinite(A(:,cEPS)) & isfinite(A(:,cETA));
    idxB = find(isB);
    if isempty(idxB)
        % No support to fill from in this beta slice
        continue;
    end

    anis_max_bistable = max(anis(idxB));
    if ~isfinite(anis_max_bistable)
        continue;
    end

    % Rule 3: mark far-from-center rows as outliers, keep NaN permanently
    isOut = inBeta & (anis > outlierTol * anis_max_bistable);
    idxOut = find(isOut);
    if ~isempty(idxOut)
        A(idxOut,cEPS) = NaN;
        A(idxOut,cETA) = NaN;
        removed_outliers_idx = [removed_outliers_idx; idxOut(:)]; %#ok<AGROW>
    end

    % Refresh support after outlier marking
    isSupport = inBeta & isfinite(A(:,cEPS)) & isfinite(A(:,cETA));
    idxSupport = find(isSupport);
    if numel(idxSupport) < 3
        continue;
    end

    % Build interpolation domain in (a2,a3)
    xS = A(idxSupport,cA2);
    yS = A(idxSupport,cA3);
    DT = delaunayTriangulation(xS, yS);
    hullIdx = convhull(xS, yS);

    % Missing rows in this beta slice
    isMissing = inBeta & (~isfinite(A(:,cEPS)) | ~isfinite(A(:,cETA)));
    idxMiss = find(isMissing);

    for im = 1:numel(idxMiss)
        iRow = idxMiss(im);

        % Rule 4: only candidate bistable if anis <= anis_max_bistable
        if anis(iRow) > anis_max_bistable + angTol
            skipped_high_anis_idx = [skipped_high_anis_idx; iRow]; %#ok<AGROW>
            continue;
        end

        xq = A(iRow,cA2);
        yq = A(iRow,cA3);

        % Rule 5: no extrapolation outside support domain
        triId = pointLocation(DT, xq, yq);
        if isnan(triId)
            inHull = inpolygon(xq, yq, xS(hullIdx), yS(hullIdx));
            if ~inHull
                skipped_outside_hull_idx = [skipped_outside_hull_idx; iRow]; %#ok<AGROW>
                continue;
            end
        end

        % Local nearest neighbors in (a1,a2,a3)
        d = sqrt((A(idxSupport,cA1)-A(iRow,cA1)).^2 + ...
                 (A(idxSupport,cA2)-A(iRow,cA2)).^2 + ...
                 (A(idxSupport,cA3)-A(iRow,cA3)).^2);
        [~, ord] = sort(d, 'ascend');
        kEff = min(K, numel(idxSupport));
        nb = idxSupport(ord(1:kEff));

        Tlocal = table(A(nb,cA1), A(nb,cA2), A(nb,cA3), A(nb,cEPS), A(nb,cETA), ...
            'VariableNames', {'a1','a2','a3','eps_bist','eta_val'});

        % Keep all local support points; do not impose external eta threshold
        etaLocal = Tlocal.eta_val(isfinite(Tlocal.eta_val));
        if isempty(etaLocal)
            continue;
        end
        etaThrLocal = min(etaLocal) - 1e-12;

        Tdense = interpolate_ternary_smooth(Tlocal, nGridLocal, smoothLocal, etaThrLocal);
        if isempty(Tdense)
            continue;
        end

        % Use nearest dense prediction to requested point
        dDense = sqrt((Tdense.a1 - A(iRow,cA1)).^2 + ...
                      (Tdense.a2 - A(iRow,cA2)).^2 + ...
                      (Tdense.a3 - A(iRow,cA3)).^2);
        [~, iNear] = min(dDense);

        ePred = Tdense.eps_bist(iNear);
        hPred = Tdense.eta_val(iNear);

        % Rule 6: sanity constraints
        if ~(isfinite(ePred) && isfinite(hPred))
            continue;
        end
        if ePred < 0 || hPred < 0
            continue;
        end

        % Discard extreme predictions relative to local neighborhood spread
        eMax = max(Tlocal.eps_bist);
        hMax = max(Tlocal.eta_val);
        if ePred > (1.5 * max(eMax, 0) + 1e-12) || hPred > (1.5 * max(hMax, 0) + 1e-12)
            continue;
        end

        A(iRow,cEPS) = ePred;
        A(iRow,cETA) = hPred;
        filled_idx = [filled_idx; iRow]; %#ok<AGROW>
    end
end

% Build report
report = struct();
report.removed_outliers_idx = unique(removed_outliers_idx(:));
report.filled_idx = unique(filled_idx(:));
report.skipped_high_anis_idx = unique(skipped_high_anis_idx(:));
report.skipped_outside_hull_idx = unique(skipped_outside_hull_idx(:));

% Preserve input format
if isTable
    tblFilled = tbl;
    tblFilled.eps_bist = A(:,cEPS);
    tblFilled.eta_val = A(:,cETA);
else
    tblFilled = tbl;
    tblFilled(:,4) = A(:,cEPS);
    tblFilled(:,5) = A(:,cETA);
end
end
