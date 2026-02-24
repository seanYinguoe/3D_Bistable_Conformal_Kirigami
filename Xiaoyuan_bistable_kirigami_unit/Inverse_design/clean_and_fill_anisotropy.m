function [tblFilled, report] = clean_and_fill_anisotropy(tbl)
%CLEAN_AND_FILL_ANISOTROPY Fill missing bistable values with hex-region-aware local interpolation.
%   [tblFilled, report] = clean_and_fill_anisotropy(tbl)
%   Input columns: [a1 a2 a3 eps_bist eta_val beta] (table or numeric matrix)
%
% Strategy per beta slice:
%   1) Define support = finite eps_bist & eta_val.
%   2) Define a hex-like admissible region using distance to isotropic center:
%        dev = max(abs([a1 a2 a3] - pi/3))
%      and radius from support (robust percentile).
%   3) Fill candidates only if:
%        - satisfy local anisotropy candidate rule
%        - inside hex-like region.
%   4) Prefer same-a1 lane 1D interpolation; fallback to local smooth KNN.
%   5) Iterate a few passes to fill contiguous interior nodes conservatively.

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
    report = struct('filled_idx', [], ...
                    'filled_method', {{}}, ...
                    'skipped_not_candidate_idx', [], ...
                    'skipped_no_support_idx', []);
    return;
end

% Column map
cA1 = 1; cA2 = 2; cA3 = 3; cEPS = 4; cETA = 5; cBETA = 6;

% Settings
tol = 1e-6;
a1Tol = 1e-4;
K = 8;
nGridLocal = 24;
smoothLocal = 0.55;
distTol = 0.40;       % sparse-support local gate in (a1,a2,a3)
maxPass = 6;          % iterative filling to complete compact region
hexQ = 95;            % support percentile for hex radius
hexExpand = 1.18;     % expansion to reduce interior holes
holePass = 2;         % extra passes for interior hole closing
holeMinNbr = 3;       % require enough nearby filled neighbors
lastHoleMinNbr = 5;   % stricter final single-hole closure

anis = max(A(:,[cA1 cA2 cA3]), [], 2) - min(A(:,[cA1 cA2 cA3]), [], 2);
devIso = max(abs(A(:,[cA1 cA2 cA3]) - pi/3), [], 2); % hex-like metric

filled_idx = [];
filled_method = {};
skipped_not_candidate_idx = [];
skipped_no_support_idx = [];

betaVals = unique(A(:,cBETA));

for ib = 1:numel(betaVals)
    beta0 = betaVals(ib);
    inBeta = abs(A(:,cBETA) - beta0) <= tol;

    % iterative fill passes (support updates as nodes are filled)
    for pass = 1:maxPass
        isSupport = inBeta & isfinite(A(:,cEPS)) & isfinite(A(:,cETA));
        idxSupportAll = find(isSupport);
        if numel(idxSupportAll) < 3
            break;
        end

        % Hex-like admissible radius from support
        devSupport = devIso(idxSupportAll);
        if numel(devSupport) >= 5
            rHex = local_percentile(devSupport, hexQ) * hexExpand;
        else
            rHex = max(devSupport) * 1.02;
        end

        % Optional dense-support in-shape gate in (a2,a3)
        useGlobalShape = false;
        shp = [];
        if numel(idxSupportAll) >= 18 && exist('alphaShape', 'class') == 8
            try
                xAll = A(idxSupportAll,cA2);
                yAll = A(idxSupportAll,cA3);
                shp = alphaShape(xAll, yAll);
                useGlobalShape = true;
            catch
                useGlobalShape = false;
            end
        end

        isMissing = inBeta & (~isfinite(A(:,cEPS)) | ~isfinite(A(:,cETA)));
        idxMiss = find(isMissing);
        if isempty(idxMiss)
            break;
        end

        nFilledThisPass = 0;

        for im = 1:numel(idxMiss)
            iRow = idxMiss(im);
            a1q = A(iRow,cA1);
            a2q = A(iRow,cA2);
            a3q = A(iRow,cA3);
            anisQ = anis(iRow);

            % Hex-like region gate
            if devIso(iRow) > (rHex + tol)
                skipped_not_candidate_idx = [skipped_not_candidate_idx; iRow]; %#ok<AGROW>
                continue;
            end

            % Same-a1 local support for candidate rule and lane interpolation
            isLane = isSupport & (abs(A(:,cA1) - a1q) <= a1Tol);
            idxLane = find(isLane);

            if ~isempty(idxLane)
                anisRef = anis(idxLane);
            else
                dAllRef = sqrt((A(idxSupportAll,cA1)-a1q).^2 + ...
                               (A(idxSupportAll,cA2)-a2q).^2 + ...
                               (A(idxSupportAll,cA3)-a3q).^2);
                [~, ordRef] = sort(dAllRef, 'ascend');
                kRef = min(K, numel(idxSupportAll));
                idxRef = idxSupportAll(ordRef(1:kRef));
                anisRef = anis(idxRef);
            end

            % Local candidate rule
            if isempty(anisRef) || ~(anisQ <= max(anisRef) + tol)
                skipped_not_candidate_idx = [skipped_not_candidate_idx; iRow]; %#ok<AGROW>
                continue;
            end

            didFill = false;

            % ---------------- Priority 1: same-a1 lane1D ----------------
            if numel(idxLane) >= 2
                a2Lane = A(idxLane,cA2);
                lower = find(a2Lane < (a2q - tol));
                upper = find(a2Lane > (a2q + tol));

                if ~isempty(lower) && ~isempty(upper)
                    [~, il] = max(a2Lane(lower));
                    [~, iu] = min(a2Lane(upper));
                    iL = idxLane(lower(il));
                    iU = idxLane(upper(iu));

                    a2L = A(iL,cA2); a2U = A(iU,cA2);
                    if abs(a2U - a2L) > eps
                        w = (a2q - a2L) / (a2U - a2L);
                        ePred = A(iL,cEPS) + w * (A(iU,cEPS) - A(iL,cEPS));
                        hPred = A(iL,cETA) + w * (A(iU,cETA) - A(iL,cETA));

                        % conservative clamp to lane local range
                        eLoc = A(idxLane,cEPS);
                        hLoc = A(idxLane,cETA);
                        eMin = min(eLoc); eMax = max(eLoc);
                        hMin = min(hLoc); hMax = max(hLoc);
                        ePad = max(1e-8, 0.05*(eMax - eMin));
                        hPad = max(1e-8, 0.05*(hMax - hMin));
                        ePred = min(max(ePred, eMin - ePad), eMax + ePad);
                        hPred = min(max(hPred, hMin - hPad), hMax + hPad);

                        if isfinite(ePred) && isfinite(hPred)
                            A(iRow,cEPS) = max(ePred, 0);
                            A(iRow,cETA) = max(hPred, 0);
                            filled_idx = [filled_idx; iRow]; %#ok<AGROW>
                            filled_method{end+1,1} = 'lane1D'; %#ok<AGROW>
                            didFill = true;
                            nFilledThisPass = nFilledThisPass + 1;
                        end
                    end
                end
            end

            if didFill
                continue;
            end

            % ------------- Priority 2: local smooth fallback -------------
            dAll = sqrt((A(idxSupportAll,cA1)-a1q).^2 + ...
                        (A(idxSupportAll,cA2)-a2q).^2 + ...
                        (A(idxSupportAll,cA3)-a3q).^2);
            [dSort, ordAll] = sort(dAll, 'ascend');
            kEff = min(K, numel(idxSupportAll));
            nb = idxSupportAll(ordAll(1:kEff));

            if useGlobalShape
                try
                    inShp = inShape(shp, a2q, a3q);
                catch
                    inShp = true;
                end
                if ~inShp
                    skipped_no_support_idx = [skipped_no_support_idx; iRow]; %#ok<AGROW>
                    continue;
                end
            else
                if isempty(dSort) || dSort(1) > distTol
                    skipped_no_support_idx = [skipped_no_support_idx; iRow]; %#ok<AGROW>
                    continue;
                end
            end

            if isempty(nb)
                skipped_no_support_idx = [skipped_no_support_idx; iRow]; %#ok<AGROW>
                continue;
            end

            Tlocal = table(A(nb,cA1), A(nb,cA2), A(nb,cA3), A(nb,cEPS), A(nb,cETA), ...
                'VariableNames', {'a1','a2','a3','eps_bist','eta_val'});

            etaLocal = Tlocal.eta_val(isfinite(Tlocal.eta_val));
            if isempty(etaLocal)
                skipped_no_support_idx = [skipped_no_support_idx; iRow]; %#ok<AGROW>
                continue;
            end
            etaThrLocal = min(etaLocal) - 1e-12;

            Tdense = interpolate_ternary_smooth(Tlocal, nGridLocal, smoothLocal, etaThrLocal);
            if isempty(Tdense)
                skipped_no_support_idx = [skipped_no_support_idx; iRow]; %#ok<AGROW>
                continue;
            end

            dDense = sqrt((Tdense.a1 - a1q).^2 + ...
                          (Tdense.a2 - a2q).^2 + ...
                          (Tdense.a3 - a3q).^2);
            [~, iNear] = min(dDense);

            ePred = Tdense.eps_bist(iNear);
            hPred = Tdense.eta_val(iNear);
            if ~(isfinite(ePred) && isfinite(hPred))
                skipped_no_support_idx = [skipped_no_support_idx; iRow]; %#ok<AGROW>
                continue;
            end

            % conservative local clamp
            eLoc = Tlocal.eps_bist; hLoc = Tlocal.eta_val;
            eMin = min(eLoc); eMax = max(eLoc);
            hMin = min(hLoc); hMax = max(hLoc);
            ePad = max(1e-8, 0.05*(eMax - eMin));
            hPad = max(1e-8, 0.05*(hMax - hMin));
            ePred = min(max(ePred, eMin - ePad), eMax + ePad);
            hPred = min(max(hPred, hMin - hPad), hMax + hPad);

            A(iRow,cEPS) = max(ePred, 0);
            A(iRow,cETA) = max(hPred, 0);
            filled_idx = [filled_idx; iRow]; %#ok<AGROW>
            filled_method{end+1,1} = 'localSmooth'; %#ok<AGROW>
            nFilledThisPass = nFilledThisPass + 1;
        end

        if nFilledThisPass == 0
            break;
        end
    end

    % Final conservative interior-hole fill:
    % fill only NaNs that are surrounded by nearby filled supports.
    idxBeta = find(inBeta);
    if numel(idxBeta) >= 4
        Xbeta = A(idxBeta, [cA1 cA2 cA3]);
        G = sum(Xbeta.^2, 2);
        D2 = max(G + G' - 2*(Xbeta*Xbeta'), 0);
        Dbeta = sqrt(D2);
        Dbeta(1:numel(idxBeta)+1:end) = inf;
        dMinBeta = min(Dbeta, [], 2);
        dMinBeta = dMinBeta(isfinite(dMinBeta) & dMinBeta > 0);
        if ~isempty(dMinBeta)
            holeTol = 1.35 * median(dMinBeta);
        else
            holeTol = distTol;
        end

        for hp = 1:holePass
            isSupport = inBeta & isfinite(A(:,cEPS)) & isfinite(A(:,cETA));
            idxSupportAll = find(isSupport);
            if numel(idxSupportAll) < holeMinNbr
                break;
            end

            devSupport = devIso(idxSupportAll);
            if numel(devSupport) >= 5
                rHex = local_percentile(devSupport, hexQ) * hexExpand;
            else
                rHex = max(devSupport) * 1.02;
            end

            isMissing = inBeta & (~isfinite(A(:,cEPS)) | ~isfinite(A(:,cETA)));
            idxMiss = find(isMissing);
            if isempty(idxMiss)
                break;
            end

            nHoleFilled = 0;
            for im = 1:numel(idxMiss)
                iRow = idxMiss(im);
                if devIso(iRow) > (rHex + tol)
                    continue;
                end

                d = sqrt((A(idxSupportAll,cA1)-A(iRow,cA1)).^2 + ...
                         (A(idxSupportAll,cA2)-A(iRow,cA2)).^2 + ...
                         (A(idxSupportAll,cA3)-A(iRow,cA3)).^2);
                nb = find(d <= holeTol);
                if numel(nb) < holeMinNbr
                    continue;
                end

                idxNb = idxSupportAll(nb);
                eLoc = A(idxNb,cEPS);
                hLoc = A(idxNb,cETA);
                w = 1 ./ (d(nb) + 1e-12);
                w = w / sum(w);

                ePred = sum(w .* eLoc);
                hPred = sum(w .* hLoc);

                eMin = min(eLoc); eMax = max(eLoc);
                hMin = min(hLoc); hMax = max(hLoc);
                ePad = max(1e-8, 0.03*(eMax - eMin));
                hPad = max(1e-8, 0.03*(hMax - hMin));
                ePred = min(max(ePred, eMin - ePad), eMax + ePad);
                hPred = min(max(hPred, hMin - hPad), hMax + hPad);

                A(iRow,cEPS) = max(ePred, 0);
                A(iRow,cETA) = max(hPred, 0);
                filled_idx = [filled_idx; iRow]; %#ok<AGROW>
                filled_method{end+1,1} = 'holeFill'; %#ok<AGROW>
                nHoleFilled = nHoleFilled + 1;
            end

            if nHoleFilled == 0
                break;
            end
        end

        % Final strict one-node hole closure:
        % fill only deeply interior NaNs that are strongly surrounded.
        isSupport = inBeta & isfinite(A(:,cEPS)) & isfinite(A(:,cETA));
        idxSupportAll = find(isSupport);
        if numel(idxSupportAll) >= lastHoleMinNbr
            devSupport = devIso(idxSupportAll);
            if numel(devSupport) >= 5
                rHex = local_percentile(devSupport, hexQ) * hexExpand;
            else
                rHex = max(devSupport) * 1.02;
            end

            isMissing = inBeta & (~isfinite(A(:,cEPS)) | ~isfinite(A(:,cETA)));
            idxMiss = find(isMissing);
            holeTol2 = 1.60 * holeTol;

            for im = 1:numel(idxMiss)
                iRow = idxMiss(im);

                % Strictly interior only (avoid boundary growth)
                if devIso(iRow) > (0.92 * rHex + tol)
                    continue;
                end

                d = sqrt((A(idxSupportAll,cA1)-A(iRow,cA1)).^2 + ...
                         (A(idxSupportAll,cA2)-A(iRow,cA2)).^2 + ...
                         (A(idxSupportAll,cA3)-A(iRow,cA3)).^2);
                nb = find(d <= holeTol2);
                if numel(nb) < lastHoleMinNbr
                    continue;
                end

                idxNb = idxSupportAll(nb);
                eLoc = A(idxNb,cEPS);
                hLoc = A(idxNb,cETA);
                w = 1 ./ (d(nb) + 1e-12);
                w = w / sum(w);

                ePred = sum(w .* eLoc);
                hPred = sum(w .* hLoc);

                eMin = min(eLoc); eMax = max(eLoc);
                hMin = min(hLoc); hMax = max(hLoc);
                ePad = max(1e-8, 0.02*(eMax - eMin));
                hPad = max(1e-8, 0.02*(hMax - hMin));
                ePred = min(max(ePred, eMin - ePad), eMax + ePad);
                hPred = min(max(hPred, hMin - hPad), hMax + hPad);

                A(iRow,cEPS) = max(ePred, 0);
                A(iRow,cETA) = max(hPred, 0);
                filled_idx = [filled_idx; iRow]; %#ok<AGROW>
                filled_method{end+1,1} = 'holeFillFinal'; %#ok<AGROW>
            end
        end
    end
end

report = struct();
report.filled_idx = filled_idx(:);
report.filled_method = filled_method;
report.skipped_not_candidate_idx = unique(skipped_not_candidate_idx(:));
report.skipped_no_support_idx = unique(skipped_no_support_idx(:));

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
