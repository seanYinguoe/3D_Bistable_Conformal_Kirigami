function weights = normalizeWeights(weights)
% Normalize the weight values so they sum to 1.
% If all weights are zero, return equal weights.
% 
% Input:
% - weights: Vector of unnormalized weights
% 
% Output:
% - weights: Normalized weights summing to 1
% %

    % Avoid division by zero by checking if the sum is nonzero
    total = sum(weights);
    if total > 0
        weights = weights / total;
    else
        % If all weights are zero, assign equal weights
        weights = ones(size(weights)) / numel(weights);
    end
end
