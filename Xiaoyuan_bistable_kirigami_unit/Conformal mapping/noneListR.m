function [temp, idxD] = noneListR(a)
%{
    Remove NaN values from either a numeric list or a cell array and maintain an index mapping.
    
    a     - Numeric array or cell array containing numbers (faces, indices, etc.)
    temp  - Filtered version of `a`, with NaN values removed
    idxD  - Dictionary (Map) that maps old indices to new ones
%}

idxD = containers.Map('KeyType', 'double', 'ValueType', 'double'); % Initialize map
idx = 0; % New index counter

if iscell(a)
    % Handle case where input is a cell array
    temp = {}; % Initialize empty cell array
    for i = 1:numel(a)
        ai = a{i}; % Get current element
        if ~(iscell(ai) && all(cellfun(@(x) isequaln(x, NaN), ai))) && ~isequaln(ai, NaN)
            temp{end+1} = ai; % Append non-NaN values
            idxD(i) = idx; % Store index mapping
            idx = idx + 1;
        end
    end
    temp = temp';
else
    % Handle case where input is a numeric array
    mask = ~isnan(a); % Logical mask for non-NaN elements
    temp = a(mask); % Filtered numeric array
    idxD = find(mask) - 1; % Store original indices (MATLAB is 1-based)
end
end
