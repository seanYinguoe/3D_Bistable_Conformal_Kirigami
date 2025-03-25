function [f_grid, idxn] = removeId(f_grid, idE)
% Replace entire rows specified by idE with NaNs in f_grid
% f_grid: cell array (multiple rows, 3 columns)
% idE: indices of rows to replace entirely with NaNs
% idxn: indices of modified rows

idxn = []; % Initialize indices of modified rows

for i = 1:numel(idE)
    rowIdx = idE(i);
    if rowIdx <= size(f_grid, 1)
        f_grid(rowIdx, :) = NaN; % Replace entire row with NaNs
        idxn = [idxn; rowIdx]; % Record modified row index
    end
end

end
