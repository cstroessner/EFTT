function y = applyToRows(func,matrix)

applyToGivenRow = @(func, matrix) @(row) func(matrix(row, :));
y = arrayfun(applyToGivenRow(func, matrix), 1:size(matrix,1))';
end