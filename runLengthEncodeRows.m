function [lengths, values] = runLengthEncodeRows(data)
% from https://stackoverflow.com/questions/28791773/find-number-of-consecutive-ones-in-binary-array

n_rows = size(data, 1);
n_cols = size(data, 2);

lengths = cell(n_rows, 1);
values = cell(n_rows, 1);
for row_idx = 1:n_rows
    startPos = find(diff([data(row_idx, 1)-1, data(row_idx, :)]));
    lengths{row_idx} = diff([startPos, n_cols + 1]);
    values{row_idx} = data(row_idx, startPos);
end