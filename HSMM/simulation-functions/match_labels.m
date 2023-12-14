function [row_ind, col_ind, cost] = match_labels(cost_mat)
% Greedy algorithm to find bijection between two finite sets of cardinality
% n, based on an n x n cost matrix.
% 
% 
% PARAMETERS
% ----------
% cost_mat : n x n cost matrix, where the rows index the vertices of the
%            first set, and the columns the vertices of the second. The
%            i_th, j_th element is then the cost of matching the i_th
%            vertex from the first set to the j_th vertex of the second
%            set. If set 1 consists of the true behavioral states, and set
%            2 the inferred states, the i_th, j_th element is then the cost
%            of matching the i_th behavioral state to the j_th inferred
%            state.
%
% RETURNS
% -------
% row_ind : Vector whose elements consist of 1:n, where n is defined as
%           above. These are the sorted (ascending) row indices of the
%           optimal pairs. They exist in an ordered one to one
%           correspondence with elements of col_ind.
% col_ind : Vector whose elements consist of the column indices of the
%           optimal pairing, ordered to match the elements of row_ind one
%           to one.
% cost    : Scalar equal to the sum of the cost of each of the optimal
%           pairings, giving one overall cost for the entire match.
% 
% Author: Jonathan Chien

n_pairs = size(cost_mat, 1);

pair_ind = NaN(n_pairs, 1);
cost_by_pair = NaN(n_pairs, 1);
for i_pair = 1:n_pairs
    % Get optimal pairing among options.
    [cost_by_pair(i_pair), pair_ind(i_pair)] = min(cost_mat, [], 'all');
    
    % Remove the pairing just found and repeat. 
    [i_row, i_col] = ind2sub([n_pairs n_pairs], pair_ind(i_pair));
    cost_mat = cover_pair(cost_mat, i_row, i_col);
end

% Get optimal pairing.
[row_ind, col_ind] = ind2sub([n_pairs n_pairs], pair_ind);

% Sort row labels and apply same order to column labels.
[row_ind, sort_ind] = sort(row_ind);
col_ind = col_ind(sort_ind);

% Calculate total cost.
cost = sum(cost_by_pair);

end


% --------------------------------------------------
function cost_mat = cover_pair(cost_mat, i_row, i_col)
% Takes in an n x n cost matrix and converts the k_th row and column to NaN
% values.

% n_pairs = n for n x n cost matrix.
n_pairs = size(cost_mat, 1);

% Assertion not currently necessary but may be helpful if this local
% routine is ever made its own function.
assert(i_row <= n_pairs);
assert(i_col <= n_pairs);

% Subscripts of row to be covered.
row.row_ind = repelem(i_row, n_pairs);
row.col_ind = 1:n_pairs;

% Subscripts of column to be covered.
col.row_ind = setdiff(1:n_pairs, i_row);
col.col_ind = repelem(i_col, n_pairs-1);

% Convert subscripts to linear indices.
linear_ind = [sub2ind([n_pairs n_pairs], row.row_ind, row.col_ind) ...
              sub2ind([n_pairs n_pairs], col.row_ind, col.col_ind)];

% Set to NaN and reshape to recover original array shape.
cost_mat(linear_ind) = NaN;
cost_mat = reshape(cost_mat, n_pairs, n_pairs);

end
