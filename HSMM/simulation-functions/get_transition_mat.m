function transition_mat ...
    = get_transition_mat(bhv_labels, bhvs_to_exclude, bhv_name_codes)
% Generate ground truth transition matrix based on empirical behavioral
% timeseries data.
%
%
% PARAMETERS
% ----------
% bhv_labels      : n_timepoints x 1 or 1 x n_timepoints array whose i_th
%                   element is the behavioral label of the i_th timepoint.
% bhvs_to_exclude : Cell array whose i_th element contains the name of one
%                   of the behaviors that we would like to exclude. Note
%                   that it is possible to include behaviors that are not
%                   actually present in a given session's data (they will
%                   simply be ignored). 
% bhv_name_codes  : Vector whose i_th element is the number code associated
%                   with the i_th element of bhvs_to_exclude. As with
%                   bhvs_to_exclude, extra codes will be ignored.
%
% RETURNS
% -------
% transition_mat : An n_bhvs x n_bhvs transition matrix whose i_th,
%                  j_th element denotes the probability of transitioning
%                  from behavior i to behavior j.
%
% Author: Jonathan Chien


% Process inputs/set up matrix indexing behaviors across time. Order of
% unique_bhvs matches that of aux_bhv_info.bhv_list for the corresponding
% session.
unique_bhvs = unique(bhv_labels, 'sorted');
n_bhvs = length(unique_bhvs);
if isrow(bhv_labels), bhv_labels = bhv_labels'; end
if iscolumn(unique_bhvs), unique_bhvs = unique_bhvs'; end
bhv_ind = bhv_labels == unique_bhvs;

% Preallocate counts matrix.
counts = NaN(n_bhvs, n_bhvs);

% Fill counts matrix.
for i_bhv = 1:n_bhvs
    % Get indices marking occurence of current behavior.
    curr_bhv_ind = bhv_ind(:,i_bhv);

    % If this behavior was the last one annotated (a special edge case),
    % set last element to false before applying circshift.
    if curr_bhv_ind(end) == 1, curr_bhv_ind(end) = false; end
    
    % Circshift by 1 means that indices now identify timepoints following
    % the timepoints where behavior of interest occured.
    next_ind = circshift(curr_bhv_ind, 1);

    % Get rows of bhv_ind indexed by next_ind, use find function (returning
    % linear indices on non-vector/scalar arrays), and convert back to
    % subscripts (need column subscript).
    [~, col_ind] = ind2sub([sum(next_ind==1) n_bhvs], find(bhv_ind(next_ind,:)));

    counts(i_bhv,:) = histcounts(col_ind, 1:(n_bhvs + 1));
end

% Get column indices corresponding to behaviors to be excluded.
col_ind = bhv_name_to_remove_ind(bhvs_to_exclude, bhv_name_codes, unique_bhvs);

% Convert column indices of behaviors to be removed to linear indices of
% these columns (and their corresponding rows) and delete these matrix
% elements.
linear_ind = col_ind_to_row_col_linear_ind(col_ind, n_bhvs);
if ~isempty(linear_ind) % linear_ind will be empty if no bhv was to be excluded or that behavior not present in a given session
    counts(linear_ind) = [];
    counts = reshape(counts, [sqrt(length(counts)), sqrt(length(counts))]);
end

% Normalize so that rows sum to 1.
transition_mat = counts ./ sum(counts, 2);

end


% --------------------------------------------------
function linear_ind = col_ind_to_row_col_linear_ind(col_ind, n)
% For a given column index, i, in {1, 2, ..., n}, where n is the
% total number of behaviors, we want to remove both the i_th
% column and the i_th row from the count matrix. Given i, (or a set
% of column indices, i, ii, iii, ...) this subroutine finds and
% returns the linear indices corresponding to elements of the i_th
% row and i_th column.

% n_cols = number of columns/behaviors to exclude.
n_cols = length(col_ind);

% Get linear indices for each behavior.
linear_ind = [];
for i_col = 1:n_cols
    % i is the column index for current behavior.
    i = col_ind(i_col);

    % Prepare row and column subscripts for the i_th row.
    row.row_subscr = repelem(i, n);
    row.col_subscr = 1:n;

    % Prepare row and column subscripts for the i_th column. Note
    % that we use setdiff/n-1 to avoid a repeat index for the
    % intersection of the i_th row and i_th column (the i_th
    % diagonal element). Also could just grab indices as before and
    % call union on resulting linear indices to remove redundancy.
    col.row_subscr = setdiff(1:n, i);
    col.col_subscr = repelem(i, n-1);

    % Convert to linear indices.
    row.linear_ind = sub2ind([n n], row.row_subscr, row.col_subscr);
    col.linear_ind = sub2ind([n n], col.row_subscr, col.col_subscr);

    % Concatenate to growing container of linear indices.
    linear_ind = unique([linear_ind row.linear_ind col.linear_ind]);
end

end
