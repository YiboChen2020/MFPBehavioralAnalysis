function [col_ind, row_membership, row_membership_unsorted, max_overlaps_2, B_] = all_state_match(x, y)
% For a set of hidden states and behavioral labels, all hidden states are
% matched to a behavioral label (useful when the number of hidden states
% exceeds the number of behavioral labels).
%
% PARAMETERS
% ----------
% x : k-vector, each of whose elements is a member of {r_1, r_2, ...,
%     r_n}, where r_i is a real number associated with the i_th behavior
%     (out of n behaviors).
% y : k-vector, each of whose elements is a member of {q_1, q_2, ...,
%     q_m}, where q_j is a real number associated with the j_th state (out
%     of m states).
%
% 
% RETURNS
% -------
% col_ind                 : For d hidden states, this is a reordering
%                           vector  of length d. The p_th element is the
%                           index of the hidden state (under the original
%                           ordering based on ascending numeric order of
%                           the states) that will be assigned to the p_th
%                           element of the hidden state list under the new
%                           ordering. This ordering rearranges the hidden
%                           states so that all hidden states matched with
%                           behavioral label 1 come first, followed by all
%                           hidden states matched with behavioral label 2,
%                           then all states matched with behavioral label 3
%                           etc. All hidden states matched with a single
%                           behavioral label are reordered according to the
%                           strength of their overlap with that label (in
%                           descending order). Note that behavioral labels
%                           don't necessarily need to be {1, 2, ...}. As
%                           long as they are real numbers, the function
%                           will operate as described based on the
%                           ascending numeric order of the behavioral
%                           labels.
% row_membership          : For d hidden states, this is a d-vector whose
%                           p_th element is the index of the behavior to
%                           which the p_th state is best matched. Note that
%                           this is under the reordering encoded by col_ind
%                           above.
% row_membership_unsorted : For d hidden states, this is a d-vector whose
%                           p_th element is the index of the behavior to
%                           which the p_th state is best matched. Note that
%                           this is for the hidden states sorted based on
%                           the ascending numeric order of the hidden
%                           states.
% max_overlaps_2          : For d hidden states, this is a d-vector whose
%                           p_th element is the normalized overlap
%                           between the p_th hidden state (under the
%                           reordering encoded in col_ind) and the
%                           behavioral label with which this state is best
%                           matched (normalized count is the raw count
%                           normalized by the number of occurences of the
%                           best matched behavioral label).
% B_                      : Overlap matrix (by count), reordered according
%                           to the optimal ordering in col_ind.
%
% Author: Jonathan Chien


% Get overlap matrix and for each y state determine the x state with which
% it maximally overlaps. max_ind is the prototypical row_membership since
% its values are the x state that is the best match for each y state.
B = get_overlap_mat(x, y, 'count');
B_norm = B ./ sum(B, 2);

% For each column, find row with max overlap; note overlap and row index.
[max_overlaps, max_ind] = max(B_norm, [], 1);
row_membership_unsorted = max_ind;

% Two rounds of sorting. First groups y states so that they are grouped by
% the x state with which they best match. 
[row_membership, sort_ind_1] = sort(max_ind, 'ascend');
max_overlaps_1 = max_overlaps(sort_ind_1);

% The second round sorts within groups based on overlap with the matched x
% state.
max_overlaps_2 = [];
sort_ind_2 = [];
for i_x = 1:length(unique(x))
    [max_vals, max_ind] = sort(max_overlaps_1(row_membership == i_x), 'descend');
    max_overlaps_2 = [max_overlaps_2, max_vals];
    sort_ind_2 = [sort_ind_2, length(sort_ind_2) + max_ind];
end

% Initialize column indices. Apply reorderings to col_ind to get final
% ordering.
col_ind = 1:length(y);
col_ind = col_ind(sort_ind_1);
col_ind = col_ind(sort_ind_2);

% Apply reordering to overlap matrix.
B_ = B(:, col_ind);

end
