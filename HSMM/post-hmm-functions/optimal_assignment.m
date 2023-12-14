function [reorder_col_ind, col_ind, extra_ind] = optimal_assignment(hidden_states, labels)
% Wraper function that prepares a cosine similarity based cost matrix for a
% call to matchpairs.m, MATLAB's function for solving the linear assignment
% problem. Note that the optimal assignment will feature k = min(m, n)
% pairs for m states and y behaviors/labels. Number of hidden states is
% expected to be greater than or equal to number of behavioral labels.
% 
% PARAMETERS
% ----------
% hidden_states : k-vector (for k timepoints), each of whose elements is a
%                 member of {r_1, r_2, ..., r_n}, where r_i is a real
%                 number associated with the i_th state (out of n states).
% labels        : k-vector (for k timepoints), each of whose elements is a
%                 member of {q_1, q_2, ..., q_m}, where q_j is a real
%                 number associated with the j_th behavioral label (out of
%                 m behaviors).
% 
% RETURNS
% -------
% reorder_col_ind : concatenation of col_ind and extra_ind
% col_ind         : For c hidden states and d behavioral labels, col_ind is
%                   a d-vector whose i_th element is a member of {1, ...,
%                   c} indicating the index of the hidden state to which
%                   the i_th behavior is best matched.
% extra_ind       : (c-d)-vector whose elements are the indices of the
%                   hidden states not matched to a behavior.
%
% Author: Jonathan Chien


hidden_states_list = unique(hidden_states);
labels_list = unique(labels);

% Match pairs by solving linear assignment problem.
cost = nan(length(labels_list), length(hidden_states_list));
for i_label = 1:length(labels_list)
for i_state = 1:length(hidden_states_list)
    cost(i_label, i_state) ...
        = cos_sim(double(labels == labels_list(i_label)), ...
                  double(hidden_states == hidden_states_list(i_state)));
end
end
cost = (1 - cost) / 2;

% matches matrix has second column sorted by default, sort by first column
% (bhv_labels) instead.
matches = matchpairs(cost, 2);
[~, sort_ind] = sort(matches(:, 1), 'ascend');
matches = matches(sort_ind, :);

% Get reordering indices.
col_ind = matches(:, 2);
extra_ind = setdiff(1:length(hidden_states_list), col_ind);
reorder_col_ind = [col_ind' extra_ind];

end
