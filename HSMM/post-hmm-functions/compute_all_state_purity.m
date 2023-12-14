function [purity, purity_by_bhv] ...
    = compute_all_state_purity(row_membership, occupancy, skip_final)
% Computes purity scores for the case where all hidden states are matched
% to a behavior (the case where the number of hidden states and labels is
% equal should be a subcase but have not tested).
%
% PARAMETERS
% ----------
% row_membership : For n hidden states, this is an n element vector whose
%                  j_th element is a member of {1, ..., n} denoting which
%                  behavioral label (hence which row of the occupancy
%                  matrix) the j_th hidden state is matched to.
% occupancy      : m x n occupancy matrix. The i_th, j_th element is the
%                  number of timepoints of the i_th behavior assigned to
%                  the j_th hidden state.
% skip_final     : {1 | 0}, if true, omits the final behavior when
%                  computing the purity score. TODO: this was a quick and
%                  hacky method to exclude the "Other" condition by placing
%                  it last); it should perhaps be replaced with something
%                  like a parameter such as "i_other" that specifies the
%                  index of the "Other" label etc.
% 
% RETURNS
% -------
% purity        : Weighted mean, based on occupancy of each behavior, of
%                 purity scores for each behavior (see below).
% purity_by_bhv : Vector of length c, where c = number of behaviors. The
%                 i_th element is the purity score for the i_th behavior
%                 based on the combined column purities of all hidden
%                 states matched to this behavior. More specifically, for
%                 the i_th behavior, the column purity of each of the
%                 states matched to this behavior is first computed (for
%                 each state matched to the i_th behavior, this is the
%                 overlap count for the hidden state and the behavior,
%                 divided by the total number of timepoints in that hidden
%                 state; e.g., if states 5, 6, 7, and 8 were matched to
%                 behavior 3, we would here take the elements of the
%                 occupancy matrix indexed by (3,5:8) and divide them
%                 elementwise by the sums of columns 5-8 of the occupancy
%                 matrix). Next, we combine these individual purity
%                 scores into one score for the behavior by taking a
%                 weighted mean based on occupancy in each of the states
%                 mathced to the i_th bhv (weight of the j_th state is n_j
%                 / m, where n_j is the number of timepoints overlapping
%                 between the j_th state and the i_th behavior, and m is
%                 the sum across n_j for all hidden states matched to the
%                 i_th behavior).   
%
% Author: Jonathan Chien


% Other must be placed last. Then if skip_final = true/1, this will cause
% it to be omitted from purity calculations.
assert(ismember(skip_final, [0 1]));

% Purity scores based on greedy assignment.
purity_by_bhv = nan(size(occupancy, 1) - skip_final, 1);
for i_bhv = 1:(size(occupancy, 1) - skip_final)
    % Get count of each hidden state matched to current behavior.
    occupancy_by_state_within_bhv = occupancy(i_bhv, row_membership == i_bhv);

    % Get sum of timepoints in each hidden state (across behaviors) matched
    % to current behavior.
    occupancy_by_state_across_bhv ...
        = sum(occupancy(1:end-skip_final, row_membership == i_bhv), 1);
    
    % Compute the column purity of each hidden state that was matched with
    % the current behavior.
    within_bhv_purity ...
        = occupancy_by_state_within_bhv ./ occupancy_by_state_across_bhv;

    % Take weighted mean of the column purity of each hidden state (weight
    % of the i_th state is n_i / m, where n_i is the number of timepoints
    % overlapping between that state and the current behavior, and m is the
    % sum across n_i for all hidden states matched to the current behavior.
    purity_by_curr_bhv ...
        = sum(ensure_col(within_bhv_purity) ...
              .* ensure_col(occupancy_by_state_within_bhv / sum(occupancy_by_state_within_bhv)));
    
    if isempty(purity_by_curr_bhv), purity_by_curr_bhv = nan; end

    % Assign into container.
    purity_by_bhv(i_bhv) = purity_by_curr_bhv;
end

% Compute final purity score.
purity = sum(ensure_col(purity_by_bhv) ...
             .* (sum(occupancy(1:end-skip_final,:), 2) ...
                 ./ sum(occupancy(1:end-skip_final,:), 'all')), ...
             'omitnan');

end
