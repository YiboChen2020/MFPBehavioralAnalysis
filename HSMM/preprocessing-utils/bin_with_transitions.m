function [X_binned, bin_labels, transitions, mixed_transitions] ...
    = bin_with_transitions(X, labels, bin_width, step_size, nv)
%
% Bin 2D data along the second dim (can be easily extended to binning along
% any dim of an n-D array) where each column of the array has a
% corresponding label. This function keeps track of transitions between
% labels, noting the bin index (as well as other relevant info) of where
% the transitions occured/what labels are involved. Will also ensure that
% only one transition occurs in each bin (exception is thrown otherwise).
%
% PARAMETERS
% ----------
% X         : m x n array of data to be binned along the second dimension.
% labels    : n-vector of labels, i_th label corresponds to i_th column of
%             X.
% bin_width : Scalar size of each bin along the second dimension of X.
% step_size : Scalar size of difference between the column indices (in X)
%             of the first element of the i_th bin and the first element of
%             the i-1_th bin.
% Name-value pair(s)
%   'make_call' : {true (default) | false}. % If a transition is contained
%                 within a bin, optionally assign the whole bin to the
%                 state/behavior that comprises the majority of timepoints
%                 in the bin (recall that there should be only two states
%                 at most within a bin, as will be checked below). In the
%                 event of a tie, assign the bin to the state from which we
%                 transitioned.
%
% RETURNS
% -------
% X_binned          : m x k array of data (X binned along the 2nd dim).
% bin_labels        : k-vector of labels whose i_th element is the label of
%                     the i_th bin.
% transitions       : j x 5 array whose i_th row corresponds to the i_th
%                     transition (an event where consecutive labels
%                     differ). In the i_th row, the first ands second
%                     column elements are respectively the state labels
%                     from and to which the i_th transition occured. The
%                     third column element is a boolean, true if the
%                     timepoint preceding the transition coincides with a
%                     bin end (the state transition happens to align with a
%                     bin transition) and false otherwise (the state
%                     transition is wholly contained by a bin). The 4th
%                     column element is the index of the bin containing the
%                     previous state (from which the transition occured).
%                     The fifth column element is the index of the bin
%                     containing the new state (to which we transitioned).
%                     (i,4) == (i,5) iff (i,3) is true.
% mixed_transitions : j x 1 cell array (where j is the number of
%                     transitions, see above) where the i_th cell contains
%                     the first and second column elements of the i_th row
%                     of the matrix transitions (see above).


arguments
    X
    labels
    bin_width
    step_size
    nv.make_call = true
end

n_timepoints = length(labels);

% Bin data.
X_binned = bin_nd_data(X, 2, bin_width, step_size, ...
                       'excess_at', 'back', ...
                       'bin_excess', 'excess_only');
n_bins = size(X_binned, 2);

% Setup bin indices and then place this index on each original timepoint.
bin_indices = 1:n_bins;
orig_timepts_labeled_by_bin = repelem(bin_indices, bin_width);
orig_timepts_labeled_by_bin =  orig_timepts_labeled_by_bin(1:n_timepoints); % In case last bin was shorter

% Get bin starts and bin ends.
bin_starts = 1 : bin_width : n_timepoints;
bin_ends = bin_width : bin_width : n_timepoints;
if bin_ends(end) ~= n_timepoints, bin_ends = [bin_ends n_timepoints]; end

% Get indices (wrt to original timepoints) where a transition between
% labeled states occurs.
trans_labels_ind = find([0 diff(labels)]);
n_trans = length(trans_labels_ind);

% Preallocate n_trans x 5 array (see documentation for "transitions" under
% RETURNS).
transitions = nan(n_trans, 5);
mixed_transitions = cell(n_trans, 1);

if n_trans == 0
    % If there was only one label.
    bin_labels = repelem(labels(1), n_bins, 1);
else
    % Preallocate for bin labels. We will construct this as we go along
    % parsing transitions.
    bin_labels = nan(n_bins, 1);
    
    % Looping should be sufficently efficient/more readable.
    for i_trans = 1:n_trans
        % Get transition info ---------------------------------------------
        
        % Get index of timepoints before and after transition.
        i_from = trans_labels_ind(i_trans) - 1;
        i_to = trans_labels_ind(i_trans);
    
        % Store into transitions array.  
        transitions(i_trans, 1) = nan;
        transitions(i_trans, 1) = labels(i_from);
        transitions(i_trans, 2) = labels(i_to);
    
        % Check if timepoint before sits at the end of a bin.
        if ismember(i_from, bin_ends)
            assert(ismember(i_to, bin_starts));
            transitions(i_trans, 3) = 1;       
        else
            transitions(i_trans, 3) = 0;
        end
    
        % Get index of the bin containing the last timepoint for the state
        % from which we just transitioned. This covers both the case where
        % this final timepoint in the preceding state coincides with the
        % final timepoint in the previous bin, as well as the case where
        % the entire transition is contained by a single bin.
        transitions(i_trans, 4) = orig_timepts_labeled_by_bin(i_from);
    
        % Get index of timepoint containing the last timepoint for the
        % state from which we just transitioned.
        transitions(i_trans, 5) = i_from;
    
        % Construct labels vector for bins --------------------------------
        % The transition index attaches to the transition following our the
        % bins we want to label.
    
        % Get starting bin index for current epoch/state. Handle
        % leading/starting edge case.
        if i_trans == 1
            i_start = 1; 
        else 
            i_start = transitions(i_trans-1,4) + 1; 
        end
    
        % Set bin labels.
        i_end = transitions(i_trans, 4);  
        if transitions(i_trans, 3)
            bin_labels(i_start:i_end) = transitions(i_trans,1);
        elseif ~transitions(i_trans,3)
            bin_labels(i_start:i_end-1) = transitions(i_trans,1);
            bin_labels(i_end) = 777;
            mixed_transitions{i_trans} = transitions(i_trans, 1:2);
    
            % See documentation for 'make_call'.
            if nv.make_call               
                mixed_bin_timepoints ...
                    = labels(orig_timepts_labeled_by_bin == transitions(i_trans,4));                
                n_from = sum(mixed_bin_timepoints == transitions(i_trans,1));
                n_to = sum(mixed_bin_timepoints == transitions(i_trans,2));
                if n_from > n_to || n_from == n_to
                    bin_labels(i_end) = transitions(i_trans,1);
                elseif n_from < n_to
                    bin_labels(i_end) = transitions(i_trans,2);
                else
                    error("Unexpected results. n_from should be less than, " + ...
                          "equal to, or greater than n_to.")
                end
            end
        end
    end
    
    % Handle leading/trailing edge case of final epoch/state. Note that in
    % the above loop,  at the i_transition, we construct bin labels for the
    % preceding epoch/state. Thus, here we need to construct bin labels for
    % the final epoch/state following the final transition.
    bin_labels(i_end+1:end) = labels(end);
end


% Check that no bins contained three states (e.g., a given bin has state A,
% then state B, then state A again); throw error if so, mainly because this
% is problematic when trying to calculate transition stats, but also
% because it suggests that the bin size may have been chosen too large.
if any(diff(transitions(:,4)) == 0)
    e = MException('bin_with_transitions:multiple_transitions_within_bin', ...
                   'At least one bin contains 3 states. The bin size may be too large.');
    throw(e);
end

end
