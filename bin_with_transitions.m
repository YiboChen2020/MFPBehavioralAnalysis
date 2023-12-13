function [X_binned, bin_labels, transitions, mixed_transitions] ...
    = bin_with_transitions(X, labels, bin_width, step_size, nv)

% TODO: get labels on bins, with special marker to mark bins containing a
% transition.


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

% Preallocate n_trans x 5 array where the (i,1) and (i,2) elements are
% respectively the state labels from and to which the i_th transition
% occured. The (i,3) element is a boolean, true if the timepoint preceding
% the transition coincides with a bin end (the state transition happens to
% align with a bin transition) and 0 otherwise (the state transition is
% wholly contained by a bin). The (i,4) element is the index of the
% bin containing the previous state (from which the transition occured).
% The (i,5) element is the index of the bin containing the new state (to
% which we transitioned). (i,4) == (i,5) iff (i,3) is true.
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
    
            % If a transition is contained within a bin, optionally assign
            % the whole bin to the state/behavior that comprises the
            % majority of timepoints in the bin (recall that there should
            % be only two states at most within a bin, as will be checked
            % below). In the event of a tie, assign the bin to the state
            % from which we transitioned.
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
    
    % Handle leading/trailing edge case of final epoch/state. Note that the
    % above loop operates such that at the i_transition, we construct bin
    % labels for the preceding epoch/state. Thus, here we need to construct
    % bin labels for the final epoch/state following the final transition.
    bin_labels(i_end+1:end) = labels(end);
end


% Check that if no bins contained three states (e.g., a given bin has state
% A, then state B, then state A again); throw error if so, mainly because
% this is annoying to deal with when trying to calculate transition stats,
% but also because it suggests that the bin size may have been chosen too
% large.
if any(diff(transitions(:,4)) == 0)
    e = MException('bin_with_transitions:multiple_transitions_within_bin', ...
                   'At least one bin contains 3 states. The bin size may be too large.');
    throw(e);
end

end
