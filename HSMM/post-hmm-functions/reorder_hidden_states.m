function [hidden_states_, hidden_states_list_] ...
    = reorder_hidden_states(hidden_states, reorder_ind)
% Takes a hidden state trajectories and rearranges the hidden state labels
% according to the order specified in reorder_ind. XXX: see RETURNS below.
%
% PARAMETERS
% ----------
% hidden_states : k-vector whose p_th element is the hidden state label for
%                 the p_th timepoint.
% reorder_ind   : d-vector (for d as the number of hidden states) whose
%                 i_th element is the index of the hidden state (under the
%                 original ordering) that will be assigned to the i_th 
%
% RETURNS
% -------
% hidden_states_      : Hidden state trajectory after rearrangement. XXX:
%                       the hidden state labels in the rearranged
%                       trajectory will be elements of {1, ..., d}, so if
%                       set of original hidden state labels does not match
%                       this (e.g., if labels were 1, 3, 4, 7), this will
%                       cause the returned trajectory to have different
%                       hidden state labels than what was passed in. This
%                       was a quick and hacky way of doing things, need to
%                       check with others before changing to avoid creating
%                       compatability issues.
% hidden_states_list_ : Hidden states list after rearrangement. XXX: unlike
%                       with the trajectory label, the set of hidden state
%                       labels here matches the original state labels (they
%                       are reordered).
%
% Author: Jonathan Chien


% Reorder hidden state labels. First reorder list.
hidden_states_list = unique(hidden_states);
hidden_states_list_ = hidden_states_list(reorder_ind);

% Reorder trajectory. 
hidden_states_ = nan(size(hidden_states));
for i_state = 1:length(reorder_ind)
    % XXX: new hidden state labels will come from {1, ..., d}. 
    hidden_states_(hidden_states == hidden_states_list(i_state)) ...
        = find(hidden_states_list_ == hidden_states_list(i_state));
end

end
