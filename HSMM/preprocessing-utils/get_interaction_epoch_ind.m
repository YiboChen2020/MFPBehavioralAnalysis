function [intro_epochs, rmv_epochs] = get_interaction_epoch_ind(session_raw, second_subjects)
% For a given session, and a set of labels for the second subject (e.g.,
% {'M', 'F', 'Pup', 'Toy'} or some subset of these), identify the indices
% of the intro and removal of these second subjects.
%
% PARAMETERS
% ----------
% session_raw    : Scalar struct output (possibly within cell array) from
%                  extract_single_subj. Needed here are the fields
%                  .behavior (list of all behaviors annotated within
%                  session) and .Fstart and .Fstop (containing
%                  respectively, the start and stop timepoints of the
%                  annotated behavior).
% second_subject : Cell array whose elements contain the string label of
%                  a second subject. Must be {'M', 'F', 'Pup', 'Toy'} or
%                  some subset thereof. Note that these strings are
%                  appended to 'Intro_' or 'Rmv_' to generate the correct
%                  label, e.g., 'Intro_' + 'M' = 'Intro_M'.
%
% RETURNS
% -------
% intro_epochs : n_interaction_epochs x 1 vector whose i_th element is the
%                index (among annotations) of the i_th intro epoch.
% rmv_epochs   : n_interaction_epochs x 1 vector whose i_th element is the
%                index (among annotations) of the i_th remove epoch.
%
% Originally code from isolate_interactions; moved to its own function here
% on 9/15/22.

% Get indices of epochs for intro and removal of specified second subjects.
intro_epochs = [];
rmv_epochs = [];
for i_second_subj = 1:length(second_subjects)
    intro_epochs ...
        = [intro_epochs; ...
           ensure_is_col(find_epochs(session_raw.behaviors, ...
                                     sprintf('Intro_%s', ...
                                             second_subjects{i_second_subj}), ...
                                     'integer'))];
    rmv_epochs ...
        = [rmv_epochs; ...
           ensure_is_col(find_epochs(session_raw.behaviors, ...
                                     sprintf('Rmv_%s', ...
                                             second_subjects{i_second_subj}), ...
                                     'integer'))];
end

assert(length(intro_epochs) == length(rmv_epochs));

% Otherwise, will group interactions together by type and in order passed
% in. E.g., if second_subjects = {'M', 'F'}, and true order is 'M', 'F',
% 'M', we would get M1, M2, and F3, if we didn't sort. This edit was made
% on 5/16/23, but filenames for results from 4/20/23 were modified on
% 5/16/23.
intro_epochs = sort(intro_epochs);
rmv_epochs = sort(rmv_epochs);


end


% --------------------------------------------------
function epoch_ind = find_epochs(behaviors, str, ind_type)
% For an n x 1 (or 1 x n) cell array where each cell contains a string,
% return the indices of cells whose strings contain the specified pattern,
% str.

epoch_ind = ~cellfun(@isempty, regexp(behaviors, str));

if strcmp(ind_type, 'integer')
    epoch_ind = find(epoch_ind);
else
    assert(strcmp(ind_type, 'logical'), ...
           "ind_type must be 'integer' or 'logical'.")
end

end


% --------------------------------------------------
function v = ensure_is_col(v)
% Transpose a vector v iff it is a row vector.

if isrow(v), v = v'; end

end
