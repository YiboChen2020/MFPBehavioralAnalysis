function all_behaviors = get_all_behaviors(path_to_data, subj_ind)
% Search all sessions from all subjects and get a list of all behaviors
% that occur.
%
% PARAMETERS
% ----------
% path_to_data : Path to directory containing data.
% ind_range    : n-vector of indices to try, where the i_th element is the
%                i_th subject index to try to load.
%
% RETURNS
% -------
% all_behaviors : b x 1 cell array, where b = the number of behaviors
%                 identified, and the i_th cell contains the name of the
%                 i_th behavior.
%
% Author: Jonathan Chien

addpath(path_to_data)

% Initialize list variable.
all_behaviors = [];

% Search across all subjects and sessions.
for i_subject = subj_ind
    % Load subject's data.
    try
        data = load(sprintf('MFP-ERa-GC%d.mat', i_subject));
    catch
        continue
    end

    % For current subject, determine number of sessions and get session
    % IDs.
    fnames = fieldnames(data.Raw);
    n_sessions = length(fnames) - 6; % There are six fields before the sessions
    all_session_ids = fnames(7:end);
    
    % For each session, search and get all behaviors.
    for i_session = 1:n_sessions
        all_behaviors = vertcat(all_behaviors, ...
                                unique(data.Raw.(all_session_ids{i_session}).behaviors));
    end
end

% Elminate redundancies across subjects/sessions.
all_behaviors = unique(all_behaviors);

end
